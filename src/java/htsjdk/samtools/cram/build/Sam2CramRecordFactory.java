/*******************************************************************************
 * Copyright 2013 EMBL-EBI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package htsjdk.samtools.cram.build;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.cram.encoding.read_features.BaseQualityScore;
import htsjdk.samtools.cram.encoding.read_features.Deletion;
import htsjdk.samtools.cram.encoding.read_features.HardClip;
import htsjdk.samtools.cram.encoding.read_features.InsertBase;
import htsjdk.samtools.cram.encoding.read_features.Padding;
import htsjdk.samtools.cram.encoding.read_features.ReadFeature;
import htsjdk.samtools.cram.encoding.read_features.RefSkip;
import htsjdk.samtools.cram.encoding.read_features.SoftClip;
import htsjdk.samtools.cram.encoding.read_features.Substitution;
import htsjdk.samtools.cram.mask.RefMaskUtils;
import htsjdk.samtools.cram.structure.CramCompressionRecord;
import htsjdk.samtools.cram.structure.ReadTag;
import htsjdk.samtools.util.Log;

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class Sam2CramRecordFactory {

    public static final String UNKNOWN_READ_GROUP_ID = "UNKNOWN";
    public static final String UNKNOWN_READ_GROUP_SAMPLE = "UNKNOWN";

    private enum TREAT_TYPE {
        IGNORE, ALIGNMENT, INSERTION, SOFT_CLIP
    }

    /**
     * Reserved for later use.
     */
    private final TREAT_TYPE treatSoftClipsAs = TREAT_TYPE.SOFT_CLIP;

    public final static byte QS_asciiOffset = 33;
    public final static byte unsetQualityScore = 32;
    public final static byte ignorePositionsWithQualityScore = -1;

    private byte[] refBases;
    private byte[] refSNPs;
    private RefMaskUtils.RefMask refPile;

    public boolean captureUnmappedBases = true;
    public boolean captureUnmappedScores = false;

    private final ByteBuffer insertionBuf = ByteBuffer.allocate(1024);

    private static final Log log = Log.getInstance(Sam2CramRecordFactory.class);

    private final Map<String, Integer> readGroupMap = new HashMap<String, Integer>();

    private long landedRefMaskScores = 0;
    private long landedPiledScores = 0;
    private long landedTotalScores = 0;

    private boolean captureInsertScores = false;
    private boolean captureSubtitutionScores = false;
    private int uncategorisedQualityScoreCutoff = 0;
    public boolean captureAllTags = false;

    public boolean preserveReadNames = false;
    public Set<String> captureTags = new TreeSet<String>();
    public Set<String> ignoreTags = new TreeSet<String>();

    {
        ignoreTags.add(SAMTag.NM.name());
        ignoreTags.add(SAMTag.MD.name());
        ignoreTags.add(SAMTag.RG.name());
    }

    public boolean losslessQS = false;

    private final List<ReadTag> readTagList = new ArrayList<ReadTag>();

    private long baseCount = 0;
    private long featureCount = 0;

    public Sam2CramRecordFactory(final byte[] refBases, final SAMFileHeader samFileHeader) {
        this.refBases = refBases;

        final List<SAMReadGroupRecord> readGroups = samFileHeader.getReadGroups();
        for (int i = 0; i < readGroups.size(); i++) {
            final SAMReadGroupRecord readGroupRecord = readGroups.get(i);
            readGroupMap.put(readGroupRecord.getId(), i);
        }
    }

    public CramCompressionRecord createCramRecord(final SAMRecord record) {
        final CramCompressionRecord cramRecord = new CramCompressionRecord();
        cramRecord.setForcePreserveQualityScores(captureUnmappedScores);
        if (record.getReadPairedFlag()) {
            cramRecord.mateAlignmentStart = record.getMateAlignmentStart();
            cramRecord.setMateUmapped(record.getMateUnmappedFlag());
            cramRecord.setMateNegativeStrand(record.getMateNegativeStrandFlag());
            cramRecord.mateSequnceID = record.getMateReferenceIndex();
            cramRecord.setDetached(true);
        } else cramRecord.mateSequnceID = -1;
        cramRecord.sequenceId = record.getReferenceIndex();
        cramRecord.readName = record.getReadName();
        cramRecord.alignmentStart = record.getAlignmentStart();

        cramRecord.setMultiFragment(record.getReadPairedFlag());
        cramRecord.setProperPair(record.getReadPairedFlag()
                && record.getProperPairFlag());
        cramRecord.setSegmentUnmapped(record.getReadUnmappedFlag());
        cramRecord.setNegativeStrand(record.getReadNegativeStrandFlag());
        cramRecord.setFirstSegment(record.getReadPairedFlag()
                && record.getFirstOfPairFlag());
        cramRecord.setLastSegment(record.getReadPairedFlag()
                && record.getSecondOfPairFlag());
        cramRecord.setSecondaryAlignment(record.getNotPrimaryAlignmentFlag());
        cramRecord.setVendorFiltered(record.getReadFailsVendorQualityCheckFlag());
        cramRecord.setDuplicate(record.getDuplicateReadFlag());

        cramRecord.readLength = record.getReadLength();
        cramRecord.mappingQuality = record.getMappingQuality();
        cramRecord.setDuplicate(record.getDuplicateReadFlag());

        cramRecord.templateSize = record.getInferredInsertSize();

        final SAMReadGroupRecord readGroup = record.getReadGroup();
        if (readGroup != null)
            cramRecord.readGroupID = readGroupMap.get(readGroup.getId());
        else
            cramRecord.readGroupID = -1;

        if (!record.getReadPairedFlag())
            cramRecord.setLastSegment(false);
        else {
            if (record.getFirstOfPairFlag())
                cramRecord.setLastSegment(false);
            else if (record.getSecondOfPairFlag())
                cramRecord.setLastSegment(true);
            else
                cramRecord.setLastSegment(true);
        }

        if (!cramRecord.isSegmentUnmapped()) {
            cramRecord.readFeatures = checkedCreateVariations(cramRecord,
                    record);
        }

        cramRecord.readBases = record.getReadBases();
        cramRecord.qualityScores = record.getBaseQualities();
        landedTotalScores += cramRecord.readLength;

        readTagList.clear();
        if (captureAllTags) {
            final List<SAMTagAndValue> attributes = record.getAttributes();
            for (final SAMTagAndValue tv : attributes) {
                if (ignoreTags.contains(tv.tag))
                    continue;
                readTagList.add(ReadTag.deriveTypeFromValue(tv.tag, tv.value));
            }
        } else {
            if (!captureTags.isEmpty()) {
                final List<SAMTagAndValue> attributes = record.getAttributes();
                cramRecord.tags = new ReadTag[attributes.size()];

                for (final SAMTagAndValue tv : attributes) {
                    if (captureTags.contains(tv.tag)) {
                        readTagList.add(ReadTag.deriveTypeFromValue(tv.tag, tv.value));
                    }
                }
            }
        }
        cramRecord.tags = readTagList
                .toArray(new ReadTag[readTagList.size()]);

        cramRecord.setVendorFiltered(record.getReadFailsVendorQualityCheckFlag());

        if (preserveReadNames)
            cramRecord.readName = record.getReadName();

        return cramRecord;
    }

    /**
     * A wrapper method to provide better diagnostics for
     * ArrayIndexOutOfBoundsException.
     *
     * @param cramRecord
     * @param samRecord
     * @return
     */
    private List<ReadFeature> checkedCreateVariations(final CramCompressionRecord cramRecord,
                                                      final SAMRecord samRecord) {
        try {
            return createVariations(cramRecord, samRecord);
        } catch (final ArrayIndexOutOfBoundsException e) {
            log.error("Reference bases array length=" + refBases.length);
            log.error("Offensive CRAM record: " + cramRecord.toString());
            log.error("Offensive SAM record: " + samRecord.getSAMString());
            throw e;
        }
    }

    private List<ReadFeature> createVariations(final CramCompressionRecord cramRecord,
                                               final SAMRecord samRecord) {
        final List<ReadFeature> features = new LinkedList<ReadFeature>();
        int zeroBasedPositionInRead = 0;
        int alignmentStartOffset = 0;
        int cigarElementLength;

        final List<CigarElement> cigarElements = samRecord.getCigar()
                .getCigarElements();

        final byte[] bases = samRecord.getReadBases();
        final byte[] qualityScore = samRecord.getBaseQualities();

        for (final CigarElement cigarElement : cigarElements) {
            cigarElementLength = cigarElement.getLength();
            final CigarOperator operator = cigarElement.getOperator();

            switch (operator) {
                case D:
                    features.add(new Deletion(zeroBasedPositionInRead + 1,
                            cigarElementLength));
                    break;
                case N:
                    features.add(new RefSkip(zeroBasedPositionInRead + 1,
                            cigarElementLength));
                    break;
                case P:
                    features.add(new Padding(zeroBasedPositionInRead + 1,
                            cigarElementLength));
                    break;
                case H:
                    features.add(new HardClip(zeroBasedPositionInRead + 1,
                            cigarElementLength));
                    break;
                case S:
                    addSoftClip(features, zeroBasedPositionInRead,
                            cigarElementLength, bases, qualityScore);
                    break;
                case I:
                    addInsertion(features, zeroBasedPositionInRead,
                            cigarElementLength, bases, qualityScore);
                    break;
                case M:
                case X:
                case EQ:
                    addSubstitutionsAndMaskedBases(cramRecord, features,
                            zeroBasedPositionInRead, alignmentStartOffset,
                            cigarElementLength, bases, qualityScore);
                    break;
                default:
                    throw new IllegalArgumentException(
                            "Unsupported cigar operator: "
                                    + cigarElement.getOperator());
            }

            if (cigarElement.getOperator().consumesReadBases())
                zeroBasedPositionInRead += cigarElementLength;
            if (cigarElement.getOperator().consumesReferenceBases())
                alignmentStartOffset += cigarElementLength;
        }

        this.baseCount += bases.length;
        this.featureCount += features.size();

        return features;
    }

    private void addSoftClip(final List<ReadFeature> features,
                             final int zeroBasedPositionInRead, final int cigarElementLength, final byte[] bases,
                             final byte[] scores) {
        final byte[] insertedBases = Arrays.copyOfRange(bases,
                zeroBasedPositionInRead, zeroBasedPositionInRead
                        + cigarElementLength);

        final SoftClip v = new SoftClip(zeroBasedPositionInRead + 1, insertedBases);
        features.add(v);
    }

    private void addHardClip(final List<ReadFeature> features,
                             final int zeroBasedPositionInRead, final int cigarElementLength, final byte[] bases,
                             final byte[] scores) {
        final byte[] insertedBases = Arrays.copyOfRange(bases,
                zeroBasedPositionInRead, zeroBasedPositionInRead
                        + cigarElementLength);

        final HardClip v = new HardClip(zeroBasedPositionInRead + 1, insertedBases.length);
        features.add(v);
    }

    private void addInsertion(final List<ReadFeature> features,
                              final int zeroBasedPositionInRead, final int cigarElementLength, final byte[] bases,
                              final byte[] scores) {
        final byte[] insertedBases = Arrays.copyOfRange(bases,
                zeroBasedPositionInRead, zeroBasedPositionInRead
                        + cigarElementLength);

        for (int i = 0; i < insertedBases.length; i++) {
            // single base insertion:
            final InsertBase ib = new InsertBase();
            ib.setPosition(zeroBasedPositionInRead + 1 + i);
            ib.setBase(insertedBases[i]);
            features.add(ib);
            if (losslessQS || scores == null || scores.length < bases.length)
                continue;
            final boolean qualityMasked = (scores[i] < uncategorisedQualityScoreCutoff);
            if (captureInsertScores || qualityMasked) {
                final byte score = (byte) (QS_asciiOffset + scores[zeroBasedPositionInRead
                        + i]);
                // if (score >= QS_asciiOffset) {
                features.add(new BaseQualityScore(zeroBasedPositionInRead + 1
                        + i, score));
                landedTotalScores++;
                // }
            }
        }
    }

    private void addSubstitutionsAndMaskedBases(final CramCompressionRecord cramRecord,
                                                final List<ReadFeature> features, final int fromPosInRead,
                                                final int alignmentStartOffset, final int nofReadBases, final byte[] bases,
                                                final byte[] qualityScore) {
        int oneBasedPositionInRead;
        final boolean noQS = (qualityScore.length == 0);

        int i = 0;
        boolean qualityAdded;
        boolean qualityMasked;
        byte refBase;
        for (i = 0; i < nofReadBases; i++) {
            oneBasedPositionInRead = i + fromPosInRead + 1;
            final int refCoord = (cramRecord.alignmentStart + i + alignmentStartOffset) - 1;
            qualityAdded = false;
            if (refCoord >= refBases.length)
                refBase = 'N';
            else
                refBase = refBases[refCoord];
            refBase = Utils.normalizeBase(refBase);

            if (bases[i + fromPosInRead] != refBase) {
                final Substitution sv = new Substitution();
                sv.setPosition(oneBasedPositionInRead);
                sv.setBase(bases[i + fromPosInRead]);
                sv.setRefernceBase(refBase);
                sv.setBaseChange(null);

                features.add(sv);

                if (losslessQS || noQS)
                    continue;

                if (captureSubtitutionScores) {
                    final byte score = (byte) (QS_asciiOffset + qualityScore[i
                            + fromPosInRead]);
                    features.add(new BaseQualityScore(oneBasedPositionInRead,
                            score));
                    qualityAdded = true;
                }
            }

            if (noQS)
                continue;

            if (!qualityAdded && refSNPs != null) {
                final byte snpOrNot = refSNPs[refCoord];
                if (snpOrNot != 0) {
                    final byte score = (byte) (QS_asciiOffset + qualityScore[i
                            + fromPosInRead]);
                    features.add(new BaseQualityScore(oneBasedPositionInRead,
                            score));
                    qualityAdded = true;
                    landedRefMaskScores++;
                }
            }

            if (!qualityAdded && refPile != null) {
                if (refPile.shouldStore(refCoord, refBase)) {
                    final byte score = (byte) (QS_asciiOffset + qualityScore[i
                            + fromPosInRead]);
                    features.add(new BaseQualityScore(oneBasedPositionInRead,
                            score));
                    qualityAdded = true;
                    landedPiledScores++;
                }
            }

            qualityMasked = (qualityScore[i + fromPosInRead] < uncategorisedQualityScoreCutoff);
            if (!qualityAdded && qualityMasked) {
                final byte score = (byte) (QS_asciiOffset + qualityScore[i
                        + fromPosInRead]);
                features.add(new BaseQualityScore(oneBasedPositionInRead, score));
                qualityAdded = true;
            }

            if (qualityAdded)
                landedTotalScores++;
        }
    }

    public boolean isCaptureInsertScores() {
        return captureInsertScores;
    }

    public void setCaptureInsertScores(final boolean captureInsertScores) {
        this.captureInsertScores = captureInsertScores;
    }

    public boolean isCaptureSubtitutionScores() {
        return captureSubtitutionScores;
    }

    public void setCaptureSubtitutionScores(final boolean captureSubtitutionScores) {
        this.captureSubtitutionScores = captureSubtitutionScores;
    }

    public int getUncategorisedQualityScoreCutoff() {
        return uncategorisedQualityScoreCutoff;
    }

    public void setUncategorisedQualityScoreCutoff(
            final int uncategorisedQualityScoreCutoff) {
        this.uncategorisedQualityScoreCutoff = uncategorisedQualityScoreCutoff;
    }

    public long getLandedRefMaskScores() {
        return landedRefMaskScores;
    }

    public long getLandedPiledScores() {
        return landedPiledScores;
    }

    public long getLandedTotalScores() {
        return landedTotalScores;
    }

    public boolean isCaptureUnmappedBases() {
        return captureUnmappedBases;
    }

    public void setCaptureUnmappedBases(final boolean captureUnmappedBases) {
        this.captureUnmappedBases = captureUnmappedBases;
    }

    public boolean isCaptureUnmappedScores() {
        return captureUnmappedScores;
    }

    public void setCaptureUnmappedScores(final boolean captureUnmappedScores) {
        this.captureUnmappedScores = captureUnmappedScores;
    }

    public byte[] getRefBases() {
        return refBases;
    }

    public void setRefBases(final byte[] refBases) {
        this.refBases = refBases;
    }

    public byte[] getRefSNPs() {
        return refSNPs;
    }

    public void setRefSNPs(final byte[] refSNPs) {
        this.refSNPs = refSNPs;
    }

    public RefMaskUtils.RefMask getRefPile() {
        return refPile;
    }

    public Map<String, Integer> getReadGroupMap() {
        return readGroupMap;
    }

    public void setRefPile(final RefMaskUtils.RefMask refPile) {
        this.refPile = refPile;
    }


    public long getBaseCount() {
        return baseCount;
    }

    public long getFeatureCount() {
        return featureCount;
    }

    public boolean isPreserveReadNames() {
        return preserveReadNames;
    }

    public void setPreserveReadNames(final boolean preserveReadNames) {
        this.preserveReadNames = preserveReadNames;
    }

}
