package htsjdk.samtools2;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.ValidationStringency;

import java.util.List;

/**
 * An IntelliJ extraction of the SAMRecord interface which can now gradually be morphed into the interface
 * described in Nils' ReadRecord2.java.
 */
public interface ReadRecord extends Cloneable {
    /**
     * Alignment score for a good alignment, but where computing a Phred-score is not feasible.
     */
    int UNKNOWN_MAPPING_QUALITY = 255;
    /**
     * Alignment score for an unaligned read.
     */
    int NO_MAPPING_QUALITY = 0;
    /**
     * If a read has this reference name, it is unaligned, but not all unaligned reads have
     * this reference name (see above).
     */
    String NO_ALIGNMENT_REFERENCE_NAME = "*";
    /**
     * If a read has this reference index, it is unaligned, but not all unaligned reads have
     * this reference index (see above).
     */
    int NO_ALIGNMENT_REFERENCE_INDEX = -1;
    /**
     * Cigar string for an unaligned read.
     */
    String NO_ALIGNMENT_CIGAR = "*";
    /**
     * If a read has reference name "*", it will have this value for position.
     */
    int NO_ALIGNMENT_START = GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
    /**
     * This should rarely be used, since a read with no sequence doesn't make much sense.
     */
    byte[] NULL_SEQUENCE = new byte[0];
    String NULL_SEQUENCE_STRING = "*";
    /**
     * This should rarely be used, since all reads should have quality scores.
     */
    byte[] NULL_QUALS = new byte[0];
    String NULL_QUALS_STRING = "*";
    /**
     * abs(insertSize) must be <= this
     */
    int MAX_INSERT_SIZE = 1<<29;

    String getReadName();

    int getReadNameLength();

    void setReadName(String value);

    String getReadString();

    void setReadString(String value);

    byte[] getReadBases();

    void setReadBases(byte[] value);

    int getReadLength();

    String getBaseQualityString();

    void setBaseQualityString(String value);

    byte[] getBaseQualities();

    void setBaseQualities(byte[] value);

    byte[] getOriginalBaseQualities();

    void setOriginalBaseQualities(byte[] oq);

    String getReferenceName();

    void setReferenceName(String value);

    Integer getReferenceIndex();

    void setReferenceIndex(int referenceIndex);

    String getMateReferenceName();

    void setMateReferenceName(String mateReferenceName);

    Integer getMateReferenceIndex();

    void setMateReferenceIndex(int referenceIndex);

    int getAlignmentStart();

    void setAlignmentStart(int value);

    int getAlignmentEnd();

    int getUnclippedStart();

    int getUnclippedEnd();

    int getReferencePositionAtReadPosition(int offset);

    void setAlignmentEnd(int value);

    int getMateAlignmentStart();

    void setMateAlignmentStart(int mateAlignmentStart);

    int getInferredInsertSize();

    void setInferredInsertSize(int inferredInsertSize);

    int getMappingQuality();

    void setMappingQuality(int value);

    String getCigarString();

    void setCigarString(String value);

    Cigar getCigar();

    int getCigarLength();

    void setCigar(Cigar cigar);

    SAMReadGroupRecord getReadGroup();

    int getFlags();

    void setFlags(int value);

    boolean isPaired();

    boolean getProperPairFlag();

    boolean getReadUnmappedFlag();

    boolean getMateUnmappedFlag();

    boolean getReadNegativeStrandFlag();

    boolean getMateNegativeStrandFlag();

    boolean getFirstOfPairFlag();

    boolean getSecondOfPairFlag();

    boolean getNotPrimaryAlignmentFlag();

    boolean getSupplementaryAlignmentFlag();

    boolean getReadFailsVendorQualityCheckFlag();

    boolean getDuplicateReadFlag();

    void setReadPairedFlag(boolean flag);

    void setProperPairFlag(boolean flag);

    void setReadUmappedFlag(boolean flag);

    void setReadUnmappedFlag(boolean flag);

    void setMateUnmappedFlag(boolean flag);

    void setReadNegativeStrandFlag(boolean flag);

    void setMateNegativeStrandFlag(boolean flag);

    void setFirstOfPairFlag(boolean flag);

    void setSecondOfPairFlag(boolean flag);

    void setNotPrimaryAlignmentFlag(boolean flag);

    void setSupplementaryAlignmentFlag(boolean flag);

    void setReadFailsVendorQualityCheckFlag(boolean flag);

    void setDuplicateReadFlag(boolean flag);

    boolean isSecondaryOrSupplementary();

    ValidationStringency getValidationStringency();

    void setValidationStringency(ValidationStringency validationStringency);

    Object getAttribute(String tag);

    Integer getIntegerAttribute(String tag);

    Short getShortAttribute(String tag);

    Byte getByteAttribute(String tag);

    String getStringAttribute(String tag);

    Character getCharacterAttribute(String tag);

    Float getFloatAttribute(String tag);

    byte[] getByteArrayAttribute(String tag);

    byte[] getUnsignedByteArrayAttribute(String tag);

    byte[] getSignedByteArrayAttribute(String tag);

    short[] getUnsignedShortArrayAttribute(String tag);

    short[] getSignedShortArrayAttribute(String tag);

    int[] getUnsignedIntArrayAttribute(String tag);

    int[] getSignedIntArrayAttribute(String tag);

    float[] getFloatArrayAttribute(String tag);

    boolean isUnsignedArrayAttribute(String tag);

    Object getAttribute(short tag);

    void setAttribute(String tag, Object value);

    void setUnsignedArrayAttribute(String tag, Object value);

    void clearAttributes();

    List<SAMTagAndValue> getAttributes();

    SAMFileHeader getHeader();

    void setHeader(SAMFileHeader header);

    byte[] getVariableBinaryRepresentation();

    int getAttributesBinarySize();

    String format();

    List<AlignmentBlock> getAlignmentBlocks();

    List<SAMValidationError> validateCigar(long recordNumber);

    @Override
    boolean equals(Object o);

    @Override
    int hashCode();

    List<SAMValidationError> isValid();

    SAMFileSource getFileSource();

    public Object clone() throws CloneNotSupportedException;

    @Override
    String toString();

    String getSAMString();

    /**
     * Replace any existing attributes with the given linked item.
     */
    void setAttributes(final SAMBinaryTagAndValue attributes);

    /**
     * @return Pointer to the first of the tags.  Returns null if there are no tags.
     */
    SAMBinaryTagAndValue getBinaryAttributes();

    /**
     * Tag name and value of an attribute, for getAttributes() method.
     */
    public static class SAMTagAndValue {
        public final String tag;
        public final Object value;

        public SAMTagAndValue(final String tag, final Object value) {
            this.tag = tag;
            this.value = value;
        }
    }
}
