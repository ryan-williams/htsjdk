package htsjdk.samtools;

import htsjdk.samtools.util.StringUtil;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.List;

/**
 * BAMRecord which is created with just a raw byte[] array representing both the fixed and
 * variable-length sections of a BAM record.
 */
public class FastBAMRecord extends AbstractReadRecord {

    private final SAMFileHeader header;
    private byte[] record = null;

    private static final int referenceIdOffset = 0;
    private static final int coordinateOffset = referenceIdOffset + 4;
    private static final int readNameLengthOffset = coordinateOffset + 4;
    private static final int mappingQualityOffset = readNameLengthOffset + 1;
    private static final int binOffset = mappingQualityOffset + 1;
    private static final int cigarLenOffset = binOffset + 2;
    private static final int flagsOffset = cigarLenOffset + 2;
    private static final int readLenOffset = flagsOffset + 2;
    private static final int mateReferenceIdOffset = readLenOffset + 4;
    private static final int mateCoordinateOffset = mateReferenceIdOffset + 4;
    private static final int insertSizeOffset = mateCoordinateOffset + 4;
    private static final int readNameOffset = insertSizeOffset + 4;

    private static final byte[] NULL_SEQUENCE = new byte[0];
    private static final String NULL_SEQUENCE_STRING = "*";

    public FastBAMRecord(final SAMFileHeader header, final byte[] record) {
        this.header = header;
        this.record = record;

        getFlags(); // decode flags
    }


    private static final int decodeUShort(final byte[] buffer, final int offset) {
        final byte b1 = buffer[offset];
        final byte b2 = buffer[offset+1];
        return ((b2 & 0xFF) << 8) | (b1 & 0xFF);
    }

    private static final int decodeInt(final byte[] buffer, final int offset) {
        final byte b1 = buffer[offset];
        final byte b2 = buffer[offset+1];
        final byte b3 = buffer[offset+2];
        final byte b4 = buffer[offset+3];
        return ((b4 & 0xFF) << 24) | ((b3 & 0xFF) << 16) | ((b2 & 0xFF) << 8) | (b1 & 0xFF);
    }

    /*
    final int referenceID = this.binaryCodec.readInt();
    final int coordinate = this.binaryCodec.readInt() + 1;
    final short readNameLength = this.binaryCodec.readUByte();
    final short mappingQuality = this.binaryCodec.readUByte();
    final int bin = this.binaryCodec.readUShort();
    final int cigarLen = this.binaryCodec.readUShort();
    final int flags = this.binaryCodec.readUShort();
    final int readLen = this.binaryCodec.readInt();
    final int mateReferenceID = this.binaryCodec.readInt();
    final int mateCoordinate = this.binaryCodec.readInt() + 1;
    final int insertSize = this.binaryCodec.readInt();
     */

    /*
        getReadName();
        getCigar();
        getReadBases();
        getBaseQualities();
        getBinaryAttributes();
     */

    private String decodeReadName() {
        // Don't include terminating null
        return StringUtil.bytesToString(record, readNameOffset, getReadNameLength() - 1);
    }

    private String decodeReadBases() {
        final int readLength = getReadLength();
        if (readLength == 0) {
            return NULL_SEQUENCE_STRING;
        }
        final int basesOffset = readNameOffset + getReadNameLength() + cigarBufferSize();
        //SAMUtils.compressedBasesToBytes(readLength, record, basesOffset)
        return StringUtil.bytesToString(record, basesOffset, basesBufferSize());
    }


    private String readName = null;
    @Override
    public String getReadName() {
        if(readName == null) {
            readName = decodeReadName();
        }
        return readName;
    }

    @Override
    public int getReadNameLength() {
        return record[readNameLengthOffset];
    }



    private int cigarLength = -1;
    @Override
    public int getCigarLength() {
        if(cigarLength == -1) {
            cigarLength = decodeUShort(record, cigarLenOffset);
        }
        return cigarLength;
    }


    private String readString;
    @Override
    public String getReadString() {
        if(readString == null) {
            readString = decodeReadBases();
        }
        return readString;
    }




    private int readLength = -1;
    @Override
    public int getReadLength() {
        if(readLength == -1) {
            readLength = decodeInt(record, readLenOffset);
        }
        return readLength;
    }

    @Override
    public int getMappingQuality() {
        return record[mappingQualityOffset];
    }


    private int alignmentStart = -1;
    @Override
    public int getAlignmentStart() {
        if(alignmentStart == -1) {
            alignmentStart = decodeInt(record, coordinateOffset);
        }

        return alignmentStart;
    }


    private Cigar cigar = null;
    @Override
    public Cigar getCigar() {
        if(cigar == null) {
            /*
            System.out.println("Record length: " + getReadName());
            System.out.println("Record length: " + record.length);
            System.out.println("Read Name Length: " + getReadNameLength());
            System.out.println("Start: " + (readNameOffset + getReadNameLength()));
            System.out.println("Cigar Length: " + getCigarLength());
            */
            final ByteBuffer byteBuffer  = ByteBuffer.wrap(record, readNameOffset + getReadNameLength(), cigarBufferSize());
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            cigar = BinaryCigarCodec.getSingleton().decode(byteBuffer);

        }
        return cigar;
    }

    private List<AlignmentBlock> alignmentBlocks = null;
    @Override
    public List<AlignmentBlock> getAlignmentBlocks() {
        if (alignmentBlocks == null) {
            alignmentBlocks = SAMUtils.getAlignmentBlocks(getCigar(), getAlignmentStart(), "read cigar");
        }
        return alignmentBlocks;
    }



    private int cigarBufferSize() {
        return getCigarLength() * 4;
    }

    private int basesBufferSize() {
        return (readLength + 1)/2;
    }



    @Override
    public byte[] getReadBases() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getBaseQualityString() {
        throw new RuntimeException("Not implemented");
    }



    private byte[] baseQualities = null;
    @Override
    public byte[] getBaseQualities() {
        //mHeader.getSequence(referenceIndex).getSequenceName();
        if(baseQualities == null) {
            baseQualities = new byte[getReadLength()];
            final int baseQualityOffset = readNameOffset + getReadNameLength() + cigarBufferSize() + basesBufferSize();

            System.arraycopy(record, baseQualityOffset, baseQualities, 0, baseQualities.length);
        }
        return baseQualities;
    }



    @Override
    public byte[] getOriginalBaseQualities() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getReferenceName() {
        throw new RuntimeException("Not implemented");
    }


    private int referenceIndex = NO_ALIGNMENT_REFERENCE_INDEX;
    @Override
    public Integer getReferenceIndex() {
        if(referenceIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
            referenceIndex = decodeInt(record, referenceIdOffset);
        }
        return referenceIndex;
    }


    @Override
    public String getMateReferenceName() {
        throw new RuntimeException("Not implemented");
    }


    @Override
    public Integer getMateReferenceIndex() {
        throw new RuntimeException("Not implemented");
    }




    @Override
    public int getAlignmentEnd() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getUnclippedStart() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getUnclippedEnd() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getReferencePositionAtReadPosition(final int offset) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getMateAlignmentStart() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getInferredInsertSize() {
        throw new RuntimeException("Not implemented");
    }


    @Override
    public String getCigarString() {
        throw new RuntimeException("Not implemented");
    }






    @Override
    public SAMReadGroupRecord getReadGroup() {
        throw new RuntimeException("Not implemented");
    }

    private int flags = -1;
    @Override
    public int getFlags() {
        if(flags == -1) {
            flags = decodeUShort(record, flagsOffset);
        }
        return flags;
    }

    @Override
    public void setFlags(final int value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public boolean getReadPairedFlag() {
        return (flags & READ_PAIRED_FLAG) != 0;
    }

    @Override
    public boolean getProperPairFlag() {
        return (flags & PROPER_PAIR_FLAG) != 0;
    }

    @Override
    public boolean getReadUnmappedFlag() {
        return (flags & READ_UNMAPPED_FLAG) != 0;
    }

    @Override
    public boolean getMateUnmappedFlag() {
        return (flags & MATE_UNMAPPED_FLAG) != 0;
    }

    @Override
    public boolean getReadNegativeStrandFlag() {
        return (flags & READ_STRAND_FLAG) != 0;
    }

    @Override
    public boolean getMateNegativeStrandFlag() {
        return (flags & MATE_STRAND_FLAG) != 0;
    }

    @Override
    public boolean getFirstOfPairFlag() {
        return (flags & FIRST_OF_PAIR_FLAG) != 0;
    }

    @Override
    public boolean getSecondOfPairFlag() {
        return (flags & SECOND_OF_PAIR_FLAG) != 0;
    }

    @Override
    public boolean getNotPrimaryAlignmentFlag() {
        return (flags & NOT_PRIMARY_ALIGNMENT_FLAG) != 0;
    }

    @Override
    public boolean getSupplementaryAlignmentFlag() {
        return (flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    @Override
    public boolean getReadFailsVendorQualityCheckFlag() {
        return (flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    @Override
    public boolean getDuplicateReadFlag() {
        return (flags & DUPLICATE_READ_FLAG) != 0;
    }

    @Override
    public void setReadPairedFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setProperPairFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadUmappedFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadUnmappedFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setMateUnmappedFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadNegativeStrandFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setMateNegativeStrandFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setFirstOfPairFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setSecondOfPairFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setNotPrimaryAlignmentFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setSupplementaryAlignmentFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadFailsVendorQualityCheckFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setDuplicateReadFlag(final boolean flag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public boolean isSecondaryOrSupplementary() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public ValidationStringency getValidationStringency() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setValidationStringency(final ValidationStringency validationStringency) {
        //throw new RuntimeException("Not implemented");
        //System.err.println("WARNING: setValidationStringency(" + validationStringency + ") ignored.");
    }

    @Override
    public Object getAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Integer getIntegerAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Short getShortAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Byte getByteAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getStringAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Character getCharacterAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Float getFloatAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public byte[] getByteArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public byte[] getUnsignedByteArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public byte[] getSignedByteArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public short[] getUnsignedShortArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public short[] getSignedShortArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int[] getUnsignedIntArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int[] getSignedIntArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public float[] getFloatArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public boolean isUnsignedArrayAttribute(final String tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setAttribute(final short tag, final Object value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Object getAttribute(final short tag) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setAttribute(final String tag, final Object value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setUnsignedArrayAttribute(final String tag, final Object value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void clearAttributes() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public List<SAMTagAndValue> getAttributes() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public SAMBinaryTagAndValue getBinaryAttributes() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Integer getIndexingBin() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setIndexingBin(final Integer mIndexingBin) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int computeIndexingBin() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public SAMFileHeader getHeader() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setHeader(final SAMFileHeader header) {
        throw new RuntimeException("Not implemented");
    }


    @Override
    public String getSAMString() {
        throw new RuntimeException("Not implemented");
    }












    @Override
    public SAMFileSource getFileSource() {
        throw new RuntimeException("Not implemented");
    }






    @Override
    public byte[] getVariableBinaryRepresentation() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getAttributesBinarySize() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String format() {
        throw new RuntimeException("Not implemented");
    }



    @Override
    public List<SAMValidationError> validateCigar(final long recordNumber) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public List<SAMValidationError> isValid() {
        throw new RuntimeException("Not implemented");
    }





    @Override
    public void setFileSource(final SAMFileSource fileSource) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void eagerDecode() {
        throw new RuntimeException("Not implemented");
    }



    @Override
    public void setReadName(final String value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadString(final String value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReadBases(final byte[] value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setBaseQualityString(final String value) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setBaseQualities(final byte[] value) {
        throw new RuntimeException("Not implemented");
    }


    @Override
    public void setOriginalBaseQualities(final byte[] oq) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setReferenceName(final String value) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setReferenceIndex(final int referenceIndex) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setMateReferenceName(final String mateReferenceName) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setMateReferenceIndex(final int referenceIndex) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setAlignmentStart(final int value) {
        throw new RuntimeException("Not implemented");
    }
    @Override
    public void setAlignmentEnd(final int value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setMateAlignmentStart(final int mateAlignmentStart) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setInferredInsertSize(final int inferredInsertSize) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setMappingQuality(final int value) {
        throw new RuntimeException("Not implemented");
    }


    @Override
    public void setCigarString(final String value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setCigar(final Cigar cigar) {
        throw new RuntimeException("Not implemented");
    }
}
