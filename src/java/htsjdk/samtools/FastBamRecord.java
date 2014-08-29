package htsjdk.samtools;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.NoSuchElementException;


/**
 * BAMRecord that keeps its data in a single raw byte[] array which
 * encodes both the fixed and variable-length sections of the BAM record
 * as they would be in a BAM file on disk.
 *
 * TODO how does this interact with BAM index?
 * TODO replace char attribute type with a type enum for the known attributes?
 */
public class FastBAMRecord extends AbstractReadRecord
{

    //==================================================================
    //=========================== STATIC OFFSETS  ======================
    //==================================================================

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
    private static final int startOfVariableRegion = insertSizeOffset + 4;
    static {
        assert startOfVariableRegion == BAMFileConstants.FIXED_BLOCK_SIZE;
    }

    //==================================================================
    //=========================== ENCODE / DECODE ======================
    //==================================================================

    private static final short decodeShort(final byte[] buffer, final int offset) {
        final byte b1 = buffer[offset];
        final byte b2 = buffer[offset+1];
        return (short) (((b2 & 0xFF) << 8) | (b1 & 0xFF));
    }

    private static final int decodeUnsignedShort(final byte[] buffer, final int offset) {
        return decodeShort(buffer, offset) & 0xFFFF;
    }

    private static final int decodeInt(final byte[] buffer, final int offset) {
        final byte b1 = buffer[offset];
        final byte b2 = buffer[offset+1];
        final byte b3 = buffer[offset+2];
        final byte b4 = buffer[offset+3];
        return ((b4 & 0xFF) << 24) | ((b3 & 0xFF) << 16) | ((b2 & 0xFF) << 8) | (b1 & 0xFF);
    }

    private static final long decodeUnsignedInt(final byte[] buffer, final int offset) {
        return decodeInt(buffer, offset) & 0xFFFFFFFFL;
    }

    private static final float decodeFloat(final byte[] buffer, final int offset) {
        throw new IllegalStateException("Not implemented");
        //return Float.intBitsToFloat(decodeInt(buffer, offset));
    }

    private static final void encodeShort(final byte[] buffer, final int offset, short value) {
        buffer[offset] = (byte) value;
        buffer[offset+1] = (byte) (value >> 8);
    }

    private static final void encodeInt(final byte[] buffer, final int offset, int value) {
        buffer[offset] = (byte) value;
        buffer[offset+1] = (byte) (value >> 8);
        buffer[offset+2] = (byte) (value >> 16);
        buffer[offset+3] = (byte) (value >> 24);
    }

    //==================================================================
    //=========================== INSTANCE VARS   ======================
    //==================================================================

    //the underlying byte array containing all data for this record
    private byte[] record = null;

    //caches for decoded values
    protected String mReadName = null;
    protected byte[] mReadBases = null;
    protected byte[] mBaseQualities = null;

    private static final int UNINITIALIZED = -1;
    private int n_cigar_ops = UNINITIALIZED;
    private int mReadLength = UNINITIALIZED;
    private int referenceIndex = UNINITIALIZED;
    protected int mAlignmentStart = UNINITIALIZED;
    protected int mAlignmentEnd = UNINITIALIZED;

    protected String referenceName = null;
    protected String mateReferenceName = null;

    private LinkedHashMap<Short, Integer> tagNameToByteOffsets = null;
    private boolean needToRecomputeTagNameToByteOffsets = true;


    public FastBAMRecord(final SAMFileHeader header, final byte[] record) {
        this.mHeader = header;
        this.record = record;

        getFlags(); // decode flags so that flag bit getters don't have to check whether flags have been decoded
    }

    @Override
    public FastBAMRecord clone() throws CloneNotSupportedException {
        return new FastBAMRecord(mHeader, record.clone());
    }

    @Override
    public String toString() {
        return getSAMString();
    }



    //==================================================================
    //=========================== PRIVATE UTILS   ======================
    //==================================================================


    //offsets into record byte array
    private int readNameOffset() {
        return startOfVariableRegion;
    }

    private int cigarOffset() {
        return readNameOffset() + readNameSizeInBytes();
    }

    private int basesOffset() {
        return cigarOffset() + cigarSizeInBytes();
    }

    private int qualsOffset() {
        return basesOffset() + basesSizeInBytes();
    }

    private int tagsOffset() {
        return qualsOffset() + qualsSizeInBytes();
    }

    //size of variable-length record fields
    private int readNameSizeInBytes() {
        return getReadNameLength() + 1;
    }

    private int cigarSizeInBytes() {
        return getCigarLength() * 4;
    }

    private static int computeBaseSizeInBytes(final int readLength) {
        return (readLength + 1)/2;
    }

    private int basesSizeInBytes() {
        return computeBaseSizeInBytes(getReadLength());
    }

    private int qualsSizeInBytes() {
        return getReadLength();
    }

    private int tagsSizeInBytes() {
        return record.length - tagsOffset();
    }



    //==================================================================
    //================== PUBLIC GETTERS & SETTERS ======================
    //==================================================================

    public byte[] getRecord() {
        return record;
    }

    @Override
    public int getFlags() {
        if (mFlags == UNINITIALIZED)
            mFlags = decodeUnsignedShort(record, flagsOffset);

        return mFlags;
    }

    @Override
    public String getReadName() {
        if(mReadName == null)
            mReadName = new String(record, readNameOffset(), getReadNameLength());
        return mReadName;
    }

    @Override
    public void setReadName(final String value) {
        if(getReadNameLength() != value.length()) {
            resizeRecordField(readNameOffset(), readNameSizeInBytes(), value.length()+1);
            setReadNameLength((byte) value.length());
        }
        System.arraycopy(value.getBytes(), 0, record, readNameOffset(), value.length());
        record[readNameOffset() + value.length()] = 0; //null terminator

        mReadName = null;
    }

    @Override
    public int getReadNameLength() {
        final int length = record[readNameLengthOffset] & 0xFF;
        return length - 1; //subtract 1 to exclude null terminator
    }

    public void setReadNameLength(final int length) {
        if(length < 0 || length > 255) //max = 255 (unsigned byte)
            throw new IllegalArgumentException("invalid read name length: " + length);

        record[readNameLengthOffset] = (byte) (length + 1); //add 1 for null terminator
    }

    //number of cigar elements (number + operator) in the cigar string
    @Override
    public int getCigarLength() {
        if(n_cigar_ops == UNINITIALIZED) {
            n_cigar_ops = decodeUnsignedShort(record, cigarLenOffset);
        }
        return n_cigar_ops;
    }

    private void setCigarLength(final int n_cigar_ops) {
        this.n_cigar_ops = n_cigar_ops;
        encodeShort(record, cigarLenOffset, (short) n_cigar_ops);
    }

    @Override
    public int getMappingQuality() {
        return record[mappingQualityOffset] & 0xFF;
    }

    @Override
    public void setMappingQuality(final int value) {
        record[mappingQualityOffset] = (byte) value;
    }

    @Override
    public int getReadLength() {
        if(mReadLength == UNINITIALIZED)
            mReadLength = decodeInt(record, readLenOffset);

        return mReadLength;
    }

    private void setReadLength(final int readLength) {
        encodeInt(record, readLenOffset, readLength);
        this.mReadLength = readLength;
        mBaseQualities = null;
        mReadBases = null;
    }

    @Override
    public int getAlignmentStart() {
        if(mAlignmentStart == UNINITIALIZED)
            mAlignmentStart = decodeInt(record, coordinateOffset) + 1; //add 1 to convert from 0-based to 1-based coords

        return mAlignmentStart;
    }

    /**
     * @param value 1-based inclusive leftmost position of the clipped sequence, or 0 if there is no position.
     */
    @Override
    public void setAlignmentStart(final int value) {
        mAlignmentStart = value;
        encodeInt(record, coordinateOffset, value == -1? -1 : value - 1); //subtract 1 to convert from 0-based to 1-based coords
        mAlignmentEnd = UNINITIALIZED;
        // Change to cigar could change alignmentEnd, and thus indexing bin
        setIndexingBin(null);
    }

    @Override
    public int getAlignmentEnd() {
        if (getReadUnmappedFlag()) {
            return NO_ALIGNMENT_START;
        }
        else if (this.mAlignmentEnd == UNINITIALIZED) {
            this.mAlignmentEnd = mAlignmentStart + getCigarReferenceLength() - 1;
        }
        return this.mAlignmentEnd;
    }

    @Override
    public void setAlignmentEnd(final int value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public int getMateAlignmentStart() {
        return decodeInt(record, mateCoordinateOffset) + 1; //add 1 to convert from 0-based to 1-based coords
    }

    @Override
    public void setMateAlignmentStart(final int mateAlignmentStart) {
        encodeInt(record, mateCoordinateOffset, mateAlignmentStart - 1);
    }

    @Override
    public Integer getReferenceIndex() {
        if(referenceIndex == UNINITIALIZED)
            referenceIndex = decodeInt(record, referenceIdOffset);

        return referenceIndex;
    }

    @Override
    public Integer getMateReferenceIndex() {
        return decodeInt(record, mateReferenceIdOffset);
    }

    @Override
    public void setReferenceIndex(final int referenceIndex) {
        this.referenceName = null;
        encodeInt(record, referenceIdOffset, referenceIndex);
    }

    @Override
    public void setMateReferenceIndex(final int mateReferenceIndex) {
        this.mateReferenceName = null;
        encodeInt(record, mateReferenceIdOffset, mateReferenceIndex);
    }

    @Override
    public String getReferenceName() {
        if(referenceName == null) {
            final int refIndex = getReferenceIndex();
            if(refIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
                referenceName = NO_ALIGNMENT_REFERENCE_NAME;
            } else {
                try {
                    referenceName = mHeader.getSequence(refIndex).getSequenceName();
                } catch (final NullPointerException e) {
                    throw new IllegalArgumentException("Reference index " + getReferenceIndex() + " not found in sequence dictionary.", e);
                }
            }
        }

        return referenceName;
    }

    @Override
    public void setReferenceName(final String value) {
        //TODO make sure logic here is correct
        if (NO_ALIGNMENT_REFERENCE_NAME.equals(value)) {
            setReferenceIndex(NO_ALIGNMENT_REFERENCE_INDEX);
            referenceName = NO_ALIGNMENT_REFERENCE_NAME;
            return;
        } else if (mHeader != null) {
            final int referenceIndex = mHeader.getSequenceIndex(value);
            if (referenceIndex != NO_ALIGNMENT_REFERENCE_INDEX) {
                setReferenceIndex(referenceIndex);
                referenceName = value;
                return;
            }
        }
        // Drop through from above if nothing done.
        referenceName = value.intern();
        referenceIndex = UNINITIALIZED;
    }

    @Override
    public String getMateReferenceName() {
        if(mateReferenceName == null) {
            final int refIndex = getMateReferenceIndex();
            if(refIndex == NO_ALIGNMENT_REFERENCE_INDEX) {
                mateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
            } else {
                try {
                    mateReferenceName = mHeader.getSequence(refIndex).getSequenceName();
                } catch (final NullPointerException e) {
                    throw new IllegalArgumentException("Reference index " + getReferenceIndex() + " not found in sequence dictionary.", e);
                }
            }
        }

        return mateReferenceName;
    }

    @Override
    public void setMateReferenceName(final String value) {
        //TODO make sure logic here is correct
        if (NO_ALIGNMENT_REFERENCE_NAME.equals(value)) {
            setMateReferenceIndex(NO_ALIGNMENT_REFERENCE_INDEX);
            mateReferenceName = NO_ALIGNMENT_REFERENCE_NAME;
            return;
        } else if (mHeader != null) {
            final int index = mHeader.getSequenceIndex(value);
            if (index != NO_ALIGNMENT_REFERENCE_INDEX) {
                setMateReferenceIndex(index);
                mateReferenceName = value;
                return;
            }
        }
        // Drop through from above if nothing done.
        mateReferenceName = value.intern();
    }

    @Override
    public int getInferredInsertSize() {
        return decodeInt(record, cigarLenOffset);
    }

    @Override
    public void setInferredInsertSize(final int inferredInsertSize) {
        encodeInt(record, insertSizeOffset, inferredInsertSize);
    }

    @Override
    public void setReadString(final String value) {
        if (NULL_SEQUENCE_STRING.equals(value)) {
            setReadBases(NULL_SEQUENCE);
        } else {
            final byte[] bases = value.getBytes();
            SAMUtils.normalizeBases(bases);
            setReadBases(bases);
        }
    }

    @Override
    public void setBaseQualityString(final String value) {
        throw new RuntimeException("Not implemented");
    }

    private void resizeRecordField(final int fieldOffset, final int oldFieldSize, final int newFieldSize) {
        if(oldFieldSize == newFieldSize)
            throw new IllegalArgumentException("fieldSize == newFieldSize");

        final int oldRecordLength = record.length;
        final byte[] newRecord = new byte[oldRecordLength - oldFieldSize + newFieldSize];
        if(fieldOffset > 0)
            System.arraycopy(record, 0, newRecord, 0, fieldOffset);
        final int oldEndOfField = fieldOffset + oldFieldSize;
        if(oldEndOfField < oldRecordLength)
            System.arraycopy(record, oldEndOfField, newRecord, fieldOffset + newFieldSize, oldRecordLength - oldEndOfField);
        record = newRecord;
    }

    @Override
    public String getReadString() {
        if (getReadLength() == 0)
            return NULL_SEQUENCE_STRING;

        return StringUtil.bytesToString(getReadBases());
    }

    @Override
    public String getBaseQualityString() {
        if(!hasBaseQualities())
            return NULL_QUALS_STRING;


        return SAMUtils.phredToFastq(record, qualsOffset(), qualsSizeInBytes());
    }

    //==================================================================
    //=========================== CIGAR METHODS   ======================
    //==================================================================

    private int computeCigarOpOffset(final int i) {
        if(i < 0 || i >= getCigarLength())
            throw new NoSuchElementException(i + " is >= the number of operators in the cigar string: " + getCigarLength());

        return cigarOffset() + i*4;
    }

    public CigarOperator[] getCigarOps() {
        final CigarOperator[] cigarOps = new CigarOperator[getCigarLength()];
        for(int i = 0, cigarOffset = cigarOffset(); i < getCigarLength(); i++, cigarOffset += 4)
            cigarOps[i] = CigarOperator.binaryToEnum(record[cigarOffset] & 0xF);
        return cigarOps;
    }

    public CigarOperator getCigarOp(int i) {
        final int offset = computeCigarOpOffset(i);
        return CigarOperator.binaryToEnum(record[offset] & 0xF);
    }

    public int getCigarOpLength(final int i) {
        final int offset = computeCigarOpOffset(i);
        final int cigarette = decodeInt(record, offset);
        final int length = cigarette >>> 4; //TODO BinaryCigarCodec might have a bug here
        return length;
    }

    public void copyCigarFrom(final FastBAMRecord rec) {
        if (cigarSizeInBytes() != rec.cigarSizeInBytes()) {
            resizeRecordField(cigarOffset(), cigarSizeInBytes(), rec.cigarSizeInBytes());
            setCigarLength(rec.getCigarLength());
        }

        System.arraycopy(rec.record, rec.cigarOffset(), record, cigarOffset(), cigarSizeInBytes());
    }

    public void setCigarOp(final int i, final CigarOperator op) {
        final int offset = computeCigarOpOffset(i);
        record[offset] &= 0xF0; // clear the lowest 4 bits
        record[offset] |= CigarOperator.enumToBinary(op); //set the lowest 4 bits to op
        // Change to cigar could change alignmentEnd, and thus indexing bin
        setIndexingBin(null);
    }

    public void setCigarOpLength(final int i, int opLength) {
        final int offset = computeCigarOpOffset(i);
        opLength <<= 4;
        opLength |= record[offset] & 0xF;
        encodeInt(record, offset, opLength);
        // Change to cigar could change alignmentEnd, and thus indexing bin
        setIndexingBin(null);
    }

    /**
     * @return The number of reference bases that the read covers, excluding padding.
     */
    public int getCigarReferenceLength() {
        int length = 0;
        for(int i = 0; i < getCigarLength(); i++) {
            switch(getCigarOp(i))
            {
                case M:
                case D:
                case N:
                case EQ:
                case X:
                    length += getCigarOpLength(i);
            }
        }
        return length;
    }

    @Override
    public String getCigarString() {
        if(getCigarLength() == 0) {
            return ReadRecord.NO_ALIGNMENT_CIGAR;
        }
        StringBuffer sb = new StringBuffer(getCigarLength()*4);
        for(int i = 0; i < getCigarLength(); i++) {
            sb.append(getCigarOpLength(i));
            sb.append(getCigarOp(i));
        }
        return sb.toString();
    }



    //==================================================================
    //=========================== BASES, QUALS =========================
    //==================================================================

    @Override
    public byte[] getReadBases() {
        if (getReadLength() == 0)
            return NULL_SEQUENCE;

        if (mReadBases == null) {
            mReadBases = SAMUtils.compressedBasesToBytes(getReadLength(), record, basesOffset());
        }
        return mReadBases;
    }

    @Override
    public void setReadBases(final byte[] bases) {
        if(bases.length != getReadLength()) {
            resizeBasesAndQuals(bases.length);
        }

        mReadBases = bases;
        final byte[] compressedBytes = SAMUtils.bytesToCompressedBases(bases);
        System.arraycopy(compressedBytes, 0, record, basesOffset(), basesSizeInBytes());
    }

    private boolean hasBaseQualities() {
        // BAM files store missing qualities as an array of 0xFF bytes.
        // 0xFF is an illegal quality score value (it cannot be encoded in SAM)
        // and so the first byte is a suitable marker.
        // We hide this quirk of the BAM encoding so that the BAM interface looks the same as SAM.
        return getReadLength() > 0 && record[qualsOffset()] != (byte) 0xFF;
    }

    @Override
    public byte[] getBaseQualities() {
        if(mBaseQualities == null) {
            mBaseQualities = new byte[getReadLength()];
            System.arraycopy(record, qualsOffset(), mBaseQualities, 0, mBaseQualities.length);
        }
        return mBaseQualities;
    }

    @Override
    public void setBaseQualities(final byte[] quals) {
        if(quals.length != getReadLength())
            resizeBasesAndQuals(quals.length);

        System.arraycopy(quals, 0, record, basesOffset(), quals.length);
        mBaseQualities = quals;
    }

    // changes the width of the bases and quals fields simultaneously
    public void resizeBasesAndQuals(final int newReadLength) {
        resizeRecordField(basesOffset(), basesSizeInBytes()+qualsSizeInBytes(), computeBaseSizeInBytes(newReadLength)+newReadLength);
        setReadLength(newReadLength);
    }

    private static void reverseBytes(byte[] record, int firstByteOffset, int lastByteOffset) {
        while(firstByteOffset < lastByteOffset) {
            // exchange the first and last
            byte temp = record[firstByteOffset];
            record[firstByteOffset++] = record[lastByteOffset];
            record[lastByteOffset--] = temp;
        }
    }

    private static void reverseComplementBytes(byte[] record, int firstByteOffset, int lastByteOffset) {
        while(firstByteOffset < lastByteOffset) {
            // exchange and complement the first and last
            final byte temp = record[firstByteOffset];
            record[firstByteOffset++] = SequenceUtil.complement(record[lastByteOffset]);
            record[lastByteOffset--] = SequenceUtil.complement(temp);
        }
        if(firstByteOffset == lastByteOffset) {
            record[firstByteOffset] = SequenceUtil.complement(record[firstByteOffset]);
        }
    }

    public void reverseBaseQualities() {
        final int qualsOffset = qualsOffset();
        reverseBytes(record, qualsOffset, qualsOffset + qualsSizeInBytes() - 1);
        mBaseQualities = null;
    }

    private static final byte COMPRESSED_A_LOW = 1;
    private static final byte COMPRESSED_C_LOW = 2;
    private static final byte COMPRESSED_G_LOW = 4;
    private static final byte COMPRESSED_T_LOW = 8;
    private static final byte COMPRESSED_A_HIGH = COMPRESSED_A_LOW << 4;
    private static final byte COMPRESSED_C_HIGH = COMPRESSED_C_LOW << 4;
    private static final byte COMPRESSED_G_HIGH = COMPRESSED_G_LOW << 4;
    private static final byte COMPRESSED_T_HIGH = (byte)(COMPRESSED_T_LOW << 4);

    static byte complementLowNibble(final byte nibble) {
        switch(nibble) {
            case COMPRESSED_A_LOW: return COMPRESSED_T_LOW;
            case COMPRESSED_C_LOW: return COMPRESSED_G_LOW;
            case COMPRESSED_G_LOW: return COMPRESSED_C_LOW;
            case COMPRESSED_T_LOW: return COMPRESSED_A_LOW;
            default: return nibble;
        }
    }

    static byte complementHighNibble(final byte nibble) {
        switch(nibble) {
            case COMPRESSED_A_HIGH: return COMPRESSED_T_HIGH;
            case COMPRESSED_C_HIGH: return COMPRESSED_G_HIGH;
            case COMPRESSED_G_HIGH: return COMPRESSED_C_HIGH;
            case COMPRESSED_T_HIGH: return COMPRESSED_A_HIGH;
            default: return nibble;
        }
    }

    /* swap the high nibbles of the 2 bytes at the given offsets and (given that they represent ACGT bases) also complement them */
    static void swapAndComplementHighNibbles(final byte[] record, final int offset1, final int offset2) {
        final byte record_offset1 = record[offset1];
        final byte record_offset2 = record[offset2];
        record[offset2] = (byte) (complementHighNibble((byte) (record_offset1 & 0xF0)) | (record_offset2 & 0x0F));
        record[offset1] = (byte) (complementHighNibble((byte) (record_offset2 & 0xF0)) | (record_offset1 & 0x0F));
    }

    /* swap the low nibbles of the 2 bytes at the given offsets and (given that they represent ACGT bases) also complement them */
    static void swapAndComplementLowNibbles(final byte[] record, final int offset1, final int offset2) {
        final byte record_offset1 = record[offset1];
        final byte record_offset2 = record[offset2];
        record[offset2] = (byte) ((record_offset2 & 0xF0) | complementLowNibble((byte) (record_offset1 & 0x0F)));
        record[offset1] = (byte) ((record_offset1 & 0xF0) | complementLowNibble((byte) (record_offset2 & 0x0F)));
    }

    /* swap the nibbles in the given byte and (given that they represent ACGT bases) also complement them */
    static byte swapAndComplementNibbles(final byte b) {
        final byte baseToBeHighNibble = complementLowNibble((byte) (b & 0x0F));
        final byte baseToBeLowNibble = complementLowNibble((byte) ((b & 0xF0) >> 4));
        return (byte) ((baseToBeHighNibble << 4) | baseToBeLowNibble);
    }

    public void reverseComplementReadBases() {
        final int firstBaseOffset = basesOffset();
        final int lastBaseOffset = firstBaseOffset + basesSizeInBytes() - 1;
        final boolean lastBaseHasOnlyOneNibble = (getReadLength() % 2 != 0);
        if(lastBaseHasOnlyOneNibble) {
            int i, j;
            for(i = firstBaseOffset, j = lastBaseOffset; i < j; i++, j--)
                swapAndComplementHighNibbles(record, i, j);
            if(i == j) {
                final byte middle_base = record[i];
                record[i] = (byte) (complementHighNibble((byte) (middle_base & 0xF0)) | (middle_base & 0x0F)); //complement the middle nibble
            }

            for(i = firstBaseOffset, j = lastBaseOffset - 1; i < j; i++, j--)
                swapAndComplementLowNibbles(record, i, j);
            if(i == j) {
                final byte middle_base = record[i];
                record[i] = (byte) ((middle_base & 0xF0) | complementLowNibble((byte) (middle_base & 0x0F))); //complement the middle nibble
            }

         } else {
            reverseBytes(record, firstBaseOffset, lastBaseOffset - (lastBaseHasOnlyOneNibble ? 1 : 0));
            //because 2 bases are packed into the 2 nibbles of every byte, the above reverseByte(..) call produces a sequence where nibbles
            //in every byte are out of order. For example, for ACGTGC, the reverseByte(..) call produces GCGTAC instead of the need CGTGCA.
            for (int i = firstBaseOffset; i < firstBaseOffset + basesSizeInBytes(); i++) //this loop changes GCGTAC to CGTGCA
                record[i] = swapAndComplementNibbles(record[i]);
        }
        mReadBases = null;
    }



    //==================================================================
    //=========================== ATTRIBUTE TAGS      ==================
    //==================================================================

    @Override
    public void clearAttributes() {
        resizeRecordField(tagsOffset(), tagsSizeInBytes(), 0);
        needToRecomputeTagNameToByteOffsets = true;
    }

    private static int getFixedAttributeSizeInBytes(char tagType) {
        switch(tagType) {
            case 'A':
            case 'c':
            case 'C':
                return 1;
            case 's':
            case 'S':
                return 2;
            case 'f':
            case 'i':
            case 'I':
                return 4;
            default:
                throw new IllegalArgumentException("Unrecognized tag type: " + tagType);
        }
    }

    /**
     * Calculates the size of an attribute's value in bytes (needed for variable-length attributes, because their size isn't explicitly
     * stored in the record array.
     */
    private static int getAttributeValueSizeInBytes(final byte[] record, final char tagType, final int offsetToStartOfValue) {
        switch(tagType) {
            case 'Z':
                int size = 0;
                while (record[offsetToStartOfValue + size] != 0)
                    size++;
                return size + 1; //add 1 to include the null terminator in size
            case 'H':
                throw new IllegalArgumentException("Not implemented");
            case 'B':
                //include arrayType byte and length integer in returned size
                final char arrayType = (char) record[offsetToStartOfValue];
                final int length = decodeInt(record, offsetToStartOfValue + 1);
                return 1 + 4 + getFixedAttributeSizeInBytes(arrayType) * length; //add 1, 4 to include arrayType and length bytes in size
            default:
                return getFixedAttributeSizeInBytes(tagType);
        }
    }

    /**
     * Map each tag name to its byte offset in the tags field (eg. relative to tagsOffset()).
     * The byte offset points to the tag's type char (rather than the 2-char tag name which is already known).
     */
    protected void computeTagNameToByteOffsets() {
        if(tagNameToByteOffsets == null)
            tagNameToByteOffsets = new LinkedHashMap<Short, Integer>();
        else
            tagNameToByteOffsets.clear();

        final int tagsOffset = tagsOffset();
        int i = tagsOffset;
        while(i < record.length) {
            final short tag = decodeShort(record, i);
            i+= 2;
            tagNameToByteOffsets.put(tag, i - tagsOffset); //save this position (the position of the tagType byte) relative to tagsOffset
            char tagType = (char) record[i];
            i++;
            i += getAttributeValueSizeInBytes(record, tagType, i);
        }
        needToRecomputeTagNameToByteOffsets = false;
    }

    public boolean hasAttribute(final short tag) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();

        return tagNameToByteOffsets.containsKey(tag);
    }

    public Iterable<Short> getAttributesIterator() {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        return tagNameToByteOffsets.keySet();
    }

    public int getNumAttributes() {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        return tagNameToByteOffsets.size();
    }

    public Short[] getAttributeTags() {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        return tagNameToByteOffsets.keySet().toArray(new Short[getNumAttributes()]);
    }

    public char getAttributeType(short tag) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            throw new IllegalArgumentException("Tag not present:  " + SAMTagUtil.getSingleton().makeStringTag(tag));

        final int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        final char tagType = (char) record[i];
        return tagType;
    }

    public void deleteAttribute(final short tag) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            return;

        final int tagOffset = tagsOffset() + tagNameToByteOffsets.get(tag) - 2;
        final char tagType = (char) record[tagOffset + 2];
        final int tagTotalSize = 3 + getAttributeValueSizeInBytes(record, tagType, tagOffset + 3); //add 3 for tagType byte and byteOffset
        resizeRecordField(tagOffset, tagTotalSize, 0);
        needToRecomputeTagNameToByteOffsets = true;
    }

    private void addAttribute(final short tag, final char tagType, final Object value) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(tagNameToByteOffsets.containsKey(tag))
            throw new IllegalArgumentException("Tag already present:  " + SAMTagUtil.getSingleton().makeStringTag(tag));

        final int oldRecordLength = record.length;

        final int newByteCount = 2 + 1; //2 for the tag + 1 for tagType
        switch(tagType) {
            //case 'H': comments in BinaryTagCode say 'H' is not supported any more
            case 'Z':
                resizeRecordField(oldRecordLength, 0, newByteCount + ((String) value).length() + 1); //add 1 for the null terminator
                //fill string with non-zero values or else the zeros will be interpreted as null terminators
                Arrays.fill(record, oldRecordLength + newByteCount, oldRecordLength + newByteCount + ((String) value).length(), (byte) 0x7E);
                break;
            case 'H':
            case 'B':
                //include arrayType byte and length integer in returned size
                //final char arrayType = (char) record[offsetToStartOfValue];
                //final int length = decodeInt(record, offsetToStartOfValue + 1);
                //return 1 + 4 + getFixedAttributeSizeInBytes(arrayType) * length; //add 1, 4 to include arrayType and length bytes in size
                throw new RuntimeException("Not implemented");

            default:
                resizeRecordField(oldRecordLength, 0, newByteCount + getFixedAttributeSizeInBytes(tagType));
                break;
        }

        //set tag name and type bytes in the new record byte array
        encodeShort(record, oldRecordLength, tag);
        record[oldRecordLength + 2] = (byte) tagType;

        //update the tagNameToByteOffsets (since this tag was added to the end, other offsets don't need to be recalculated)
        tagNameToByteOffsets.put(tag, oldRecordLength + 2 - tagsOffset());

        //set the value into the byte array
        setExistingAttribute(tag, value);

    }

    private void encodeNullTerminatedString(final byte[] record, final int offset, final short tag, final String value) {
        final int existingValueLength = getAttributeValueSizeInBytes(record, 'Z', offset);
        if(value.length() + 1 == existingValueLength) {  //add 1 to include null terminator
            System.arraycopy(value.getBytes(), 0, record, offset, value.length());
        } else {
            deleteAttribute(tag);
            addAttribute(tag, 'Z', value);
        }
    }

    public void setAttribute(final short tag, final Object value) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();

        if(!tagNameToByteOffsets.containsKey(tag)) {
            if (value == null)
                return;
            final char tagType = BinaryTagCodec.getTagValueType(value);
            addAttribute(tag, tagType, value);
        } else if(value == null) {
            deleteAttribute(tag);
        } else {
            setExistingAttribute(tag, value);
        }
    }

    private static boolean integerTagType1FitsIntoType2(final char tagType1, final char tagType2) {
        switch(tagType1) {
            case 'c':
            case 'C':
                return tagType2 == 'c' || tagType2 == 'C';
            case 's':
            case 'S':
                return tagType2 == 'c' || tagType2 == 'C' || tagType2 == 's' || tagType2 == 'S';
            case 'i':
            case 'I':
                return tagType2 == 'c' || tagType2 == 'C' || tagType2 == 's' || tagType2 == 'S' || tagType2 == 'i' || tagType2 == 'I';
            default:
                return false;
        }
    }

    //better to use this method instead of setAttribute(final short tag, Object value)
    //so that, if the attribute needs to be added, type is already known and doens't have to be infered using BinaryTagCodec.getTagValueType
    public void setAttribute(final short tag, final char tagType, final Object value) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();

        if(!tagNameToByteOffsets.containsKey(tag)) {
            addAttribute(tag, tagType, value);
        } else if(value == null) {
            deleteAttribute(tag);
        } else {
            int i = tagsOffset() + tagNameToByteOffsets.get(tag);
            final char existingTagType = (char) record[i];
            if(existingTagType != tagType && !integerTagType1FitsIntoType2(tagType, existingTagType)) {
                throw new IllegalArgumentException("tagType argument ("+tagType+") for tag " + SAMTagUtil.getSingleton().makeStringTag(tag)+
                        " doesn't match this tag's existing tagType in the record (" + getAttributeType(tag) + ")");
            }
            setExistingAttribute(tag, value);
        }
    }

    public void copyAttributeFrom(final FastBAMRecord other, final short tag) {
        Object value = other.getAttribute(tag);
        char tagType = other.getAttributeType(tag);
        this.setAttribute(tag, tagType, value);
    }

    private void setExistingAttribute(final short tag, final Object value) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            throw new IllegalArgumentException("Tag not present:  " + SAMTagUtil.getSingleton().makeStringTag(tag));

        int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        final char tagType = (char) record[i++];
        switch(tagType) {
            case 'Z':
                encodeNullTerminatedString(record, i, tag, (String) value);
                break;
            case 'A':
            case 'c':
            case 'C':
                record[i] = ((Number) value).byteValue();
                break;
            case 's':
            case 'S':
                encodeShort(record, i, ((Number) value).shortValue());
                break;
            case 'i':
            case 'I':
                encodeInt(record, i, ((Number) value).intValue());
                break;
            case 'f':
                throw new RuntimeException("Not implemented");
            case 'H':
                throw new RuntimeException("Not implemented");
            case 'B':
                throw new RuntimeException("Not implemented");
            default:
                throw new IllegalArgumentException("Unrecognized tag type " + tagType);
        }
    }

    public Object getAttribute(final short tag) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            return null;

        int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        final char tagType = (char) record[i++];
        switch(tagType) {
            case 'Z':
                int Z_length = 0;
                while(record[i + Z_length] != 0)
                    Z_length++; //compute string length
                return new String(record, i, Z_length);
            case 'A':
                return (char) record[i];
            case 'c':
                return record[i];
            case 'C':
                return record[i] & 0xFF; //unsigned
            case 's':
                return decodeShort(record, i);
            case 'S':
                return decodeUnsignedShort(record, i);
            case 'i':
                return decodeInt(record, i);
            case 'I':
                return decodeUnsignedInt(record, i);
            case 'f':
                return decodeFloat(record, i);
            case 'H':
                throw new RuntimeException("Not implemented");
            case 'B':
                throw new RuntimeException("Not implemented");
            default:
                throw new IllegalArgumentException("Unrecognized tag type " + tagType);
        }
    }

    public Integer getIntegerAttribute(final short tag) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            return null;
        final int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        return decodeInt(record, i + 1);
    }

    public void setIntegerAttribute(final short tag, final int value) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            throw new IllegalArgumentException("Tag not present:  " + SAMTagUtil.getSingleton().makeStringTag(tag));
        int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        encodeInt(record, i + 1, value);
    }

    public void reverseOrReverseComplementStringAttribute(final short tag, final boolean reverseComplement) {
        if(needToRecomputeTagNameToByteOffsets)
            computeTagNameToByteOffsets();
        if(!tagNameToByteOffsets.containsKey(tag))
            throw new IllegalArgumentException("Tag not present:  " + SAMTagUtil.getSingleton().makeStringTag(tag));
        int i = tagsOffset() + tagNameToByteOffsets.get(tag);
        if((char) record[i] != 'Z')
            throw new IllegalArgumentException("Tag type is: " + (char) record[i] + " instead of 'Z'");
        i++;
        final int size = getAttributeValueSizeInBytes(record, 'Z', i);
        final int lastByteOffset = i + size - 2; //subtract 2 to make this offset point to the byte before the null terminator
        if(reverseComplement)
            reverseComplementBytes(record, i, lastByteOffset);
        else
            reverseBytes(record, i, lastByteOffset);
    }



    //========================================================
    //======================== MISC METHODS ==================
    //========================================================

    @Override
    public List<AlignmentBlock> getAlignmentBlocks() {
        final int cigarLength = getCigarLength();
        if(cigarLength == 0)
            return Collections.emptyList();

        final List<AlignmentBlock> alignmentBlocks = new ArrayList<AlignmentBlock>();
        int readBase = 1;
        int refBase  = getAlignmentStart();

        for(int i = 0; i < cigarLength; i++) {
            final CigarOperator op = getCigarOp(i);
            switch(op) {
                case H : break; // ignore hard clips
                case P : break; // ignore pads
                case S : readBase += getCigarOpLength(i); break; // soft clip read bases
                case N : refBase += getCigarOpLength(i); break;  // reference skip
                case D : refBase += getCigarOpLength(i); break;
                case I : readBase += getCigarOpLength(i); break;
                case M :
                case EQ :
                case X :
                    final int length = getCigarOpLength(i);
                    alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length));
                    readBase += length;
                    refBase  += length;
                    break;
                default :
                    throw new IllegalStateException("Case statement didn't deal with op: " + op);
            }
        }
        return Collections.unmodifiableList(alignmentBlocks); //TODO cache this?
    }


    @Override
    public SAMReadGroupRecord getReadGroup() {
        if(!hasAttribute(SAMTagUtil.RG))
            return null;
        final String rgId = (String)getAttribute(SAMTagUtil.RG);
        return getHeader().getReadGroup(rgId);
    }

    //========================================================
    //=============== UNIMPLEMENTED METHODS ==================
    //========================================================

    @Override
    public void setValidationStringency(final ValidationStringency validationStringency) {
        //TODO validation?
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
    public Cigar getCigar() {
        throw new IllegalStateException("Not Implemented");
    }

    @Override
    public ValidationStringency getValidationStringency() {
        throw new RuntimeException("Not implemented");
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
    public void setAttribute(final String tag, final Object value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setUnsignedArrayAttribute(final String tag, final Object value) {
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
        throw new RuntimeException("Will not implemented");
    }

    @Override
    public List<SAMValidationError> isValid() {
        throw new RuntimeException("Will not be implemented");
    }

    @Override
    public byte[] getOriginalBaseQualities() {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setOriginalBaseQualities(final byte[] oq) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setFileSource(final SAMFileSource fileSource) {
        throw new RuntimeException("Will not be implemented");
    }

    @Override
    public void eagerDecode() {
        throw new RuntimeException("Will not be implemented");
    }

    @Override
    public void setCigarString(final String value) {
        throw new RuntimeException("Will not be implemented");
    }

    @Override
    public void setCigar(final Cigar cigar) {
        throw new RuntimeException("Will not be implemented");
    }
}

