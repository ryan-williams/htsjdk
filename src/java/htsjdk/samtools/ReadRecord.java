package htsjdk.samtools;

import java.util.List;

/**
 * Created by weisburd on 8/18/14.
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

    /**
     * This method is preferred over getReadName().length(), because for BAMRecord
     * it may be faster.
     * @return length not including a null terminator.
     */
    int getReadNameLength();

    void setReadName(String value);

    /**
     * @return read sequence as a string of ACGTN=.
     */
    String getReadString();

    void setReadString(String value);

    /**
     * Do not modify the value returned by this method.  If you want to change the bases, create a new
     * byte[] and call setReadBases() or call setReadString().
     * @return read sequence as ASCII bytes ACGTN=.
     */
    byte[] getReadBases();

    void setReadBases(byte[] value);

    /**
     * This method is preferred over getReadBases().length, because for BAMRecord it may be faster.
     * @return number of bases in the read.
     */
    int getReadLength();

    /**
     * @return Base qualities, encoded as a FASTQ string.
     */
    String getBaseQualityString();

    void setBaseQualityString(String value);

    /**
     * Do not modify the value returned by this method.  If you want to change the qualities, create a new
     * byte[] and call setBaseQualities() or call setBaseQualityString().
     * @return Base qualities, as binary phred scores (not ASCII).
     */
    byte[] getBaseQualities();

    void setBaseQualities(byte[] value);

    /**
     * If the original base quality scores have been store in the "OQ" tag will return the numeric
     * score as a byte[]
     */
    byte[] getOriginalBaseQualities();

    /**
     * Sets the original base quality scores into the "OQ" tag as a String.  Supplied value should be
     * as phred-scaled numeric qualities.
     */
    void setOriginalBaseQualities(byte[] oq);

    /**
     * @return Reference name, or null if record has no reference.
     */
    String getReferenceName();

    void setReferenceName(String value);

    /**
     * @return index of the reference sequence for this read in the sequence dictionary, or -1
     * if read has no reference sequence set, or if a String reference name is not found in the sequence index..
     */
    Integer getReferenceIndex();

    /**
     * @param referenceIndex Must either equal -1 (indicating no reference), or exist in the sequence dictionary
     * in the header associated with this record.
     */
    void setReferenceIndex(int referenceIndex);

    /**
     * @return index of the reference sequence for this read's mate in the sequence dictionary, or -1
     * if mate has no reference sequence set.
     */
    Integer getMateReferenceIndex();

    /**
     * @param referenceIndex Must either equal -1 (indicating no reference), or exist in the sequence dictionary
     * in the header associated with this record.
     */
    void setMateReferenceIndex(int referenceIndex);

    /**
     * @return 1-based inclusive leftmost position of the clipped sequence, or 0 if there is no position.
     */
    int getAlignmentStart();

    /**
     * @param value 1-based inclusive leftmost position of the clipped sequence, or 0 if there is no position.
     */
    void setAlignmentStart(int value);

    /**
     * @return 1-based inclusive rightmost position of the clipped sequence, or 0 read if unmapped.
     */
    int getAlignmentEnd();

    /**
     * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     *
     * Invalid to call on an unmapped read.
     */
    int getUnclippedStart();

    /**
     * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     *
     * Invalid to call on an unmapped read.
     */
    int getUnclippedEnd();

    /**
     * @return 1-based inclusive reference position of the unclipped sequence at a given offset,
     *         or 0 if there is no position.
     *         For example, given the sequence NNNAAACCCGGG, cigar 3S9M, and an alignment start of 1,
     *         and a (1-based)offset 10 (start of GGG) it returns 7 (1-based offset starting after the soft clip.
     *         For example: given the sequence AAACCCGGGTTT, cigar 4M1D6M, an alignment start of 1,
     *         an offset of 4 returns reference position 4, an offset of 5 returns reference position 6.
     *         Another example: given the sequence AAACCCGGGTTT, cigar 4M1I6M, an alignment start of 1,
     *         an offset of 4 returns reference position 4, an offset of 5 returns 0.
     * @offset 1-based location within the unclipped sequence
     */
    int getReferencePositionAtReadPosition(int offset);

    /**
     * Unsupported.  This property is derived from alignment start and CIGAR.
     */
    void setAlignmentEnd(int value);

    /**
     * @return 1-based inclusive leftmost position of the clipped mate sequence, or 0 if there is no position.
     */
    int getMateAlignmentStart();

    void setMateAlignmentStart(int mateAlignmentStart);

    /**
     * @return insert size (difference btw 5' end of read & 5' end of mate), if possible, else 0.
     * Negative if mate maps to lower position than read.
     */
    int getInferredInsertSize();

    void setInferredInsertSize(int inferredInsertSize);

    /**
     * @return phred scaled mapping quality.  255 implies valid mapping but quality is hard to compute.
     */
    int getMappingQuality();

    void setMappingQuality(int value);

    String getCigarString();

    void setCigarString(String value);

    /**
     * Do not modify the value returned by this method.  If you want to change the Cigar, create a new
     * Cigar and call setCigar() or call setCigarString()
     * @return Cigar object for the read, or null if there is none.
     */
    Cigar getCigar();

    /**
     * This method is preferred over getCigar().getNumElements(), because for BAMRecord it may be faster.
     * @return number of cigar elements (number + operator) in the cigar string.
     */
    int getCigarLength();

    void setCigar(Cigar cigar);

    /**
     * Get the SAMReadGroupRecord for this SAMRecord.
     * @return The SAMReadGroupRecord from the SAMFileHeader for this SAMRecord, or null if
     * 1) this record has no RG tag, or 2) the header doesn't contain the read group with
     * the given ID.
     * @throws NullPointerException if this.getHeader() returns null.
     * @throws ClassCastException if RG tag does not have a String value.
     */
    SAMReadGroupRecord getReadGroup();

    /**
     * It is preferable to use the get*Flag() methods that handle the flag word symbolically.
     */
    int getFlags();

    void setFlags(int value);

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    boolean getReadPairedFlag();

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    boolean getProperPairFlag();

    /**
     * the query sequence itself is unmapped.
     */
    boolean getReadUnmappedFlag();

    /**
     * the mate is unmapped.
     */
    boolean getMateUnmappedFlag();

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    boolean getReadNegativeStrandFlag();

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    boolean getMateNegativeStrandFlag();

    /**
     * the read is the first read in a pair.
     */
    boolean getFirstOfPairFlag();

    /**
     * the read is the second read in a pair.
     */
    boolean getSecondOfPairFlag();

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    boolean getNotPrimaryAlignmentFlag();

    /**
     * the alignment is supplementary (TODO: further explanation?).
     */
    boolean getSupplementaryAlignmentFlag();

    /**
     * the read fails platform/vendor quality checks.
     */
    boolean getReadFailsVendorQualityCheckFlag();

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    boolean getDuplicateReadFlag();

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    void setReadPairedFlag(boolean flag);

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    void setProperPairFlag(boolean flag);

    /**
     * the query sequence itself is unmapped.  This method name is misspelled.
     * Use setReadUnmappedFlag instead.
     * @deprecated
     */
    void setReadUmappedFlag(boolean flag);

    /**
     * the query sequence itself is unmapped.
     */
    void setReadUnmappedFlag(boolean flag);

    /**
     * the mate is unmapped.
     */
    void setMateUnmappedFlag(boolean flag);

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    void setReadNegativeStrandFlag(boolean flag);

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    void setMateNegativeStrandFlag(boolean flag);

    /**
     * the read is the first read in a pair.
     */
    void setFirstOfPairFlag(boolean flag);

    /**
     * the read is the second read in a pair.
     */
    void setSecondOfPairFlag(boolean flag);

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    void setNotPrimaryAlignmentFlag(boolean flag);

    /**
     * the alignment is supplementary (TODO: further explanation?).
     */
    void setSupplementaryAlignmentFlag(boolean flag);

    /**
     * the read fails platform/vendor quality checks.
     */
    void setReadFailsVendorQualityCheckFlag(boolean flag);

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    void setDuplicateReadFlag(boolean flag);

    /**
     * Tests if this record is either a secondary and/or supplementary alignment;
     * equivalent to {@code (getNotPrimaryAlignmentFlag() || getSupplementaryAlignmentFlag())}.
     */
    boolean isSecondaryOrSupplementary();

    ValidationStringency getValidationStringency();

    /**
     * Control validation of lazily-decoded elements.
     */
    void setValidationStringency(ValidationStringency validationStringency);

    /**
     * Get the value for a SAM tag.
     * WARNING: Some value types (e.g. byte[]) are mutable.  It is dangerous to change one of these values in
     * place, because some SAMRecord implementations keep track of when attributes have been changed.  If you
     * want to change an attribute value, call setAttribute() to replace the value.
     *
     * @param tag Two-character tag name.
     * @return Appropriately typed tag value, or null if the requested tag is not present.
     */
    Object getAttribute(String tag);

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into an Integer if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in an Integer.
     */
    Integer getIntegerAttribute(String tag);

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into a Short if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in a Short.
     */
    Short getShortAttribute(String tag);

    /**
     * Get the tag value and attempt to coerce it into the requested type.
     * @param tag The requested tag.
     * @return The value of a tag, converted into a Byte if possible.
     * @throws RuntimeException If the value is not an integer type, or will not fit in a Byte.
     */
    Byte getByteAttribute(String tag);

    String getStringAttribute(String tag);

    Character getCharacterAttribute(String tag);

    Float getFloatAttribute(String tag);

    /** Will work for signed byte array, unsigned byte array, or old-style hex array */
    byte[] getByteArrayAttribute(String tag);

    byte[] getUnsignedByteArrayAttribute(String tag);

    /** Will work for signed byte array or old-style hex array */
    byte[] getSignedByteArrayAttribute(String tag);

    short[] getUnsignedShortArrayAttribute(String tag);

    short[] getSignedShortArrayAttribute(String tag);

    int[] getUnsignedIntArrayAttribute(String tag);

    int[] getSignedIntArrayAttribute(String tag);

    float[] getFloatArrayAttribute(String tag);

    /**
     * @return True if this tag is an unsigned array, else false.
     * @throws htsjdk.samtools.SAMException if the tag is not present.
     */
    boolean isUnsignedArrayAttribute(String tag);

    /**
     * @see htsjdk.samtools.SAMRecord#getAttribute(String)
     * @param tag Binary representation of a 2-char String tag as created by SAMTagUtil.
     */
    Object getAttribute(short tag);

    /**
     * Set a named attribute onto the SAMRecord.  Passing a null value causes the attribute to be cleared.
     * @param tag two-character tag name.  See http://samtools.sourceforge.net/SAM1.pdf for standard and user-defined tags.
     * @param value Supported types are String, Char, Integer, Float, byte[], short[]. int[], float[].
     * If value == null, tag is cleared.
     *
     * Byte and Short are allowed but discouraged.  If written to a SAM file, these will be converted to Integer,
     * whereas if written to BAM, getAttribute() will return as Byte or Short, respectively.
     *
     * Long with value between 0 and MAX_UINT is allowed for BAM but discouraged.  Attempting to write such a value
     * to SAM will cause an exception to be thrown.
     *
     * To set unsigned byte[], unsigned short[] or unsigned int[] (which is discouraged because of poor Java language
     * support), setUnsignedArrayAttribute() must be used instead of this method.
     *
     * String values are not validated to ensure that they conform to SAM spec.
     */
    void setAttribute(String tag, Object value);

    /**
     * Because Java does not support unsigned integer types, we think it is a bad idea to encode them in SAM
     * files.  If you must do so, however, you must call this method rather than setAttribute, because calling
     * this method is the way to indicate that, e.g. a short array should be interpreted as unsigned shorts.
     * @param value must be one of byte[], short[], int[]
     */
    void setUnsignedArrayAttribute(String tag, Object value);

    /**
     * Removes all attributes.
     */
    void clearAttributes();

    /**
     * @return list of {tag, value} tuples
     */
    List<SAMTagAndValue> getAttributes();

    SAMFileHeader getHeader();

    /**
     * Setting header into SAMRecord facilitates conversion btw reference sequence names and indices
     * @param header contains sequence dictionary for this SAMRecord
     */
    void setHeader(SAMFileHeader header);

    /**
     * If this record has a valid binary representation of the variable-length portion of a binary record stored,
     * return that byte array, otherwise return null.  This will never be true for SAMRecords.  It will be true
     * for BAMRecords that have not been eagerDecoded(), and for which none of the data in the variable-length
     * portion has been changed.
     */
    byte[] getVariableBinaryRepresentation();

    /**
     * Depending on the concrete implementation, the binary file size of attributes may be known without
     * computing them all.
     * @return binary file size of attribute, if known, else -1
     */
    int getAttributesBinarySize();

    /**
     *
     * @return String representation of this.
     * @deprecated This method is not guaranteed to return a valid SAM text representation of the SAMRecord.
     * To get standard SAM text representation, use htsjdk.samtools.SAMRecord#getSAMString().
     */
    String format();

    /**
     * Returns blocks of the read sequence that have been aligned directly to the
     * reference sequence. Note that clipped portions of the read and inserted and
     * deleted bases (vs. the reference) are not represented in the alignment blocks.
     */
    List<AlignmentBlock> getAlignmentBlocks();

    /**
     * Run all validations of CIGAR.  These include validation that the CIGAR makes sense independent of
     * placement, plus validation that CIGAR + placement yields all bases with M operator within the range of the reference.
     * @param recordNumber For error reporting.  -1 if not known.
     * @return List of errors, or null if no errors.
     */
    List<SAMValidationError> validateCigar(long recordNumber);

    @Override
    boolean equals(Object o);

    @Override
    int hashCode();

    /**
     * Perform various validations of SAMRecord.
     * Note that this method deliberately returns null rather than Collections.emptyList() if there
     * are no validation errors, because callers tend to assume that if a non-null list is returned, it is modifiable.
     * @return null if valid.  If invalid, returns a list of error messages.
     */
    List<SAMValidationError> isValid();

    /**
     * Gets the source of this SAM record -- both the reader that retrieved the record and the position on disk from
     * whence it came.
     * @return The file source.  Note that the reader will be null if not activated using SAMFileReader.enableFileSource().
     */
    SAMFileSource getFileSource();

    /**
     * Note that this does a shallow copy of everything, except for the attribute list, for which a copy of the list
     * is made, but the attributes themselves are copied by reference.  This should be safe because callers should
     * never modify a mutable value returned by any of the get() methods anyway.
     */
    @Override
    Object clone() throws CloneNotSupportedException;

    /** Simple toString() that gives a little bit of useful info about the read. */
    @Override
    String toString();

    /**
     Returns the record in the SAM line-based text format.  Fields are
     separated by '\t' characters, and the String is terminated by '\n'.
     */
    String getSAMString();

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
