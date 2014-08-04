package htsjdk.samtools;
 
import java.util.List;
 
/**
 * Other classes to be turned into Interfaces:
 *    1) Cigar (and friends)
 *    2) SAMFileHeader -> ReadSetHeader(?) and related objects (e.g. SAMReadGroupRecord)
 *
 * @author nhomer
 */
public interface ReadRecord {
    /**
     * Alignment score for a good alignment, but where computing a Phred-score is not feasible.
     */
    int MAPQ_UNKNOWN = 255;
    /**
     * If a read has this reference name, it is unaligned, but not all unaligned reads have
     * this reference name (see above).
     */
    String REFERENCE_NAME_NO_ALIGNMENT = "*";
    /**
     * If a read has this reference index, it is unaligned, but not all unaligned reads have
     * this reference index (see above).
     */
    int REFERENCE_INDEX_NO_ALIGNMENT = -1;
    /**
     * Cigar string for an unaligned read.
     */
    String CIGAR_NO_ALIGNMENT = "*";
    /**
     * If a read has reference name "*", it will have this value for position.
     */
    int NO_ALIGNMENT_START = GenomicIndexUtil.UNSET_GENOMIC_LOCATION;
 
    ///////////////////////////////////////////////////////////////////////////
    // Core read fields
    ///////////////////////////////////////////////////////////////////////////
    String getReadName();
    void setReadName(String value); // NH: OK
 
    byte[] getReadBases(); // NH: OK
    void setReadBases(byte[] value); // NH: OK
    int getReadLength(); // NH: OK
 
    byte[] getBaseQualities(); // NH: rename getBaseQualities
    void setBaseQualities(byte[] value); // NH: rename setBaseQualities
 
    String getReferenceName(); // NH: OK
    void setReferenceName(String value); // NH: OK
 
    int getReferenceIndex(); // NH: OK, but relies on tying this to a header
    void setReferenceIndex(int referenceIndex); // NH: OK
 
    int getAlignmentStart(); //  NH: OK
    void setAlignmentStart(int value);// NH: OK
 
    int getAlignmentEnd(); //  NH: OK
    void setAlignmentEnd(int value); // NH: OK
 
    int getInsertSize(); // NH: OK
    void setInsertSize(int insertSize);// NH: OK
 
    int getMappingQuality(); // NH: OK
    void setMappingQuality(int value); // NH: OK
    
    /** Returns true if the record has a valid reference sequence and alignment start position; false otherwise. */
    boolean hasPositionInformation();
 
    /** Returns true if the read is aligned and has a mapping quality value other than MAPPING_QUALITY_UNKNOWN. */
    boolean hasMappingQuality();
 
    Cigar getCigar();
    void setCigar(Cigar cigar);
    int getCigarLength();
    CigarElement getCigarOperator(int i);
    int getCigarOperatorLength(int i);
 
    ///////////////////////////////////////////////////////////////////////////
    // Mate fields
    ///////////////////////////////////////////////////////////////////////////
    String getMateReferenceName(); // NH: OK, but do we want to change the mate to something else?
    void setMateReferenceName(String mateReferenceName); // NH: OK
 
    int getMateReferenceIndex(); // NH: OK
    void setMateReferenceIndex(int referenceIndex); // NH: OK
 
    int getMateAlignmentStart(); // NH: OK
    void setMateAlignmentStart(int mateAlignmentStart); // NH: OK
 
    /** Returns true if the record has a valid mate reference sequence and mate alignment start position; false otherwise. */
    boolean hasMatePositionInformation();
 
 
 
    SAMReadGroupRecord getReadGroup(); // NH: too specific, remove
 
 
 
    boolean isPaired();
    void   setPaired(boolean paired);
 
    boolean isProperlyPaired();
    void   setProperlyPaired(boolean properlyPaired);
 
    boolean isMapped();
    void   setMapped(boolean mapped);
 
    boolean isMateMapped();
    void   setMateMapped(boolean mateMapped);
 
    boolean isPositiveStrand();
    void   setPositiveStrand(boolean positiveStrand);
 
    boolean isMatePositiveStrand();
    void   setMatePositiveStrand(boolean matePositiveStrand);
 
    boolean isFirstOfPair();
    void   setFirstOfPair(boolean firstOfPair);
 
    boolean isSecondOfPair();
    void   setSecondOfPair(boolean secondOfPair);
    
    /** Returns true if the read is either the sole read from a template OR is the first read from a template with multiple reads. */
    boolean isFirstOrFrag(); // Decide if we need a better name
 
    boolean isPrimaryAlignment();
    void   setPrimaryAlignment(boolean primaryAlignment);
 
    boolean isSupplementaryRecord();
    void   setSupplementaryRecord(boolean supplementary);
 
    boolean passesFilter();
    void setPassesFilter(boolean passesFilter);
 
    boolean isDuplicate();
    void setDuplicate(boolean duplicate);
 
    /** Returns true if the record represents either a non-primary alignment, or a supplementary record for the primary alignment. */
    boolean isSecondaryOrSupplementary();
 
 
    boolean hasAttribute(String tag);
 
    Integer getIntegerAttribute(String tag);
    int getIntegerAttribute(String tag, int defaultValue);
 
    Short getShortAttribute(String tag);
    short getShortAttribute(String tag, short defaultValue);
 
    Byte getByteAttribute(String tag);
    byte getByteAttribute(String tag, byte defaultValue);
 
    String getStringAttribute(String tag);
    String getStringAttribute(String tag, String defaultValue);
 
    Character getCharacterAttribute(String tag);
    char getCharacterAttribute(String tag, char defaultValue);
 
    Float getFloatAttribute(String tag);
    float getFloatAttribute(String tag, float defaultValue);
 
    // TF: do we really want to support unsigned array attributes?  Ugh.
    byte[] getByteArrayAttribute(String tag);
 
    byte[] getUnsignedByteArrayAttribute(String tag);
 
    byte[] getSignedByteArrayAttribute(String tag);
 
    short[] getUnsignedShortArrayAttribute(String tag);
 
    short[] getSignedShortArrayAttribute(String tag);
 
    int[] getUnsignedIntArrayAttribute(String tag);
 
    int[] getSignedIntArrayAttribute(String tag);
 
    float[] getFloatArrayAttribute(String tag);
 
    // TF: do we want to support the addition of arbitrary untyped attributes?  What's the contract of this method
    void setAttribute(String tag, Object value);
    Object getAttribute(short tag); // NH: remove, too specific
 
    void setUnsignedArrayAttribute(String tag, Object value);
 
    void clearAttributes();
 
    List<SAMTagAndValue> getAttributes(); // NH: make this generic?
 
    SAMFileHeader getHeader(); // NH: can we get away without a header?
 
    void setHeader(SAMFileHeader header);  // NH: OK
 
    byte[] getVariableBinaryRepresentation(); // NH: too specific, remove
 
    List<AlignmentBlock> getAlignmentBlocks(); // NH: move to utility classes
 
    List<SAMValidationError> validateCigar(long recordNumber); // NH: too specific, remove
 
    @Override
    boolean equals(Object o); // NH: OK
 
    @Override
    int hashCode(); // NH: too specific, remove
 
    List<SAMValidationError> isValid(); // NH: too specific, remove
 
    SAMFileSource getFileSource(); // NH: too specific, remove
 
    @Override
    String toString(); // NH: OK
 
    String getSAMString(); // NH: too specific, remove
 
    // NH: too specific, remove
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
