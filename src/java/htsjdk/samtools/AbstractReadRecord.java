package htsjdk.samtools;

/**
 * Abstract ReadRecord class that contains implementations of methods common to all ReadRecord subtyptes.
 */
public abstract class AbstractReadRecord implements ReadRecord {

    protected SAMFileHeader mHeader = null;
    protected int mFlags = -1;
    protected Integer mIndexingBin = null;

    @Override
    public SAMFileHeader getHeader() {
        return mHeader;
    }

    /**
     * Setting header into SAMRecord facilitates conversion btw reference sequence names and indices
     * @param header contains sequence dictionary for this SAMRecord
     */
    @Override
    public void setHeader(final SAMFileHeader header) {
        this.mHeader = header;
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        return (AbstractReadRecord) super.clone();
    }

    /**
     Returns the record in the SAM line-based text format.  Fields are
     separated by '\t' characters, and the String is terminated by '\n'.
     */
    @Override
    public String getSAMString() {
        return SAMTextWriter.getSAMString(this);
    }


    @Override
    public Integer getIndexingBin() {
        return mIndexingBin;
    }

    /**
     * Used internally when writing BAMRecords.
     * @param mIndexingBin c.f. http://samtools.sourceforge.net/SAM1.pdf
     */
    @Override
    public void setIndexingBin(final Integer mIndexingBin) {
        this.mIndexingBin = mIndexingBin;
    }

    /**
     * Does not change state of this.
     * @return indexing bin based on alignment start & end.
     */
    @Override
    public int computeIndexingBin() {
        // reg2bin has zero-based, half-open API
        final int alignmentStart = getAlignmentStart()-1;
        int alignmentEnd = getAlignmentEnd();
        if (alignmentEnd <= 0) {
            // If alignment end cannot be determined (e.g. because this read is not really aligned),
            // then treat this as a one base alignment for indexing purposes.
            alignmentEnd = alignmentStart + 1;
        }
        return GenomicIndexUtil.reg2bin(alignmentStart, alignmentEnd);
    }

    @Override
    public void setFlags(final int value) {
        mFlags = value;
        // Could imply change to readUnmapped flag, which could change indexing bin
        setIndexingBin(null);
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    @Override
    public boolean getReadPairedFlag() {
        return (mFlags & READ_PAIRED_FLAG) != 0;
    }

    private void requireReadPaired() {
        if (!getReadPairedFlag()) {
            throw new IllegalStateException("Inappropriate call if not paired read");
        }
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    @Override
    public boolean getProperPairFlag() {
        requireReadPaired();
        return getProperPairFlagUnchecked();
    }

    protected boolean getProperPairFlagUnchecked() {
        return (mFlags & PROPER_PAIR_FLAG) != 0;
    }

    /**
     * the query sequence itself is unmapped.
     */
    @Override
    public boolean getReadUnmappedFlag() {
        return (mFlags & READ_UNMAPPED_FLAG) != 0;
    }

    /**
     * the mate is unmapped.
     */
    @Override
    public boolean getMateUnmappedFlag() {
        requireReadPaired();
        return getMateUnmappedFlagUnchecked();
    }

    protected boolean getMateUnmappedFlagUnchecked() {
        return (mFlags & MATE_UNMAPPED_FLAG) != 0;
    }

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    @Override
    public boolean getReadNegativeStrandFlag() {
        return (mFlags & READ_STRAND_FLAG) != 0;
    }

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    @Override
    public boolean getMateNegativeStrandFlag() {
        requireReadPaired();
        return getMateNegativeStrandFlagUnchecked();
    }

    protected boolean getMateNegativeStrandFlagUnchecked() {
        return (mFlags & MATE_STRAND_FLAG) != 0;
    }

    /**
     * the read is the first read in a pair.
     */
    @Override
    public boolean getFirstOfPairFlag() {
        requireReadPaired();
        return getFirstOfPairFlagUnchecked();
    }

    protected boolean getFirstOfPairFlagUnchecked() {
        return (mFlags & FIRST_OF_PAIR_FLAG) != 0;
    }

    /**
     * the read is the second read in a pair.
     */
    @Override
    public boolean getSecondOfPairFlag() {
        requireReadPaired();
        return getSecondOfPairFlagUnchecked();
    }

    protected boolean getSecondOfPairFlagUnchecked() {
        return (mFlags & SECOND_OF_PAIR_FLAG) != 0;
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    @Override
    public boolean getNotPrimaryAlignmentFlag() {
        return (mFlags & NOT_PRIMARY_ALIGNMENT_FLAG) != 0;
    }

    /**
     * the alignment is supplementary (TODO: further explanation?).
     */
    @Override
    public boolean getSupplementaryAlignmentFlag() {
        return (mFlags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    /**
     * the read fails platform/vendor quality checks.
     */
    @Override
    public boolean getReadFailsVendorQualityCheckFlag() {
        return (mFlags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    @Override
    public boolean getDuplicateReadFlag() {
        return (mFlags & DUPLICATE_READ_FLAG) != 0;
    }

    /**
     * the read is paired in sequencing, no matter whether it is mapped in a pair.
     */
    @Override
    public void setReadPairedFlag(final boolean flag) {
        setFlag(flag, READ_PAIRED_FLAG);
    }

    /**
     * the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment).
     */
    @Override
    public void setProperPairFlag(final boolean flag) {
        setFlag(flag, PROPER_PAIR_FLAG);
    }

    /**
     * the query sequence itself is unmapped.  This method name is misspelled.
     * Use setReadUnmappedFlag instead.
     * @deprecated
     */
    @Override
    public void setReadUmappedFlag(final boolean flag) {
        setReadUnmappedFlag(flag);
    }

    /**
     * the query sequence itself is unmapped.
     */
    @Override
    public void setReadUnmappedFlag(final boolean flag) {
        setFlag(flag, READ_UNMAPPED_FLAG);
        // Change to readUnmapped could change indexing bin
        setIndexingBin(null);
    }

    /**
     * the mate is unmapped.
     */
    @Override
    public void setMateUnmappedFlag(final boolean flag) {
        setFlag(flag, MATE_UNMAPPED_FLAG);
    }

    /**
     * strand of the query (false for forward; true for reverse strand).
     */
    @Override
    public void setReadNegativeStrandFlag(final boolean flag) {
        setFlag(flag, READ_STRAND_FLAG);
    }

    /**
     * strand of the mate (false for forward; true for reverse strand).
     */
    @Override
    public void setMateNegativeStrandFlag(final boolean flag) {
        setFlag(flag, MATE_STRAND_FLAG);
    }

    /**
     * the read is the first read in a pair.
     */
    @Override
    public void setFirstOfPairFlag(final boolean flag) {
        setFlag(flag, FIRST_OF_PAIR_FLAG);
    }

    /**
     * the read is the second read in a pair.
     */
    @Override
    public void setSecondOfPairFlag(final boolean flag) {
        setFlag(flag, SECOND_OF_PAIR_FLAG);
    }

    /**
     * the alignment is not primary (a read having split hits may have multiple primary alignment records).
     */
    @Override
    public void setNotPrimaryAlignmentFlag(final boolean flag) {
        setFlag(flag, NOT_PRIMARY_ALIGNMENT_FLAG);
    }

    /**
     * the alignment is supplementary (TODO: further explanation?).
     */
    @Override
    public void setSupplementaryAlignmentFlag(final boolean flag) {
        setFlag(flag, SUPPLEMENTARY_ALIGNMENT_FLAG);
    }

    /**
     * the read fails platform/vendor quality checks.
     */
    @Override
    public void setReadFailsVendorQualityCheckFlag(final boolean flag) {
        setFlag(flag, READ_FAILS_VENDOR_QUALITY_CHECK_FLAG);
    }

    /**
     * the read is either a PCR duplicate or an optical duplicate.
     */
    @Override
    public void setDuplicateReadFlag(final boolean flag) {
        setFlag(flag, DUPLICATE_READ_FLAG);
    }

    private void setFlag(final boolean flag, final int bit) {
        if (flag) {
            mFlags |= bit;
        } else {
            mFlags &= ~bit;
        }
    }

    /**
     * Tests if this record is either a secondary and/or supplementary alignment;
     * equivalent to {@code (getNotPrimaryAlignmentFlag() || getSupplementaryAlignmentFlag())}.
     */
    @Override
    public boolean isSecondaryOrSupplementary() {
        return getNotPrimaryAlignmentFlag() || getSupplementaryAlignmentFlag();
    }


}
