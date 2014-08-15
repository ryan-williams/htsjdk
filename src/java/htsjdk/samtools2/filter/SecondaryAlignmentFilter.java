package htsjdk.samtools2.filter;

import htsjdk.samtools2.ReadRecord;

/**
 * SamRecordFilter that filters out secondary alignments, but not supplemental alignments.
 */
public class SecondaryAlignmentFilter implements SamRecordFilter {
    /**
     * Returns true if the read is marked as secondary.
     */
    public boolean filterOut(final ReadRecord record) { return record.getNotPrimaryAlignmentFlag(); }

    /**
     * Returns true if either read is marked as secondary.
     */
    public boolean filterOut(final ReadRecord first, final ReadRecord second) {
        return first.getNotPrimaryAlignmentFlag() || second.getNotPrimaryAlignmentFlag();
    }
}
