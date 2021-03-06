package htsjdk.samtools;

/**
 * This class helps us compute and compare duplicate scores, which are used for selecting the non-duplicate
 * during duplicate marking (see MarkDuplicates).
 * @author nhomer
 */
public class DuplicateScoringStrategy {

    public enum ScoringStrategy {
        SUM_OF_BASE_QUALITIES,
        TOTAL_MAPPED_REFERENCE_LENGTH
    }

    /** Calculates a score for the read which is the sum of scores over Q15. */
    private static short getSumOfBaseQualities(final SAMRecord rec) {
        short score = 0;
        for (final byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }

        return score;
    }

    /**
     * Returns the duplicate score computed from the given fragment.
     */
    public static short computeDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy) {
        return computeDuplicateScore(record, scoringStrategy, false);
    }

    /**
     * Returns the duplicate score computed from the given fragment.
     *
     * If true is given to assumeMateCigar, then any score that can use the mate cigar to to compute the mate's score will return the score
     * computed on both ends.
     */
    public static short computeDuplicateScore(final SAMRecord record, final ScoringStrategy scoringStrategy, final boolean assumeMateCigar) {
        short score = 0;

        switch (scoringStrategy) {
            case SUM_OF_BASE_QUALITIES:
                score += getSumOfBaseQualities(record);
                break;
            case TOTAL_MAPPED_REFERENCE_LENGTH:
                if (!record.getReadUnmappedFlag()) {
                    score += record.getCigar().getReferenceLength();
                }
                if (assumeMateCigar && record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                    score += SAMUtils.getMateCigar(record).getReferenceLength();
                }
                break;
        }
        return score;
    }

    /**
     * Compare two records based on their duplicate scores.  If the scores are equal, we break ties based on mapping quality
     * (added to the mate's mapping quality if paired and mapped), then library/read name.
     *
     * If true is given to assumeMateCigar, then any score that can use the mate cigar to to compute the mate's score will return the score
     * computed on both ends.
     *
     * We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2.
     */
    public static int compare(final SAMRecord rec1, final SAMRecord rec2, final ScoringStrategy scoringStrategy, final boolean assumeMateCigar) {
        int cmp;

        // always prefer paired over non-paired
        if (rec1.getReadPairedFlag() != rec2.getReadPairedFlag()) return rec1.getReadPairedFlag() ? 1 : -1;

        cmp = computeDuplicateScore(rec2, scoringStrategy, assumeMateCigar) - computeDuplicateScore(rec1, scoringStrategy, assumeMateCigar);

        /**
         * Finally, use library ID and read name
         * This is important because we cannot control the order in which reads appear for reads that are comparable up to now (i.e. cmp == 0).  We want to deterministically
         * choose them, and so we need this.
         */
        if (0 == cmp) cmp = SAMUtils.getCanonicalRecordName(rec1).compareTo(SAMUtils.getCanonicalRecordName(rec2));

        return cmp;
    }

    /**
     * Compare two records based on their duplicate scores.  The duplicate scores for each record is assume to be
     * pre-computed by computeDuplicateScore and stored in the "DS" tag.  If the scores are equal, we break
     * ties based on mapping quality (added to the mate's mapping quality if paired and mapped), then library/read name.
     *
     * We allow different scoring strategies. We return <0 if rec1 has a better strategy than rec2.
     */
    public static int compare(final SAMRecord rec1, final SAMRecord rec2, final ScoringStrategy scoringStrategy) {
        return compare(rec1, rec2, scoringStrategy, false);
    }

}
