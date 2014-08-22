package htsjdk.samtools;

/**
 * Factory class for generating instances of SAMRecord and BAMRecord when reading from SAM/BAM files.
 *
 * @author Tim Fennell
 */
public class SAMRecordFactory {

    private static final SAMRecordFactory INSTANCE = new SAMRecordFactory();

    public static SAMRecordFactory getInstance() {
        return INSTANCE;
    }

    /** Create a new SAMRecord to be filled in */
    public ReadRecord createSAMRecord(final SAMFileHeader header) {
        return new SAMRecord(header);
    }

    /** Create a new BAM Record. */
    public ReadRecord createBAMRecord (final SAMFileHeader header,
                                      final int referenceSequenceIndex,
                                      final int alignmentStart,
                                      final short readNameLength,
                                      final short mappingQuality,
                                      final int indexingBin,
                                      final int cigarLen,
                                      final int flags,
                                      final int readLen,
                                      final int mateReferenceSequenceIndex,
                                      final int mateAlignmentStart,
                                      final int insertSize,
                                      final byte[] variableLengthBlock) {

        return new BAMRecord(header,
                referenceSequenceIndex,
                alignmentStart,
                readNameLength,
                mappingQuality,
                indexingBin,
                cigarLen,
                flags,
                readLen,
                mateReferenceSequenceIndex,
                mateAlignmentStart,
                insertSize,
                variableLengthBlock);
    }


    /** Create a new BAM Record. */
    public ReadRecord createFastBAMRecord (final SAMFileHeader header,
                                           final byte[] record) {
        return new FastBAMRecord(header, record);
    }
}
