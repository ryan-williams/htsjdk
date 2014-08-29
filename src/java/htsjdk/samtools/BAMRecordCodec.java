/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package htsjdk.samtools;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.SortingCollection;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

/**
 * Class for translating between in-memory and disk representation of BAMRecord.
 */
public class BAMRecordCodec implements SortingCollection.Codec<ReadRecord> {
    private final BinaryCigarCodec cigarCodec = new BinaryCigarCodec();
    private final SAMFileHeader header;
    private final BinaryCodec binaryCodec = new BinaryCodec();
    private final BinaryTagCodec binaryTagCodec = new BinaryTagCodec(binaryCodec);
    private final SAMRecordFactory samRecordFactory;

    public BAMRecordCodec(final SAMFileHeader header) {
        this(header, SAMRecordFactory.getInstance());
    }

    public BAMRecordCodec(final SAMFileHeader header, final SAMRecordFactory factory) {
        this.header = header;
        this.samRecordFactory = factory;
    }

    public BAMRecordCodec clone() {
        // Do not clone the references to codecs, as they must be distinct for each instance.
        return new BAMRecordCodec(this.header, this.samRecordFactory);
    }


    /** Sets the output stream that records will be written to. */
    public void setOutputStream(final OutputStream os) {
        this.binaryCodec.setOutputStream(os);
    }

    /** Sets the output stream that records will be written to. */
    public void setOutputStream(final OutputStream os, final String filename) {
        this.binaryCodec.setOutputStream(os);
        this.binaryCodec.setOutputFileName(filename);
    }

    /** Sets the input stream that records will be read from. */
    public void setInputStream(final InputStream is) {
        this.binaryCodec.setInputStream(is);
    }

    /** Sets the input stream that records will be read from. */
    public void setInputStream(final InputStream is, final String filename) {
        this.binaryCodec.setInputStream(is);
        this.binaryCodec.setInputFileName(filename);
    }

    /**
     * Write object to OutputStream.
     * Reference and mate reference indices must be resolvable, which either means that these have been set into the
     * SAMRecord directly, or the SAMRecord must have a header assigned into it so that reference names can be
     * resolved into indices.
     *
     * @param alignment Record to be written.
     */
    public void encode(final ReadRecord alignment) {
        //TODO make sure index bin logic here is correct
        if(alignment.getIndexingBin() == null && alignment.getReferenceIndex() >= 0)
            alignment.setIndexingBin(alignment.computeIndexingBin());
        else
            alignment.setIndexingBin(0);

        byte[] record = ((FastBAMRecord) alignment).getRecord();
        this.binaryCodec.writeInt(record.length);
        this.binaryCodec.writeBytes(record);
    }

    /**
     * Read the next record from the input stream and convert into a java object.
     *
     * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
     *         a record.
     */
    public ReadRecord decode() {
        int recordLength = 0;
        try {
            recordLength = this.binaryCodec.readInt();
        }
        catch (RuntimeEOFException e) {
            return null;
        }

        if (recordLength < BAMFileConstants.FIXED_BLOCK_SIZE) {
            throw new SAMFormatException("Invalid record length: " + recordLength);
        }

        final byte[] record = new byte[recordLength];
        this.binaryCodec.readBytes(record);
        return this.samRecordFactory.createFastBAMRecord(header, record);

    }
}
