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

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

/**
 * @author alecw@broadinstitute.org
 */
public class SAMRecordUtil {

    /** List of String tags that must be reversed if present when a SAMRecord is reverseComplemented */
    private static final short[] STRING_TAGS_TO_REVERSE = {
            SAMTagUtil.U2,
            SAMTagUtil.OQ
    };

    /**
     * Reverse-complement all known sequence and base quality attributes of the SAMRecord.
     */
    public static void reverseComplement(final ReadRecord readRec) {
        final FastBAMRecord rec = (FastBAMRecord) readRec;
        rec.reverseComplementReadBases();
        rec.reverseBaseQualities();

        if(rec.hasAttribute(SAMTagUtil.getSingleton().SQ))
            throw new IllegalArgumentException("SQ tag no longer supported(?)");
            //SQTagUtil.reverseComplementSqArray(sqTagValue);
            //rec.setAttribute(SAMTagUtil.getSingleton().SQ, sqTagValue);

        if(rec.hasAttribute(SAMTagUtil.getSingleton().E2))
            rec.reverseOrReverseComplementStringAttribute(SAMTagUtil.getSingleton().E2, true);

        for (final short tagToReverse : STRING_TAGS_TO_REVERSE)
            if(rec.hasAttribute(tagToReverse))
                rec.reverseOrReverseComplementStringAttribute(tagToReverse, false);
    }
}
