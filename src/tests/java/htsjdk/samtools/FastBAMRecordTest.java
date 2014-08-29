package htsjdk.samtools;

import htsjdk.samtools.util.SequenceUtil;

import static org.testng.Assert.*;
import static org.testng.Assert.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.ArrayList;

import htsjdk.samtools.util.StringUtil;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class FastBAMRecordTest {
    private final File BASIC_BAM_FILE = new File("testdata/htsjdk/samtools/FastBAMRecordTest/basic.bam");

    private ArrayList<FastBAMRecord> records = new ArrayList<FastBAMRecord>();

    private FastBAMRecord testRecord = null;


    @BeforeMethod(alwaysRun = true)
    public void initRecords() throws Exception {
        records.clear();
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
        final SAMFileReader reader = new SAMFileReader(BASIC_BAM_FILE);
        for (ReadRecord record : reader) {
            records.add((FastBAMRecord) record);
        }
        reader.close();

        testRecord = (FastBAMRecord) records.get(records.size() - 1).clone(); //use the last record as the test record
    }

    @Test
    public void testBasicGettersAndSetters() {
        //read name
        testRecord.setReadName("bla");
        assertEquals(testRecord.getReadName(), "bla");
        testRecord.setReadName("bla2");
        assertEquals(testRecord.getReadName(), "bla2");
        assertEquals(testRecord.getReadNameLength(), 4);

        testRecord.setFlags(14);
        assertEquals(testRecord.getFlags(), 14);

        testRecord.setMappingQuality(0x80);
        assertEquals(testRecord.getMappingQuality(), 0x80);

        //set read bases
        byte[] temp = new byte[] {'A', 'C', 'G'};
        testRecord.setReadBases(temp);
        assertEquals(testRecord.getReadLength(), 3);
        assertEquals(testRecord.getReadString(), "ACG");
        assertEquals(testRecord.getReadBases(), temp);
    }

    @Test
    public void testCigar() throws CloneNotSupportedException {
        FastBAMRecord r1 = records.get(0);
        assertEquals(r1.getCigarString(), "45M56S");
        assertEquals(r1.getCigarLength(), 2);
        assertEquals(r1.getCigarReferenceLength(), 45);
        assertEquals(r1.getCigarOps(), new CigarOperator[] { CigarOperator.M, CigarOperator.S });
        assertEquals(r1.getCigarOp(0), CigarOperator.M);
        assertEquals(r1.getCigarOp(1), CigarOperator.S);
        assertEquals(r1.getCigarOpLength(0), 45);
        assertEquals(r1.getCigarOpLength(1), 56);

        FastBAMRecord r2 = records.get(1);
        assertEquals(r2.getCigarString(), "*");
        assertEquals(r2.getCigarLength(), 0);

        //test copyCigarFrom(..)
        FastBAMRecord r3 = r2.clone();
        r3.copyCigarFrom(r1);
        assertEquals(r3.getCigarString(), "45M56S");
        assertEquals(r3.getCigarLength(), 2);
        assertEquals(r3.getCigarReferenceLength(), 45);
        assertEquals(r3.getCigarOps(), new CigarOperator[] { CigarOperator.M, CigarOperator.S });
        assertEquals(r3.getCigarOp(0), CigarOperator.M);
        assertEquals(r3.getCigarOp(1), CigarOperator.S);
        assertEquals(r3.getCigarOpLength(0), 45);
        assertEquals(r3.getCigarOpLength(1), 56);

        r3.copyCigarFrom(r2);
        assertEquals(r3.getCigarString(), "*");
        assertEquals(r3.getCigarLength(), 0);
        assertEquals(r3.getCigarReferenceLength(), 0);
    }

    @Test
    public void testAttributes() throws CloneNotSupportedException {
        FastBAMRecord r = testRecord;
        r.clearAttributes();
        assertEquals(r.getNumAttributes(), 0);

        final short tag1 = SAMTagUtil.makeBinaryTag("XA");
        final short tag2 = SAMTagUtil.makeBinaryTag("XB");
        final short tag3 = SAMTagUtil.makeBinaryTag("XC");
        final short tag4 = SAMTagUtil.makeBinaryTag("XD");
        final short tag5 = SAMTagUtil.makeBinaryTag("XE");

        //set 5 tags of all types
        r.setAttribute(tag1, 'c', 0x3);
        r.setAttribute(tag2, 's', 0x3);
        r.setAttribute(tag3, 'i', 0x3);
        r.setAttribute(tag4, 'Z', "olleh");
        r.setAttribute(tag5, 'c', 0x3);

        //check values
        assertEquals(r.getNumAttributes(), 5);
        assertEquals(r.getAttribute(tag1), (byte) 3);
        assertEquals(r.getAttribute(tag2), (short) 3);
        assertEquals(r.getAttribute(tag3), 3);
        assertEquals(r.getAttribute(tag4), "olleh");
        assertEquals(r.getAttribute(tag5), (byte) 3);

        //change 5 tags of all types
        r.setAttribute(tag1, 'c', 0x8F);
        r.setAttribute(tag2, 's', 0x8FFF);
        r.setAttribute(tag3, 'i', 0x8FFFFFFF);
        r.setAttribute(tag4, 'Z', "hello");
        r.setAttribute(tag5, 'c', 0x1);

        //check values
        assertEquals(r.getNumAttributes(), 5);
        assertEquals(r.getAttribute(tag1), (byte) 0x8F);
        assertEquals(r.getAttribute(tag2), (short) 0x8FFF);
        assertEquals(r.getAttribute(tag3), 0x8FFFFFFF);
        assertEquals(r.getAttribute(tag4), "hello");
        assertEquals(r.getAttribute(tag5), (byte) 1);

        //delete 1 tag
        r.deleteAttribute(tag1);
        assertEquals(r.getNumAttributes(), 4);
        assertEquals(r.getAttribute(tag2), (short) 0x8FFF);
        assertEquals(r.getAttribute(tag3), 0x8FFFFFFF);
        assertEquals(r.getAttribute(tag4), "hello");
        assertEquals(r.getAttribute(tag5), (byte) 1);

        //test hasAttribute(..)
        assertFalse(r.hasAttribute(tag1));
        assertTrue(r.hasAttribute(tag2));

        //delete 1 more tag
        r.deleteAttribute(tag4);
        assertEquals(r.getNumAttributes(), 3);
        assertEquals(r.getAttribute(tag2), (short) 0x8FFF);
        assertEquals(r.getAttribute(tag3), 0x8FFFFFFF);
        assertEquals(r.getAttribute(tag5), (byte) 1);

        assertFalse(r.hasAttribute(tag4));
        assertTrue(r.hasAttribute(tag3));

        //delete 1 more tag
        r.deleteAttribute(tag3);
        assertEquals(r.getNumAttributes(), 2);
        assertEquals(r.getAttribute(tag2), (short) 0x8FFF);
        assertEquals(r.getAttribute(tag5), (byte) 1);

        //add the tag back
        r.setAttribute(tag3, 1);
        assertEquals(r.getAttribute(tag3), (byte) 1);

        //change the value
        r.setAttribute(tag3, 2);
        assertEquals(r.getAttribute(tag3), (byte) 2);

        assertEquals(r.getAttributeTags(), new Short[] { tag2, tag5, tag3 });

        //test reverseOrReverseComplementStringAttribute(..)
        r.setAttribute(tag4, "ACGTCAGC");
        assertEquals(r.getAttribute(tag4), "ACGTCAGC");
        r.reverseOrReverseComplementStringAttribute(tag4, true); //reverse complement
        assertEquals(r.getAttribute(tag4), "GCTGACGT");
        r.reverseOrReverseComplementStringAttribute(tag4, true); //reverse complement it back
        assertEquals(r.getAttribute(tag4), "ACGTCAGC");

        r.reverseOrReverseComplementStringAttribute(tag4, false); //just reverse it
        assertEquals(r.getAttribute(tag4), "CGACTGCA");
        r.reverseOrReverseComplementStringAttribute(tag4, false); //just reverse it
        assertEquals(r.getAttribute(tag4), "ACGTCAGC");

        //test copyAttributeFrom(..)
        FastBAMRecord r2 = r.clone();
        r.setAttribute(tag4, "A");
        assertEquals(r.getAttribute(tag4), "A");
        assertEquals(r2.getAttribute(tag4), "ACGTCAGC");

        r2.copyAttributeFrom(r, tag4);
        assertEquals(r2.getAttribute(tag4), "A");
    }


    //TODO test getAlignmentBlocks()?

    @Test
    public void testEncodeDecode() {
        BAMRecordCodec codec = new BAMRecordCodec(testRecord.getHeader());
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        codec.setOutputStream(os);
        for(ReadRecord record : records) {
            codec.encode(record);
        }
        testRecord.setReferenceIndex(-1);
        testRecord.setAlignmentStart(-1);
        codec.encode(testRecord);

        codec.setInputStream(new ByteArrayInputStream(os.toByteArray()));
        for(ReadRecord record : records) {
            assertEquals(codec.decode().toString(), record.toString());
        }
        assertEquals(codec.decode().toString(), testRecord.toString());
    }


    private static byte encodeNibbles(String bases) {
        byte[] compressedBytes = SAMUtils.bytesToCompressedBases(new byte[] { (byte) bases.charAt(0), (byte) bases.charAt(1) });
        return compressedBytes[0];
    }

    private static String decodeNibbles(byte b) {
        byte[] bases = SAMUtils.compressedBasesToBytes(2, new byte[] { b }, 0);
        return (char) bases[0] + "" + (char) bases[1];
    }

    @Test
    public void testSwapAndComplementNibbles() {
        assertEquals(decodeNibbles(encodeNibbles("CG")), "CG");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("CT"))), "AG");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("AT"))), "AT");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("CG"))), "CG");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("AG"))), "CT");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("TT"))), "AA");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("GG"))), "CC");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("CC"))), "GG");
        assertEquals(decodeNibbles(FastBAMRecord.swapAndComplementNibbles(encodeNibbles("AA"))), "TT");
    }

    @Test
    public void testReverseComplementReadBases() {
        testRecord.setReadBases("T".getBytes());
        assertEquals(testRecord.getReadString(), "T");
        testRecord.reverseComplementReadBases();
        assertEquals(testRecord.getReadString(), "A");

        testRecord.setReadBases("AG".getBytes());
        assertEquals(testRecord.getReadString(), "AG");
        testRecord.reverseComplementReadBases();
        assertEquals(testRecord.getReadString(), "CT");

        testRecord.setReadBases("AGC".getBytes());
        assertEquals(testRecord.getReadString(), "AGC");
        testRecord.reverseComplementReadBases();
        assertEquals(testRecord.getReadString(), "GCT");

        testRecord.setReadBases("CGCT".getBytes());
        assertEquals(testRecord.getReadString(), "CGCT");
        testRecord.reverseComplementReadBases();
        assertEquals(testRecord.getReadString(), "AGCG");


        testRecord.setReadBases("AGCTA".getBytes());
        assertEquals(testRecord.getReadString(), "AGCTA");
        testRecord.reverseComplementReadBases();
        assertEquals(testRecord.getReadString(), "TAGCT");
    }




    @Test(dataProvider = "testBasicGetters")
    public void testBasicGettersOnData(int recordNumber,
                                       String readName,
                                       int flags,
                                       String referenceSequence,
                                       int position,
                                       final int mapq,
                                       String cigar,
                                       String seq,
                                       String qual) {

        FastBAMRecord r = records.get(recordNumber);
        assertEquals(r.getReadName(), readName);
        assertEquals(r.getReadNameLength(), readName.length());
        assertEquals(r.getFlags(), flags);
        assertEquals(r.getMappingQuality(), mapq);
        assertEquals(r.getReadLength(), seq.length());
        assertEquals(r.getAlignmentStart(), position);
        assertEquals(r.getMappingQuality(), mapq);
        assertEquals(r.getCigarString(), cigar);
        assertEquals(r.getReadString(), seq);
        assertEquals(r.getBaseQualityString(), qual);
    }


    @Test(dataProvider = "testBasicGetters")
    public void testReverseComplementOnData(int recordNumber,
                                             String readName,
                                             int flags,
                                             String referenceSequence,
                                             int position,
                                             final int mapq,
                                             String cigar,
                                             String seq,
                                             String qual) {

        FastBAMRecord r = records.get(recordNumber);

        assertEquals(r.getReadString(), seq);
        r.reverseComplementReadBases();
        assertEquals(r.getReadString(), SequenceUtil.reverseComplement(seq));

        assertEquals(r.getBaseQualityString(), qual);
        r.reverseBaseQualities();
        assertEquals(r.getBaseQualityString(), StringUtil.reverseString(qual));

    }

    @DataProvider(name = "testBasicGetters")
    public Object[][] testBasicGettersDataProvider() {
        return new Object[][] {
            {0, "B084HABXX110425:3:1101:8656:141570", 99, "1",  9997, 0, "45M56S", "CCTAAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAAACCCCTAACCCCTAACTCTAACCCTAACCTTCTCACACAC", ">AB>@BDAABC?><6BDBDACCE@DBA<DCCA24CBDBCCECD?#########################################################"},
            {1, "B084HABXX110425:3:2106:5540:16844", 101, "1", 10002, 0,      "*", "GTTCCAGTTTACTGCTATTTTCGCAGGTTTTGAGTTTTTAACATTGCAAATTACAAAACATTTTTGCCTTAGATAAAAGGAAATTTTAATATTAAAACTCA", "BBACDDDBDCDBECDECDDDDC;CDDDCDDBBEDCDDDDDDBBCBDCC=DCC@BDEDDADCCBDCDCDECBCDCCDEDCCD@DCDCDDDDCBDCDBDADBB"},
        };
    }
}