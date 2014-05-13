/*******************************************************************************
 * Copyright 2013 EMBL-EBI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package htsjdk.samtools;

import htsjdk.samtools.cram.build.ContainerFactory;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.build.Sam2CramRecordFactory;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CramCompressionRecord;
import htsjdk.samtools.cram.structure.CramHeader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

public class CRAMFileWriter extends SAMFileWriterImpl {
    private final String fileName;
    private final List<CramCompressionRecord> records = new ArrayList<CramCompressionRecord>();
    private ContainerFactory containerFactory;
    private final int recordsPerSlice = 10000;

    private Sam2CramRecordFactory sam2CramRecordFactory;
    private final OutputStream os;

    private boolean preserveReadNames = false;

    public CRAMFileWriter(final File outputFile) throws FileNotFoundException {
        fileName = outputFile.getName();
        os = new FileOutputStream(outputFile);
    }

    boolean shouldFlushContainer() {
        final int containerSize = recordsPerSlice * 10;
        return records.size() >= containerSize;
    }

    void flushContainer() throws IllegalArgumentException,
            IllegalAccessException, IOException {
        final Container container = containerFactory.buildContainer(records);
        records.clear();
        CramIO.writeContainer(container, os);
    }

    @Override
    protected void writeAlignment(final SAMRecord alignment) {
        if (shouldFlushContainer())
            try {
                flushContainer();
            } catch (final Exception e) {
                throw new RuntimeException(e);
            }

        final CramCompressionRecord cramRecord = sam2CramRecordFactory
                .createCramRecord(alignment);
        records.add(cramRecord);
    }

    @Override
    protected void writeHeader(final String textHeader) {
        containerFactory = new ContainerFactory(getFileHeader(), recordsPerSlice,
                preserveReadNames);

        sam2CramRecordFactory = new Sam2CramRecordFactory(new byte[0], getFileHeader());

        final CramHeader cramHeader = new CramHeader(2, 0, fileName, getFileHeader());

        try {
            CramIO.writeCramHeader(cramHeader, os);
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    protected void finish() {
        if (!records.isEmpty())
            try {
                flushContainer();
            } catch (final Exception e) {
                throw new RuntimeException(e);
            }
    }

    @Override
    protected String getFilename() {
        return fileName;
    }

    public boolean isPreserveReadNames() {
        return preserveReadNames;
    }

    public void setPreserveReadNames(final boolean preserveReadNames) {
        this.preserveReadNames = preserveReadNames;
    }
}
