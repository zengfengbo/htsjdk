/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package htsjdk.tribble.index.csi;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.util.LittleEndianOutputStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

public class CoordinateSortedIndexTest {
    private static final File BCF_CSI_FILE = new File("testdata/tribble/csi/b.bcf.gz.csi");
    private static final File VCF_CSI_FILE = new File("testdata/tribble/csi/b.vcf.gz.csi");

    /**
     * Read an existing index from disk, write it to a temp file, read that in, and assert that both in-memory
     * representations are identical.  Disk representations may not be identical due to arbitrary bin order and
     * compression differences.
     */
    @Test(dataProvider = "readWriteTestDataProvider")
    public void readWriteTest(final File csiFile) throws Exception {
        final CoordinateSortedIndex index = new CoordinateSortedIndex(csiFile);
        final File indexFile = File.createTempFile("CoordinateSortedIndexTest.", CoordinateSortedIndex.CSI_EXTENSION);
        final LittleEndianOutputStream los = new LittleEndianOutputStream(new BlockCompressedOutputStream(indexFile));
        index.write(los);
        los.close();
        final CoordinateSortedIndex index2 = new CoordinateSortedIndex(indexFile);
        Assert.assertEquals(index, index2);
        // Unfortunately, can't do byte comparison of original file and temp file, because 1) different compression
        // levels; and more importantly, arbitrary order of bins in bin list.
    }

    @DataProvider(name = "readWriteTestDataProvider")
    public Object[][] readWriteTestDataProvider() {
        return new Object[][] {
                {BCF_CSI_FILE},
                {VCF_CSI_FILE}
        };
    }
}
