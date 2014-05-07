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

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class CoordinateSortedIndexWrapper implements Index {
    private final CoordinateSortedIndex index;
    private final List<String> sequenceNames;

    public CoordinateSortedIndexWrapper(final CoordinateSortedIndex index, final List<String> sequenceNames) {
        this.index = index;
        this.sequenceNames = Collections.unmodifiableList(sequenceNames);
        if (sequenceNames.size() != index.getNumSequences()) {
            throw new TribbleException(String.format("Length of sequence dictionary %d != number of references in index %d",
                    sequenceNames.size(), index.getNumSequences()));
        }
    }


    @Override
    public List<Block> getBlocks(final String chr, final int start, final int end) {
        final int sequenceIndex = sequenceNames.indexOf(chr);
        if (sequenceIndex == -1) return Collections.emptyList();
        return index.getBlocks(sequenceIndex, start, end);
    }

    @Override
    public boolean isCurrentVersion() {
        return index.isCurrentVersion();
    }

    @Override
    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    @Override
    public boolean containsChromosome(final String chr) {
        return sequenceNames.contains(chr);
    }

    /**
     * @param stream the stream to write the index to.  Caller must close after invocation.  This should be a
     *               BlockCompressedOutputStream, but does not cause an error if it is not.
     * @throws IOException
     */
    @Override
    public void write(final LittleEndianOutputStream stream) throws IOException {
        index.write(stream);
    }

    @Override
    public void writeBasedOnFeatureFile(final File featureFile) throws IOException {
        index.writeBasedOnFeatureFile(featureFile);
    }

    @Override
    public Map<String, String> getProperties() {
        return null;
    }

    @Override
    public boolean equalsIgnoreProperties(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final CoordinateSortedIndexWrapper that = (CoordinateSortedIndexWrapper) o;

        if (!index.equals(that.index)) return false;
        if (!sequenceNames.equals(that.sequenceNames)) return false;

        return true;
    }
}
