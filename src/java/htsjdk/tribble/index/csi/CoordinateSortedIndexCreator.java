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

import htsjdk.samtools.BinningIndexContent;
import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.index.BinningIndexCreator;

public class CoordinateSortedIndexCreator extends BinningIndexCreator<CoordinateSortedIndexWrapper> {
    public static final int DEFAULT_MIN_SHIFT = 14;
    private final int min_shift;
    private byte[] auxiliaryBlock = null;

    public CoordinateSortedIndexCreator(final int min_shift, final SAMSequenceDictionary sequenceDictionary) {
        super(sequenceDictionary);
        this.min_shift = min_shift;
    }

    public byte[] getAuxiliaryBlock() {
        return auxiliaryBlock;
    }

    public void setAuxiliaryBlock(final byte[] auxiliaryBlock) {
        this.auxiliaryBlock = auxiliaryBlock;
    }

    @Override
    protected CoordinateSortedIndexWrapper createIndex(final BinningIndexContent[] indices) {
        // TODO: compute max length
        int maxSequenceLength = 0;
        if (getSequenceDictionary() != null) {
            for (final SAMSequenceRecord sequenceRecord : getSequenceDictionary().getSequences()) {
                if (sequenceRecord.getSequenceLength() > maxSequenceLength) {
                    maxSequenceLength = sequenceRecord.getSequenceLength();
                }
            }
        }
        return new CoordinateSortedIndexWrapper(
                new CoordinateSortedIndex(min_shift, GenomicIndexUtil.computeDepth(min_shift, maxSequenceLength), indices, auxiliaryBlock),
                getSequenceNames());
    }
}
