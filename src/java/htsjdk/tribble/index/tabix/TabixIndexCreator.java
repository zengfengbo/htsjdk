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
package htsjdk.tribble.index.tabix;

import htsjdk.samtools.BinningIndexBuilder;
import htsjdk.samtools.BinningIndexContent;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * IndexCreator for Tabix.
 * Features are expected to be 1-based, inclusive.
 */
public class TabixIndexCreator extends htsjdk.tribble.index.BinningIndexCreator<TabixIndex> {
    private final TabixFormat formatSpec;

    /**
     * @param sequenceDictionary is not required, but if present all features added must refer to sequences in the
     *                           dictionary.  It is used to optimize the memory needed to build the index.
     */
    public TabixIndexCreator(final SAMSequenceDictionary sequenceDictionary,
                             final TabixFormat formatSpec) {
        super(sequenceDictionary);
        this.formatSpec = formatSpec.clone();
    }

    public TabixIndexCreator(final TabixFormat formatSpec) {
        this(null, formatSpec);
    }

    @Override
    protected TabixIndex createIndex(final BinningIndexContent[] indices) {
        return new TabixIndex(formatSpec, getSequenceNames(), indices);
    }
}
