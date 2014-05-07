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

import htsjdk.samtools.Bin;
import htsjdk.samtools.BinningIndexContent;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.LinearIndex;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.util.LittleEndianInputStream;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

public class CoordinateSortedIndex {
    public static final String CSI_EXTENSION = ".csi";

    private static final byte[] MAGIC = {'C', 'S', 'I', 1};
    public static final int MAGIC_NUMBER;
    static {
        final ByteBuffer bb = ByteBuffer.allocate(MAGIC.length);
        bb.put(MAGIC);
        bb.flip();
        MAGIC_NUMBER = bb.order(ByteOrder.LITTLE_ENDIAN).getInt();
    }

    private final int min_shift;
    private final int depth;
    private final BinningIndexContent[] indices;
    private final byte[] auxiliaryBlock;


    public CoordinateSortedIndex(final int min_shift, final int depth, final BinningIndexContent[] indices, final byte[] auxiliaryBlock) {
        this.min_shift = min_shift;
        this.depth = depth;
        this.indices = indices;
        this.auxiliaryBlock = auxiliaryBlock;
    }

    /**
     * @param inputStream This is expected to be buffered and be gzip-decompressing as appropriate.  Caller
     *                    should close input stream after ctor returns.
     */
    public CoordinateSortedIndex(final InputStream inputStream) throws IOException {
        this(inputStream, false);
    }

    /**
     * Convenient ctor that opens the file, wraps with with BGZF reader, and closes after reading index.
     */
    public CoordinateSortedIndex(final File csiFile) throws IOException {
        this(new BlockCompressedInputStream(csiFile), true);
    }

    private CoordinateSortedIndex(final InputStream inputStream, final boolean closeInputStream) throws IOException {
        final LittleEndianInputStream dis = new LittleEndianInputStream(inputStream);
        if (dis.readInt() != MAGIC_NUMBER) {
            throw new TribbleException(String.format("Unexpected magic number 0x%x", MAGIC_NUMBER));
        }
        min_shift = dis.readInt();
        depth = dis.readInt();
        final int auxiliaryLength = dis.readInt();
        if (auxiliaryLength > 0) {
            auxiliaryBlock = new byte[auxiliaryLength];
            final int bytesRead = dis.read(auxiliaryBlock);
            if (bytesRead != auxiliaryLength) {
                throw new IOException(String.format("Expected to read %d bytes, but only read %d", auxiliaryLength, bytesRead));
            }
        } else {
            auxiliaryBlock = null;
        }
        final int numSequences = dis.readInt();
        indices = new BinningIndexContent[numSequences];

        for (int i = 0; i < numSequences; ++i) {
            indices[i] = loadSequence(i, dis);
        }
        if (closeInputStream) CloserUtil.close(dis);
    }

    public int getNumSequences() {
        return indices.length;
    }

    private BinningIndexContent loadSequence(final int referenceSequenceIndex, final LittleEndianInputStream dis) throws IOException {
        final int numBins = dis.readInt();
        if (numBins == 0) return null;
        int nonNullBins = 0;
        final ArrayList<Bin> bins = new ArrayList<Bin>();
        final SortedMap<Integer, Long> linearEntries = new TreeMap<Integer, Long>();
        for (int i = 0; i < numBins; ++i) {
            final Map.Entry<Bin, Long> entry = loadBin(referenceSequenceIndex, dis);
            final Bin bin = entry.getKey();
            final int linearIndexOffset = GenomicIndexUtil.getLinearIndexOffsetForBin(bin.getBinNumber(), depth);
            final Long existingLinearIndexValue = linearEntries.get(linearIndexOffset);
            if (existingLinearIndexValue != null) {
                if (existingLinearIndexValue != entry.getValue()) {
                    throw new TribbleException(String.format("Linear index mismatch for bin %d, offset %d: %ld != %ld",
                            bin.getBinNumber(), linearIndexOffset, existingLinearIndexValue, entry.getValue()));
                }
            } else {
                linearEntries.put(linearIndexOffset, entry.getValue());
            }
            // File is not sparse, but array being produced is sparse, so grow array with nulls as appropriate
            // so that bin number == index into array.
            ++nonNullBins;
            if (bins.size() > bin.getBinNumber()) {
                if (bins.get(bin.getBinNumber()) != null) {
                    throw new TribbleException("Bin " + bin.getBinNumber() + " appears more than once in file");
                }
                bins.set(bin.getBinNumber(),bin);
            } else {
                // Grow bins array as needed.
                bins.ensureCapacity(bin.getBinNumber() + 1);
                while (bins.size() < bin.getBinNumber()) bins.add(null);
                bins.add(bin);
            }
        }
        final int firstLinearIndexOffset = linearEntries.firstKey();
        final long[] linearIndexArray = new long[linearEntries.lastKey() - firstLinearIndexOffset + 1];
        for (final Map.Entry<Integer, Long> linearEntry : linearEntries.entrySet()) {
            linearIndexArray[linearEntry.getKey() - firstLinearIndexOffset] = linearEntry.getValue();
        }

        return new BinningIndexContent(referenceSequenceIndex,
                new BinningIndexContent.BinList(bins.toArray(new Bin[bins.size()]), nonNullBins),
                new LinearIndex(referenceSequenceIndex, firstLinearIndexOffset, linearIndexArray));
    }

    private Map.Entry<Bin, Long> loadBin(final int referenceSequenceIndex, final LittleEndianInputStream dis) throws IOException {
        final int binNumber = dis.readInt();
        final long linearOffset = dis.readLong();
        final Bin bin = new Bin(referenceSequenceIndex, binNumber);
        final int numChunks = dis.readInt();
        final List<Chunk> chunkList = new ArrayList<Chunk>(numChunks);
        for (int i = 0; i < numChunks; ++i) {
            chunkList.add(loadChunk(dis));
        }
        bin.setChunkList(chunkList);
        return new Map.Entry<Bin, Long>() {
            @Override
            public Bin getKey() {
                return bin;
            }

            @Override
            public Long getValue() {
                return linearOffset;
            }

            @Override
            public Long setValue(final Long aLong) {
                return null;
            }
        };
    }

    private Chunk loadChunk(final LittleEndianInputStream dis) throws IOException {
        final long start = dis.readLong();
        final long end = dis.readLong();
        return new Chunk(start, end);
    }

    public List<Block> getBlocks(final int sequenceIndex, final int start, final int end) {
        if (sequenceIndex == -1 || indices[sequenceIndex] == null) {
            return Collections.emptyList();
        }
        return Block.makeBlockListFromChunkList(indices[sequenceIndex].getChunksOverlapping(start, end));
    }

    public boolean isCurrentVersion() {
        return true;
    }

    /**
     * @param los the stream to write the index to.  Caller must close after invocation.  This should be a
     *               BlockCompressedOutputStream, but does not cause an error if it is not.
     * @throws IOException
     */
    public void write(final LittleEndianOutputStream los) throws IOException {
        los.writeInt(MAGIC_NUMBER);
        los.writeInt(min_shift);
        los.writeInt(depth);
        if (auxiliaryBlock == null) {
            los.writeInt(0);
        } else {
            los.writeInt(auxiliaryBlock.length);
            los.write(auxiliaryBlock);
        }
        los.writeInt(indices.length);
        for (final BinningIndexContent index : indices) {
            writeSequence(index, los);
        }

    }

    private void writeSequence(final BinningIndexContent indexContent, final LittleEndianOutputStream los) throws IOException {
        if (indexContent == null) {
            los.writeInt(0);
        } else {
            final BinningIndexContent.BinList binList = indexContent.getBins();
            los.writeInt(binList.numberOfNonNullBins);
            for (final Bin bin : binList) {
                final long linearIndexEntry = indexContent.getLinearIndex().get(GenomicIndexUtil.getLinearIndexOffsetForBin(bin.getBinNumber(), depth));
                writeBin(bin, linearIndexEntry, los);
            }
        }
    }

    private void writeBin(final Bin bin, final long linearIndexEntry, final LittleEndianOutputStream los) throws IOException {
        los.writeInt(bin.getBinNumber());
        los.writeLong(linearIndexEntry);
        final List<Chunk> chunkList = bin.getChunkList();
        los.writeInt(chunkList.size());
        for (final Chunk chunk: chunkList) {
            los.writeLong(chunk.getChunkStart());
            los.writeLong(chunk.getChunkEnd());
        }

    }

    /**
     * Writes the index with BGZF.
     * @param indexFile Where to write the index.
     */
    public void write(final File indexFile) throws IOException {
        final LittleEndianOutputStream los = new LittleEndianOutputStream(new BlockCompressedOutputStream(indexFile));
        write(los);
        los.close();
    }



    public void writeBasedOnFeatureFile(final File featureFile) throws IOException {
        if (!featureFile.isFile()) return;
        write(new File(featureFile.getAbsolutePath() + CSI_EXTENSION));

    }

    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final CoordinateSortedIndex that = (CoordinateSortedIndex) o;

        if (min_shift != that.min_shift) return false;
        if (depth != that.depth) return false;
        if (!Arrays.equals(indices, that.indices)) return false;
        return Arrays.equals(auxiliaryBlock, that.auxiliaryBlock);
    }

    @Override
    public int hashCode() {
        int result = min_shift;
        result = 31 * result + depth;
        result = 31 * result + Arrays.hashCode(indices);
        result = 31 * result + (auxiliaryBlock != null ? Arrays.hashCode(auxiliaryBlock) : 0);
        return result;
    }
}
