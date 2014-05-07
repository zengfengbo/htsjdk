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
package htsjdk.samtools;

import java.util.BitSet;

/**
 * Constants and methods used by BAM and Tribble indices
 */
public class GenomicIndexUtil {
    /**
     * Reports the total amount of genomic data that any bin can index.
     */
    public static final int BIN_GENOMIC_SPAN = 512*1024*1024;

    /**
     * What is the starting bin for each level?
     */
    public static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    /**
     * Reports the maximum number of bins that can appear in a binning index.
     */
    public static final int MAX_BINS = 37450;   // =(8^6-1)/7+1

    public static final int MAX_LINEAR_INDEX_SIZE = MAX_BINS+1-LEVEL_STARTS[LEVEL_STARTS.length-1];

    public static final int BAI_MIN_SHIFT = 14;
    public static final int BAI_DEPTH = 5;


    /**
     * E.g. for a SAMRecord with no genomic coordinate.
     */
    public static final int UNSET_GENOMIC_LOCATION = 0;

    /**
     * Calculate the bin given an alignment covering [beg, end)
     * @param beg 0-based, inclusive
     * @param end 0-based, exclusive
     * @return bin number
     */
    static int reg2bin(final int beg, int end) {
        int ret = reg2bin_original(beg, end);
        if (ret != reg2bin(beg, end, BAI_MIN_SHIFT, BAI_DEPTH)) {
            throw new IllegalStateException(String.format("Results do not agree for beg: %d end: %d", beg, end));
        }
        return ret;
    }

    /**
     * calculate the bin given an alignment in [beg,end)
     * Copied from SAM spec.
     * @param beg 0-based start of read (inclusive)
     * @param end 0-based end of read (exclusive)
     */
    private static int reg2bin_original(final int beg, int end) {
        --end;

        if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
        if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
        if (beg>>20 == end>>20) return  ((1<<9)-1)/7 + (beg>>20);
        if (beg>>23 == end>>23) return  ((1<<6)-1)/7 + (beg>>23);
        if (beg>>26 == end>>26) return  ((1<<3)-1)/7 + (beg>>26);
        return 0;
    }

    /**
     * Get candidate bins for the specified region
     * @param startPos 1-based start of target region, inclusive.
     * @param endPos 1-based end of target region, inclusive.
     * @return bit set for each bin that may contain SAMRecords in the target region.
     */
    public static BitSet regionToBins(final int startPos, final int endPos) {
        final int adjustedEndPos;
        if (endPos > 0) adjustedEndPos = endPos;
        else adjustedEndPos = (int)computeMaxLen(BAI_MIN_SHIFT, BAI_DEPTH);
        final int adjustedStartPos = (startPos > 0? startPos - 1: 0);
        final BitSet ret = reg2bins(adjustedStartPos, adjustedEndPos, BAI_MIN_SHIFT, BAI_DEPTH);
        if (!ret.equals(regionToBins_original(startPos, endPos))) {
            throw new IllegalStateException(String.format("Results do not agree for startPos: %d endPos: %d", startPos, endPos));
        }
        return ret;
    }


    private static BitSet regionToBins_original(final int startPos, final int endPos) {
        final int maxPos = 0x1FFFFFFF;
        final int start = (startPos <= 0) ? 0 : (startPos-1) & maxPos;
        final int end = (endPos <= 0) ? maxPos : (endPos-1) & maxPos;
        if (start > end) {
            return null;
        }
        int k;
        final BitSet bitSet = new BitSet(GenomicIndexUtil.MAX_BINS);
        bitSet.set(0);
        for (k =    1 + (start>>26); k <=    1 + (end>>26); ++k) bitSet.set(k);
        for (k =    9 + (start>>23); k <=    9 + (end>>23); ++k) bitSet.set(k);
        for (k =   73 + (start>>20); k <=   73 + (end>>20); ++k) bitSet.set(k);
        for (k =  585 + (start>>17); k <=  585 + (end>>17); ++k) bitSet.set(k);
        for (k = 4681 + (start>>14); k <= 4681 + (end>>14); ++k) bitSet.set(k);
        return bitSet;
    }

    /**
     /**
     *
     * Copied from CSI spec.
     * @param beg 0-based start of read (inclusive)
     * @param end 0-based end of read (exclusive)
     * @param min_shift  number of bits for minimal interval
     * @param depth of binning index (based on max contig length)
     * @return calculate the bin given an alignment in [beg,end)
     */
    public static int reg2bin(final long beg, long end, int min_shift, int depth)
    {
        int l, s = min_shift, t = ((1<<depth*3) - 1) / 7;
        for (--end, l = depth; l > 0; --l, s += 3, t -= 1<<l*3)
            if (beg>>s == end>>s) return t + (int)(beg>>s);
        return 0;
    }

    /**
     * Get candidate bins for specified region
     * @param beg 0-based start, inclusive
     * @param end 0-based end, exclusive
     * @param min_shift number of bits for minimal interval
     * @param depth  of binning index (based on max contig length)
     * @return bit set for each bin that may contain records in the target region.
     */
    public static BitSet reg2bins(final long beg, long end, final int min_shift, int depth) {
        int l, t, n, s = min_shift + depth*3;
        // TODO: Figure out the right BitSet size based on depth (or maybe not)
        final BitSet bitSet = new BitSet(GenomicIndexUtil.MAX_BINS);
        for (--end, l = n = t = 0; l <= depth; s -= 3, t += 1<<l*3, ++l) {
            int b = t + (int)(beg>>s), e = t + (int)(end>>s);
            for (int i = b; i <= e; ++i) bitSet.set(i);
        }
        return bitSet;
    }

    /**
     * @param min_shift number of bits for minimal interval
     * @param max_len longest reference sequence to be indexed, or 0 to indicate unknown
     * @return depth to be passed to reg2bin and reg2bins
     */
    public static int computeDepth(final int min_shift, long max_len) {
        if ( max_len <= 0) max_len = ((long)1<<31) - 1;  // In case contig line is broken.
        max_len += 256;
        int depth;
        long s;
        for (depth = 0, s = 1<<min_shift; max_len > s; ++depth, s <<= 3);
        return depth;
    }

    /**
     * Inverse of computeDepth.  Note that a range of lengths will produce the same depth.  This returns the largest one.
     */
    public static long computeMaxLen(final int min_shift, final int depth) {
        long s = 1<<min_shift;
        for (int d = 0; d < depth; ++d, s <<= 3) {
            //
        }
        return s - 256;
    }

    /**
     * from hts.c
     * @param bin
     * @return parent bin number
     */
    private static int hts_bin_parent(final int bin) {
        return (bin - 1) >> 3;
    }

    /**
     * from hts.c
     * @param level
     * @return bin number of first bin at this level
     */
    private static int hts_bin_first(final int level) {
        return ((1<<((level<<1) + level)) - 1) / 7;
    }

    /**
     * from hts.c
     * @param bin
     * @param n_lvls
     * @return index into linear offset for bin
     */
    private static int hts_bin_bot(int bin, int n_lvls)
    {
        // compute the level of bin
        int level = 0;
        for (int b = bin; b > 0; ++level, b = hts_bin_parent(b));

        // Compute bin number relative to the first bin at this level,
        // then multiply by the number of slots in the linear index per bin at this level.
        return (bin - hts_bin_first(level)) << (n_lvls - level) * 3;
    }

    public static int getLinearIndexOffsetForBin(final int bin, final int depth) {
        return hts_bin_bot(bin, depth);
    }
}
