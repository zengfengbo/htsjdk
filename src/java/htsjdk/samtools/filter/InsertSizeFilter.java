package htsjdk.samtools.filter;

import htsjdk.samtools.SAMRecord;

/**
 * SAM filter for insert size range.
 */
public class InsertSizeFilter implements SamRecordFilter {
    final int minInsertSize;
    final int maxInsertSize;

    public InsertSizeFilter(final int minInsertSize, final int maxInsertSize) {
        this.minInsertSize = minInsertSize;
        this.maxInsertSize = maxInsertSize;
    }

    @Override
    public boolean filterOut(final SAMRecord rec) {
        // Treat both parameters == 0 as not filtering
        if (minInsertSize == 0 && maxInsertSize == 0) return false;

        if (rec.getReadPairedFlag()) {
            final int ins = Math.abs(rec.getInferredInsertSize());
            return ins < minInsertSize || ins > maxInsertSize;
        }

        // If the read isn't paired and either min or max is specified filter it out
        return minInsertSize != 0 || maxInsertSize != 0;
    }

    @Override
    public boolean filterOut(final SAMRecord r1, final SAMRecord r2) {
        return filterOut(r1) || filterOut(r2);
    }
}