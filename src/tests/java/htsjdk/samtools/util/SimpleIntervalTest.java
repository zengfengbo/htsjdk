package htsjdk.samtools.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SimpleIntervalTest {

    @Test
    public void testZeroLengthInterval(){
        SimpleInterval interval = new SimpleInterval("1",100,99);
        Assert.assertEquals(interval.getContig(), "1");
        Assert.assertEquals(interval.getStart(), 100);
        Assert.assertEquals(interval.getEnd(), 99);
    }

    @DataProvider(name = "badIntervals")
    public Object[][] badIntervals(){
        return new Object[][]{
                {null,1,12, "null contig"},
                {"1", 0, 10, "start==0"},
                {"1", -10, 10, "negative start"},
                {"1", 10, 8, "end < start - 1"}
        };
    }

    @Test(dataProvider = "badIntervals", expectedExceptions = IllegalArgumentException.class)
    public void badIntervals(String contig, int start, int end, String name){
        SimpleInterval interval = new SimpleInterval(contig, start, end);
    }
}