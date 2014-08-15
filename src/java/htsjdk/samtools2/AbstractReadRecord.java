package htsjdk.samtools2;

import htsjdk.samtools2.ReadRecord;

abstract class AbstractReadRecord implements ReadRecord {

    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

}