package htsjdk.samtools;

/**
 * Abstract ReadRecord class that contains implementations of methods common to all ReadRecord subtyptes.
 */
public abstract class AbstractReadRecord implements ReadRecord {

    @Override
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

}
