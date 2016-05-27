package uk.ac.liv.mzidlib.multiplesearch;

public class myDataMissingException extends Exception {

    private static final long serialVersionUID = 1L;

    public myDataMissingException(String msg) {
        super(msg);
    }

    public myDataMissingException(String msg, Throwable t) {
        super(msg, t);
    }

}
