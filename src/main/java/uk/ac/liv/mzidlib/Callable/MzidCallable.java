package uk.ac.liv.mzidlib.Callable;

import java.util.concurrent.Callable;
import uk.ac.liv.mzidlib.MzIdentMLLib;

/**
 *
 * @author Fawaz Ghali 30-Apr-2015
 */
public class MzidCallable implements Callable< Boolean> {

    String[] args;

    public MzidCallable(String[] args) {
        this.args = args;

    }

    public Boolean call() {
        try {
            MzIdentMLLib mzidLib = new MzIdentMLLib();
            mzidLib.init(args);
            return true;
        } catch (Exception e) {
            e.printStackTrace();
            return false;
        }
    }

}
