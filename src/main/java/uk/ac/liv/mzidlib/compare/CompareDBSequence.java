package uk.ac.liv.mzidlib.compare;

import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;

/**
 *
 * @author Fawaz Ghali
 */
public class CompareDBSequence {

    public static boolean compareDBSequence(DBSequence db1, DBSequence db2) {
        return db1.getAccession().equals(db2.getAccession());
    }
}
