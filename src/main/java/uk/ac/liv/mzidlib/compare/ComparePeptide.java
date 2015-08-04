package uk.ac.liv.mzidlib.compare;

import java.util.List;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;

/**
 *
 * @author Fawaz Ghali
 */
public class ComparePeptide {

    public static boolean comparePeptide(Peptide p1, Peptide p2) {
        if (!p1.getPeptideSequence().equals(p2.getPeptideSequence())) {
            return false;
        } else {

            List<Modification> m1List = p1.getModification();
            List<Modification> m2List = p2.getModification();
            if (m1List == null && m2List != null) {
                return false;
            }
            if (m1List != null && m2List == null) {
                return false;
            }
            if (m1List != null && m2List != null) {
                boolean comapreList = CompareModification.compareModification(m1List, m2List);
                if (!comapreList) {
                    return false;
                }
            }
            List<CvParam> c1List = p1.getCvParam();
            List<CvParam> c2List = p2.getCvParam();
            if (c1List == null && c2List != null) {
                return false;
            }
            if (c1List != null && c2List == null) {
                return false;
            }
            if (c1List != null && c2List != null) {
                boolean comapreList = CompareCvParam.compareCvParam(c1List, c2List);
                if (!comapreList) {
                    return false;
                }
            }

        }

        return true;
    }
}
