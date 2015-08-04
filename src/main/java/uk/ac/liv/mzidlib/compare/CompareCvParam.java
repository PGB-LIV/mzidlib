package uk.ac.liv.mzidlib.compare;

import java.util.List;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

/**
 *
 * @author Fawaz Ghali
 */
public class CompareCvParam extends CvParam {

	private static final long serialVersionUID = 1L;

	public static boolean compareCvParam(CvParam p1, CvParam p2) {
     
        if ((p1.getAccession() != null && p2.getAccession() == null) || (p1.getAccession() == null && p2.getAccession() != null)) {
            return false;
        }
        if (p1.getAccession() != null && p2.getAccession() != null && !p1.getAccession().equals(p2.getAccession())) {
            return false;
        }

        if ((p1.getValue() != null && p2.getValue() == null) || (p1.getValue() == null && p2.getValue() != null)) {
            return false;
        }
        if (p1.getValue() != null && p2.getValue() != null && !p1.getValue().equals(p2.getValue())) {
            return false;
        }
        return true;
    }

    public static boolean compareCvParam(List<CvParam> p1, List<CvParam> p2) {
        if (p1.size() != p2.size()) {
            return false;
        }
        for (int i = 0; i < p1.size(); i++) {
            CvParam cvParam1 = p1.get(i);
            boolean found = false;
            for (int j = 0; j < p2.size(); j++) {
                CvParam cvParam2 = p2.get(j);
                if (compareCvParam(cvParam1, cvParam2)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }
        return true;


    }
}
