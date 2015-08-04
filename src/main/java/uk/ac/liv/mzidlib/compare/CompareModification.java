package uk.ac.liv.mzidlib.compare;

import java.util.List;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;

/**
 *
 * @author Fawaz Ghali
 */
public class CompareModification {

    public static boolean compareModification(Modification m1, Modification m2) {
        if ((m1.getMonoisotopicMassDelta()== null && m2.getMonoisotopicMassDelta() != null)
                || (m1.getMonoisotopicMassDelta() != null && m2.getMonoisotopicMassDelta() == null)
                || (m1.getMonoisotopicMassDelta() != null && m2.getMonoisotopicMassDelta() != null && m1.getMonoisotopicMassDelta().doubleValue() != m2.getMonoisotopicMassDelta().doubleValue())) {
            return false;
        }
        if ((m1.getLocation() == null && m2.getLocation() != null)
                || (m1.getLocation() != null && m2.getLocation() == null)
                || (m1.getLocation() != null && m2.getLocation() != null && m1.getLocation().intValue() != m2.getLocation().intValue())) {
            return false;
        }

        if ((m1.getCvParam() == null && m2.getCvParam() != null)
                || (m1.getCvParam() != null && m2.getCvParam() == null)) {
            return false;
        }

        if (m1.getCvParam() != null && m2.getCvParam() != null) {
            boolean compareCvParamList = CompareCvParam.compareCvParam(m1.getCvParam(), m2.getCvParam());
            if (!compareCvParamList) {
                return false;
            }
        }

        if ((m1.getResidues() == null && m2.getResidues() != null)
                || (m1.getResidues() != null && m2.getResidues() == null)) {
            return false;
        }

        if (m1.getResidues() != null && m2.getResidues() != null) {
            boolean equalLists = m1.getResidues().size() == m2.getResidues().size() && m1.getResidues().containsAll(m2.getResidues());
            if (!equalLists) {
                return false;
            }
        }

        return true;
    }

    public static boolean compareModification(List<Modification> p1, List<Modification> p2) {
        if (p1.size() != p2.size()) {
            return false;
        }
        for (int i = 0; i < p1.size(); i++) {
            Modification modification1 = p1.get(i);
            boolean found = false;
            for (int j = 0; j < p2.size(); j++) {
                Modification modification2 = p2.get(j);
                if (compareModification(modification1, modification2)) {
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
