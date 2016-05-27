package uk.ac.liv.mzidlib.util;

import java.io.File;
import java.util.Iterator;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author jonesar
 */
public class MzidLibUtils {

    public CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
    }

    public CvParam makeCvParam(String accession, String name, Cv cv, String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setValue(value);
        return cvParam;
    }

    public UserParam makeUserParam(String name, String value) {
        UserParam userParam = new UserParam();
        userParam.setName(name);
        userParam.setValue(value);
        return userParam;
    }

    public Cv getpsiCV(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("psi")) {
                cv1 = cv;
                break;
            }
        }
        return cv1;
    }

    public Cv getUnimod(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("unimod")) {
                cv1 = cv;
                break;
            }
        }
        return cv1;
    }

    public Cv getUnit(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("unit")) {
                cv1 = cv;
                break;
            }
        }
        return cv1;
    }

    public static int safeLongToInt(long l) {
        if (l < Integer.MIN_VALUE || l > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(l + " cannot be cast to int without changing its value.");
        }
        return (int) l;
    }

}
