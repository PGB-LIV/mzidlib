
package uk.ac.liv.mzidlib.util;

import java.util.Iterator;
import java.util.List;

import de.proteinms.xtandemparser.xtandem.Domain;
import de.proteinms.xtandemparser.xtandem.FragmentIon;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.SubstitutionModification;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.ReadUnimod;
import uk.ac.liv.unimod.ModT;

/**
 *
 * @author jonesar
 */
public class MzidLibUtils {

    public static CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
    }

    public static CvParam makeCvParam(String accession, String name, Cv cv,
                                      String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setValue(value);
        return cvParam;
    }

    public static Cv makeCv(String id, String uri, String name, String version) {
        Cv retCv = new Cv();
        retCv.setId(id);
        retCv.setFullName(name);
        retCv.setUri(uri);
        retCv.setVersion(version);
        return retCv;
    }

    public static Cv makeCv(String id, String uri, String name) {
        Cv retCv = new Cv();
        retCv.setId(id);
        retCv.setFullName(name);
        retCv.setUri(uri);
        return retCv;
    }

    public static CvParam makeCvParam(String accession, String name, Cv cv,
                                      String unitAccession, String unitName,
                                      Cv alternateUnitCv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setUnitAccession(unitAccession);
        cvParam.setUnitCv(alternateUnitCv);
        cvParam.setUnitName(unitName);
        return cvParam;
    }

    public static UserParam makeUserParam(String name, String value) {
        UserParam userParam = new UserParam();
        userParam.setName(name);
        userParam.setValue(value);
        return userParam;
    }

    public static Cv getpsiCV(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(
                MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("psi")) {
                cv1 = cv;
                break;
            }
        }
        return cv1;
    }

    public static Cv getUnimod(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(
                MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("unimod")) {
                cv1 = cv;
                break;
            }
        }
        return cv1;
    }

    public static Cv getUnit(MzIdentMLUnmarshaller unmarshaller) {
        Cv cv1 = null;

        Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(
                MzIdentMLElement.CV);
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
            throw new IllegalArgumentException(l
                    + " cannot be cast to int without changing its value.");
        }
        return (int) l;
    }

    public static CvParam getCvParamWithMassUnits(boolean isDaltonUnit) {
        CvParam cvParam = new CvParam();

        //<cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
        cvParam.setCv(CvConstants.PSI_CV);
        cvParam.setUnitCv(CvConstants.UNIT_CV);

        if (isDaltonUnit) {
            cvParam.setUnitAccession("UO:0000221");
            cvParam.setUnitName("dalton");
        } else {
            cvParam.setUnitAccession("UO:0000169");
            cvParam.setUnitName("parts per million");
        }
        return cvParam;
    }

    public static Modification translateToMzidModification(
            de.proteinms.xtandemparser.interfaces.Modification reportedMod,
            Domain domain, boolean fragmentIsMono, double unimodMassError) {
        ReadUnimod unimodDoc = new ReadUnimod();
        Modification mzidMod = new Modification();
        double mass = reportedMod.getMass();
        int loc = Integer.parseInt(reportedMod.getLocation());

        CvParam modParam = new CvParam();

        if (fragmentIsMono) {
            mzidMod.setMonoisotopicMassDelta(mass);
        } else {
            mzidMod.setAvgMassDelta(mass);
        }

        int pepLoc = loc - domain.getDomainStart(); //location in Tandem is given as location within the whole protein
        char reportedAminoAcid = domain.getDomainSequence().charAt(pepLoc);
        mzidMod.setLocation(pepLoc + 1);        //mzid starts counting from 1, except for NTerm/CTerm mods which are 0 
        List<String> residueList = mzidMod.getResidues();
        residueList.add("" + reportedAminoAcid);

        //check if we can find a known modification type that matches what is reported in the given mass and aminoAcid, within
        //a given mass error tolerance: 
        ModT unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                             fragmentIsMono, reportedAminoAcid);

        //Part below is needed because it is not clear in XTandem output what are N or C-term mods, 
        //these have the same location as mods on the first aa or last aa in the peptide.
        //So if no unimod was found above, perhaps this is because it is really a N or C-term modification. 
        //Below we try to see if it fits a known N or C-term modification:
        if (unimod == null && pepLoc == 0) {
            //See if this is a possible N-terminal mod

            //System.out.println("\tModification not found in Unimod, so look to see if it is N-terminal\n");            
            unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                            fragmentIsMono, '[');
            mzidMod.setLocation(0);
        }
        if (unimod == null && loc == domain.getDomainEnd()) {
            //See if this is a possible C-terminal mod
            //System.out.println("\tNot found, so look to see if it is N-terminal\n");
            unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                            fragmentIsMono, ']');
            mzidMod.setLocation(0); //also 0?
        }

        //Set the found details to modParam. If no unimod record was found, set the modification as "unknown" 
        if (unimod != null) {
            modParam.setAccession("UNIMOD:" + unimod.getRecordId());
            modParam.setCv(CvConstants.UNIMOD_CV);
            modParam.setName(unimod.getTitle());
        } else {
            //modification with mass not recognized:
            modParam.setName("unknown modification");
            modParam.setCv(CvConstants.PSI_CV);
            modParam.setAccession("MS:1001460");
        }

        mzidMod.getCvParam().add(modParam);
        return mzidMod;
    }

    public static SubstitutionModification translateToMzidSubstitution(
            de.proteinms.xtandemparser.interfaces.Modification reportedMod,
            Domain domain, boolean fragmentIsMono) {
        SubstitutionModification mzidSubs = new SubstitutionModification();

        int loc = Integer.parseInt(reportedMod.getLocation());

        int pepLoc = loc - domain.getDomainStart(); //location in Tandem is given as location within the whole protein
        mzidSubs.setLocation(pepLoc + 1);

        //Note: in X!Tandem, like in mzIdentML: sequence reported is the original:
        char reportedAminoAcid = domain.getDomainSequence().charAt(pepLoc);
        mzidSubs.setOriginalResidue(reportedAminoAcid + "");

        mzidSubs.setReplacementResidue(reportedMod.getSubstitutedAminoAcid());

        //TODO - the current storing of the modification mass in either
        //setMonoisotopicMassDelta or setAvgMassDelta seems wrong as the massDelta is 
        //a delta on the parent ion mass which is always monoisotopic in xtandem. On the other hand, 
        //a modification has to have one or more matching fragments as well which all are 
        //shifted by x, x being the modification mass.
        //So check with X!Tandem developer whether the modification (and substitution) masses 
        //follow the fragment mass type or the parent mass type.
        double mass = reportedMod.getMass();

        if (fragmentIsMono) {
            mzidSubs.setMonoisotopicMassDelta(mass);
        } else {
            mzidSubs.setAvgMassDelta(mass);
        }

        return mzidSubs;

    }

    public static CvParam getFragmentCVParam(int iontype) {

        CvParam cvParam = new CvParam();
        cvParam.setCv(CvConstants.PSI_CV);

        //nb: don't forget the "break" when adding new ones here!
        switch (iontype) {
            case FragmentIon.A_ION:
                cvParam.setAccession("MS:1001229");
                cvParam.setName("frag: a ion");
                break;
            case FragmentIon.AH2O_ION:
                cvParam.setAccession("MS:1001234");
                cvParam.setName("frag: a ion - H2O");
                break;
            case FragmentIon.ANH3_ION:
                cvParam.setAccession("MS:1001235");
                cvParam.setName("frag: a ion - NH3");
                break;
            case FragmentIon.B_ION:
                cvParam.setAccession("MS:1001224");
                cvParam.setName("frag: b ion");
                break;
            case FragmentIon.BH2O_ION:
                cvParam.setAccession("MS:1001222");
                cvParam.setName("frag: b ion - H2O");
                break;
            case FragmentIon.BNH3_ION:
                cvParam.setAccession("MS:1001232");
                cvParam.setName("frag: b ion - NH3");
                break;
            case FragmentIon.C_ION:
                cvParam.setAccession("MS:1001231");
                cvParam.setName("frag: c ion");
                break;
            case FragmentIon.X_ION:
                cvParam.setAccession("MS:1001228");
                cvParam.setName("frag: x ion");
                break;
            case FragmentIon.Y_ION:
                cvParam.setAccession("MS:1001220");
                cvParam.setName("frag: y ion");
                break;
            case FragmentIon.YH2O_ION:
                cvParam.setAccession("MS:1001223");
                cvParam.setName("frag: y ion - H2O");
                break;
            case FragmentIon.YNH3_ION:
                cvParam.setAccession("MS:1001233");
                cvParam.setName("frag: y ion - NH3");
                break;
            case FragmentIon.Z_ION:
                cvParam.setAccession("MS:1001230");
                cvParam.setName("frag: z ion");
                break;
            case FragmentIon.MH_ION:
                cvParam.setAccession("MS:1001523");
                cvParam.setName("frag: precursor ion");
                break;
            case FragmentIon.MHNH3_ION:
                cvParam.setAccession("MS:1001522");
                cvParam.setName("frag: precursor ion - NH3");
                break;
            case FragmentIon.MHH2O_ION:
                cvParam.setAccession("MS:1001521");
                cvParam.setName("frag: precursor ion - H2O");
                break;
        }

        return cvParam;

    }

}
