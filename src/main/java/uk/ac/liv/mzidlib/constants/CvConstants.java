/*
 * Date: 21-Mar-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.constants.CvConstants.java
 *
 */

package uk.ac.liv.mzidlib.constants;

import static uk.ac.liv.mzidlib.util.MzidLibUtils.makeCvParam;

import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @since 21-Mar-2017 09:56:07
 */
public class CvConstants {

    private static final String PSI_MS_VERSION = "4.0.8";
    private static final String PSI_MS_URI
            = "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo";
    
    public static final Cv PSI_CV = makeCv("PSI-MS", PSI_MS_URI, "PSI-MS", PSI_MS_VERSION);
    
    public static final CvParam NO_SPECIAL_PROCESSING = makeCvParam("MS:1002495", "no special processing", PSI_CV);
    public static final CvParam PEPTIDE_LEVEL_SCORING = makeCvParam("MS:1002490", "peptide-level scoring", PSI_CV);
    public static final CvParam MODIFICATION_LOCOLIZATION_SCORING = makeCvParam("MS:1002491", "modification localization scoring", PSI_CV);
    public static final CvParam CONSENSUS_SCORING = makeCvParam("MS:1002492", "consensus scoring", PSI_CV);
    public static final CvParam SAMPLE_PRE_FRACTIONATION = makeCvParam("MS:1002493", "sample pre-fractionation", PSI_CV);
    public static final CvParam CROSS_LINKING_SEARCH = makeCvParam("MS:1002494", "cross-linking search", PSI_CV);
    public static final CvParam DE_NOVO_SEARCH = makeCvParam("MS:1001010", "de novo search", PSI_CV);
    public static final CvParam SPECTRAL_LIBRARY_SEARCH = makeCvParam("MS:1001031", "spectral library search", PSI_CV);
    public static final CvParam PROTEOGENOMICS_SEARCH = makeCvParam("MS:1002635", "proteogenomics search", PSI_CV);
        
    public static Cv makeCv(String id, String uri, String name, String version) {
        Cv retCv = new Cv();
        retCv.setId(id);
        retCv.setFullName(name);
        retCv.setUri(uri);
        retCv.setVersion(version);
        return retCv;
    }

}
