/*
 * Date: 21-Mar-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.constants.CvConstants.java
 *
 */

package uk.ac.liv.mzidlib.constants;

import static uk.ac.liv.mzidlib.util.MzidLibUtils.makeCv;
import static uk.ac.liv.mzidlib.util.MzidLibUtils.makeCvParam;

import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

/**
 *
 * @author Da Qi
 * University of Liverpool
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
    public static final CvParam GROUP_PSMS_BY_SEQUENCE = makeCvParam("MS:1002496", "group PSMs by sequence", PSI_CV);
    public static final CvParam MZIDLIB = makeCvParam("MS:1002237", "mzidLib", PSI_CV);
    
    //FDR calculation
    public static final CvParam XTANDEM_EXPECT = makeCvParam("MS:1001330", "X\\!Tandem:expect", PSI_CV);
    public static final CvParam MASCOT_EXPECTATION_VALUE = makeCvParam("MS:1001172", "Mascot:expectation value", PSI_CV);
    public static final CvParam SEQUEST_EXPECTATION_VALUE = makeCvParam("MS:1001159", "SEQUEST:expectation value", PSI_CV);
    public static final CvParam OMSSA_EVALUE = makeCvParam("MS:1001328", "OMSSA:evalue", PSI_CV);
    public static final CvParam PROTEINPROSPECTOR_EXPECTATION_VALUE = makeCvParam("MS:1002045", "ProteinProspector:expectation value", PSI_CV);
    public static final CvParam MSGF_EVALUE = makeCvParam("MS:1002053", "MS-GF:EValue", PSI_CV);
    //Peptide level FDR
    public static final CvParam DISTINCT_PEPTIDE_LEVEL_LOCAL_FDR = makeCvParam("MS:1002359", "distinct peptide-level local FDR", PSI_CV);
    public static final CvParam DISTINCT_PEPTIDE_LEVEL_Q_VALUE = makeCvParam("MS:1001868", "distinct peptide-level q-value", PSI_CV);
    public static final CvParam DISTINCT_PEPTIDE_LEVEL_FDRSCORE = makeCvParam("MS:1002360", "distinct peptide-level FDRScore", PSI_CV);
    //PSM level FDR
    public static final CvParam PSM_LEVEL_LOCAL_FDR = makeCvParam("MS:1002351", "PSM-level local FDR", PSI_CV);
    public static final CvParam PSM_LEVEL_Q_VALUE = makeCvParam("MS:1002354", "PSM-level q-value", PSI_CV);
    public static final CvParam PSM_LEVEL_FDRSCORE = makeCvParam("MS:1002355", "PSM-level FDRScore", PSI_CV);
    
}
