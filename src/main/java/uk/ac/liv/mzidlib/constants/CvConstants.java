/*
 * Date: 21-Mar-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.constants.CvConstants.java
 *
 */

package uk.ac.liv.mzidlib.constants;

import uk.ac.ebi.jmzidml.model.mzidml.Cv;

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

    public static Cv makeCv(String id, String uri, String name, String version) {
        Cv retCv = new Cv();
        retCv.setId(id);
        retCv.setFullName(name);
        retCv.setUri(uri);
        retCv.setVersion(version);
        return retCv;
    }

}
