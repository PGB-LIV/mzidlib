/*
 * Date: 21-Mar-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.constants.CvConstants.java
 *
 * jmzquantml is Copyright 2017 University of Liverpool.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
