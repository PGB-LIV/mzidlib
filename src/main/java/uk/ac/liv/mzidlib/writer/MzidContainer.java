/*
 * Date: 23-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.MzidContainer.java
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

package uk.ac.liv.mzidlib.writer;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.BibliographicReference;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;

/**
 * The container interface to return MzIdentML elements.
 *
 * @author Da Qi
 * University of Liverpool
 * @time 23-Aug-2017 14:26:58
 */
public interface MzidContainer {

    public AnalysisCollection getAnalysisCollection();

    //public AnalysisData getAnalysisData();

    public AnalysisProtocolCollection getAnalysisProtocolCollection();

    public AnalysisSampleCollection getAnalysisSampleCollection();

    public AnalysisSoftwareList getAnalysisSoftwareList();

    public AuditCollection getAuditCollection();

    public BibliographicReference getBibliographicReference();

    public CvList getCvList();

    //public DataCollection getDataCollection();

    public Inputs getInputs();

    public ProteinDetectionList getProteinDetectionList();

    public Provider getProvider();

    public SequenceCollection getSequenceCollection();

    public SpectrumIdentificationList getSpectrumIdentificationList();

}
