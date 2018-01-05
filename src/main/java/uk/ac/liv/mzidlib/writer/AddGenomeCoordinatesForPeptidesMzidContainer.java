/*
 * Date: 27-Sep-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.AddGenomeCoordinatesForPeptidesMzidContainer.java
 *
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

import java.io.File;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.BibliographicReference;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 27-Sep-2017 10:52:14
 */
public class AddGenomeCoordinatesForPeptidesMzidContainer implements
        MzidContainer {

    private String inputMzidFn;
    private MzIdentMLUnmarshaller um;
    private MzIdentMLVersion version;

    public AddGenomeCoordinatesForPeptidesMzidContainer(String input) {
        this.inputMzidFn = input;
        init();

    }

    @Override
    public AnalysisCollection getAnalysisCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        return um.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
    }

    @Override
    public AnalysisSampleCollection getAnalysisSampleCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        return um.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
    }

    @Override
    public AuditCollection getAuditCollection() {
        return um.unmarshal(MzIdentMLElement.AuditCollection);
    }

    @Override
    public BibliographicReference getBibliographicReference() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public CvList getCvList() {
        return this.um.unmarshal(MzIdentMLElement.CvList);
    }

    @Override
    public Inputs getInputs() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public MzIdentMLVersion getMzidVersion() {
        return this.version;
    }

    @Override
    public ProteinDetectionList getProteinDetectionList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Provider getProvider() {
        return um.unmarshal(MzIdentMLElement.Provider);
    }

    @Override
    public SequenceCollection getSequenceCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public SpectrumIdentificationList getSpectrumIdentificationList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private void init() {
        this.um = new MzIdentMLUnmarshaller(new File(this.inputMzidFn));
    }

}
