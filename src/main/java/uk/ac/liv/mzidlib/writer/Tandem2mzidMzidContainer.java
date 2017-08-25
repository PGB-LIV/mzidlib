/*
 * Date: 23-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.Tandem2mzidMzidContainer.java
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

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import de.proteinms.xtandemparser.xtandem.InputParams;
import de.proteinms.xtandemparser.xtandem.PerformParams;
import de.proteinms.xtandemparser.xtandem.XTandemFile;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.BibliographicReference;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.util.CVUtils;
import uk.ac.liv.mzidlib.util.Utils;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 23-Aug-2017 15:34:34
 */
public class Tandem2mzidMzidContainer implements MzidContainer {

    private final XTandemFile xfile;
    private MzIdentMLVersion version;
    private String databaseFileFormatId;
    private String databaseFileFormatName;
    private String massSpecFileFormatId;
    private String massSpecFileFormatName;
    private String decoyRegularExpression;
    private Pattern proteinCodeRegexPattern;
    private boolean outputFragmentation;
    private boolean isMs2SpectrumIdStartingAtZero;
    private Map<String, String> cvMap;
    
    private PerformParams tandemParams;
    private String tandemVersion;
    private String dbLocation;
    private String dbName;

    /**
     * Constructor.
     *
     * @param input                     input tandem file name
     * @param dbFileFormatId            database file format Id
     * @param msFileFormatId            mass spectrum file format Id
     * @param isMs2SpecIdStartingAtZero flag show if MS2 spectrum index starting
     *                                  at zero
     * @param decoyRegex                docoy database regular expression
     * @param proteinCodeRegex          protein code regular expression
     * @param outFragmentation          output fragmentation
     * @param ver                       Mzid file version
     *
     * @throws SAXException                 SAX exception
     * @throws ParserConfigurationException parser configuration exception
     * @throws IOException                  IO exception
     */
    public Tandem2mzidMzidContainer(String input, String dbFileFormatId,
                                    String msFileFormatId,
                                    Boolean isMs2SpecIdStartingAtZero,
                                    String decoyRegex,
                                    String proteinCodeRegex,
                                    boolean outFragmentation,
                                    MzIdentMLVersion ver)
            throws SAXException, ParserConfigurationException, IOException {
        this.xfile = new XTandemFile(input);

        init(xfile, dbFileFormatId, msFileFormatId, isMs2SpecIdStartingAtZero,
             decoyRegex, proteinCodeRegex, outFragmentation, ver);

    }

    @Override
    public AnalysisCollection getAnalysisCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisData getAnalysisData() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisSampleCollection getAnalysisSampleCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AuditCollection getAuditCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public BibliographicReference getBibliographicReference() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public CvList getCvList() {
        CvList cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        localCvList.add(CvConstants.PSI_CV);
        localCvList.add(CvConstants.UNIMOD_CV);
        localCvList.add(CvConstants.UNIT_CV);
        return cvList;
    }

    @Override
    public DataCollection getDataCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Inputs getInputs() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public ProteinDetectionList getProteinDetectionList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Provider getProvider() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public SequenceCollection getSequenceCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public SpectrumIdentificationList getSpectrumIdentificationList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private String getCvName(String cvItemId)
            throws IOException {
        //If CV map is not yet initialized, do it:
        if (this.cvMap == null) {
            this.cvMap = Utils.getInitializedCVMap();
        }

        //validate:
        if (this.cvMap.get(cvItemId) == null) {
            throw new RuntimeException(
                    "Given item not found in Controlled Vocabulary : "
                    + cvItemId);
        } else {
            return this.cvMap.get(cvItemId);
        }
    }

    private void init(XTandemFile xfile, String dbFileFormatId,
                      String msFileFormatId, Boolean isMs2SpecIdStartingAtZero,
                      String decoyRegex, String proteinCodeRegex,
                      boolean outFragmentation,
                      MzIdentMLVersion ver)
            throws IOException {
        // decoyRegularExpression
        if (decoyRegex != null && !decoyRegex.trim().isEmpty()) {
            this.decoyRegularExpression = decoyRegex;
        } else {
            this.decoyRegularExpression = null;
        }

        // proteinCodeRegexPattern
        if (proteinCodeRegex != null && !proteinCodeRegex.trim().isEmpty()) {
            //this regex should ensure the protein code is parsed from the longer 
            //string that X!Tandem is currently making of this. 
            //See also https://code.google.com/p/mzidentml-lib/issues/detail?id=14
            try {
                this.proteinCodeRegexPattern = Pattern.compile(proteinCodeRegex);
            } catch (PatternSyntaxException pe) {
                throw new RuntimeException("Given proteinCodeRegex ["
                        + proteinCodeRegex + "] is not a "
                        + "valid regular expression.", pe);
            }
        } else {
            this.proteinCodeRegexPattern = null;
        }

        // outputFragementation
        this.outputFragmentation = outFragmentation;

        // databaseFileFormatId
        if (dbFileFormatId != null) {
            this.databaseFileFormatId = dbFileFormatId;
            this.databaseFileFormatName = getCvName(databaseFileFormatId);
        } else {
            this.tandemParams = xfile.getPerformParameters();
            this.dbLocation = tandemParams.getSequenceSource_1();
            String[] cvIdAndName = CVUtils.getDatabaseFileFormat(dbLocation);
            this.databaseFileFormatId = cvIdAndName[0];
            this.databaseFileFormatName = cvIdAndName[1];
        }

        // massSpecFileFormatId
        if (msFileFormatId != null) {
            this.massSpecFileFormatId = msFileFormatId;
            this.massSpecFileFormatName = getCvName(massSpecFileFormatId);
        } else {
            //Try to infer from the file itself:
            InputParams inputParams = xfile.getInputParameters();
            //Validate: if the spectrum path is null, then we can assume all 
            //input parameters are missing as the spectrum path 
            //is a mandatory for X!tandem to run:
            String spectrumFile = inputParams.getSpectrumPath();
            if (spectrumFile == null) {
                throw new RuntimeException(
                        "Expected parameter 'spectrum, path' not found in "
                        + "X!Tandem file. Please run your X!Tandem search with "
                        + "option 'output, parameters=yes'. See "
                        + "http://thegpm.org/tandem/api/opara.html for more details.");
            }
            String[] cvIdAndName = CVUtils.getMassSpecFileFormatID(spectrumFile);
            this.massSpecFileFormatId = cvIdAndName[0];
            this.massSpecFileFormatName = cvIdAndName[1];
        }

        // isMs2SpectrumIdStartingAtZero
        if (isMs2SpecIdStartingAtZero == null) {
            //if file format is mzML (MS:1000584), then spectrum starts at 0, 
            //otherwise it starts at 1
            if (this.massSpecFileFormatId.equalsIgnoreCase("MS:1000584")) {
                this.isMs2SpectrumIdStartingAtZero = true;
            } else {
                this.isMs2SpectrumIdStartingAtZero = false;
            }
        } else {
            this.isMs2SpectrumIdStartingAtZero = isMs2SpecIdStartingAtZero;
        }

        // MzIdentMLVersion
        if (null == ver) {
            this.version = MzIdentMLVersion.Version_1_1;
        } else {
            this.version = ver;
        }
    }

}
