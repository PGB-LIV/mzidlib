/*
 * Date: 23-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.Tandem2mzidMzidContainer.java
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import javax.xml.parsers.ParserConfigurationException;

import de.proteinms.xtandemparser.xtandem.Domain;
import de.proteinms.xtandemparser.xtandem.FragmentIon;
import de.proteinms.xtandemparser.xtandem.InputParams;
import de.proteinms.xtandemparser.xtandem.PerformParams;
import de.proteinms.xtandemparser.xtandem.Protein;
import de.proteinms.xtandemparser.xtandem.Spectrum;
import de.proteinms.xtandemparser.xtandem.SupportData;
import de.proteinms.xtandemparser.xtandem.XTandemFile;
import org.xml.sax.SAXException;
import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.BibliographicReference;
import uk.ac.ebi.jmzidml.model.mzidml.ContactRole;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Enzymes;
import uk.ac.ebi.jmzidml.model.mzidml.FileFormat;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentArray;
import uk.ac.ebi.jmzidml.model.mzidml.Fragmentation;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.InputSpectra;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.Measure;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.ModificationParams;
import uk.ac.ebi.jmzidml.model.mzidml.Organization;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.Person;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.Role;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SourceFile;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIDFormat;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.SubstitutionModification;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.ReadUnimod;
import uk.ac.liv.mzidlib.util.CVUtils;
import uk.ac.liv.mzidlib.util.MzidLibUtils;
import uk.ac.liv.mzidlib.util.Utils;
import uk.ac.liv.unimod.ModT;

/**
 * Tandem2mzidMzidContainer.
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 23-Aug-2017 15:34:34
 */
public class Tandem2mzidMzidContainer implements MzidContainer {

    private final XTandemFile xfile;
    private final String inputFileName;
    private static final String XTANDEM_REVERSED_FLAG = ":reversed";
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
    private boolean fragmentIsMono;
    private Measure mzMeasure;
    private Measure intMeasure;
    private Measure errorMeasure;
    private int pepEvidCounter = 0;

    private PerformParams tandemParams;
    private String tandemVersion;
    private String dbLocation;
    private String dbName;

    private Person docOwner;
    private AnalysisSoftware analysisSoftwareXtandem;
    private SpectrumIdentificationProtocol siProtocol;
    private SearchDatabase searchDb;
    private SequenceCollection sequenceCollection;
    private SpectrumIdentificationList siList;

    private static final String ANALYSIS_SOFT_ID = "ID_software";
    private static final String SI_PROTOCOL_ID = "SearchProtocol_1";
    private static final String SOURCE_FILE_ID = "SourceFile_1";
    private static final String SEARCH_DB_ID = "SearchDB_1";
    private static final String SPECTRA_DATA_ID = "SID_1";
    private static final String SPECT_IDENT_ID = "SpectIdent_1";
    private static final String SI_LIST_ID = "SI_List_1";
    private static final String MEASURE_MZ_ID = "Measure_MZ";
    private static final String MEASURE_INT_ID = "Measure_Int";
    private static final String MEASURE_ERROR_ID = "Measure_Error";
    //TODO This parameter is hard-coded (ARJ changed from 0.001 to 0.01 - Aug2012; 
    //perhaps should be set dynamically from search params)
    private double unimodMassError = 0.01;
    private final InputParams inputParams;
    private final ReadUnimod unimodDoc;

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
        this.inputFileName = input;
        this.xfile = new XTandemFile(input);

        init(xfile, dbFileFormatId, msFileFormatId, isMs2SpecIdStartingAtZero,
             decoyRegex, proteinCodeRegex, outFragmentation, ver);

        inputParams = xfile.getInputParameters();
        //If the spectrum path is null, then we can assume all input parameters are missing
        //as the spectrum path is a mandatory for X!tandem to run:
        if (inputParams.getSpectrumPath() == null) {
            throw new RuntimeException(
                    "Expected parameter not found in X!Tandem file. "
                    + "Please run your X!Tandem search with option 'output, "
                    + "parameters=yes'. "
                    + "See http://thegpm.org/tandem/api/opara.html "
                    + "for more details.");
        }

        fragmentIsMono = isMonoFragment(inputParams);

        if (outputFragmentation) {
            mzMeasure = new Measure();
            mzMeasure.setId(MEASURE_MZ_ID);
            List<CvParam> cvParamList = mzMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001225",
                                                     "product ion m/z",
                                                     CvConstants.PSI_CV,
                                                     "MS:1000040", "m/z",
                                                     CvConstants.PSI_CV));
            intMeasure = new Measure();
            intMeasure.setId(MEASURE_INT_ID);
            cvParamList = intMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001226",
                                                     "product ion intensity",
                                                     CvConstants.PSI_CV,
                                                     "MS:1000131",
                                                     "number of counts",
                                                     CvConstants.PSI_CV));
            errorMeasure = new Measure();
            errorMeasure.setId(MEASURE_ERROR_ID);
            cvParamList = errorMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001227",
                                                     "product ion m/z error",
                                                     CvConstants.PSI_CV,
                                                     "MS:1000040", "m/z",
                                                     CvConstants.PSI_CV));
        }

        unimodDoc = new ReadUnimod();

        //create SpectrumIdentificationProtocol
        siProtocol = createSpectrumIdentificationProtocol();

        //create SearchDatabase
        searchDb = createSearchDatabase();

        createSpectrumIdentificationListAndSequenceCollection();

    }

    public Tandem2mzidMzidContainer(String input, String dbFileFormatId,
                                    String msFileFormatId,
                                    Boolean isMs2SpecIdStartingAtZero,
                                    String decoyRegex,
                                    String proteinCodeRegex,
                                    MzIdentMLVersion ver)
            throws SAXException, ParserConfigurationException, IOException {
        this(input, dbFileFormatId, msFileFormatId, isMs2SpecIdStartingAtZero,
             decoyRegex, proteinCodeRegex, true, ver);
    }

    @Override
    public AnalysisCollection getAnalysisCollection() {

        SpectrumIdentification specIdent = new SpectrumIdentification();
        specIdent.setId(SPECT_IDENT_ID);
        specIdent.setSpectrumIdentificationList(siList);
        specIdent.setSpectrumIdentificationProtocol(siProtocol);
        List<SearchDatabaseRef> searchDbRefList = specIdent
                .getSearchDatabaseRef();
        SearchDatabaseRef searchDbRef = new SearchDatabaseRef();
        searchDbRef.setSearchDatabase(searchDb);
        searchDbRefList.add(searchDbRef);

        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
        InputSpectra inputSpec = new InputSpectra();
        inputSpec.setSpectraData(
                siList.getSpectrumIdentificationResult().get(0).getSpectraData());
        inputSpecList.add(inputSpec);
        AnalysisCollection analysisCollection = new AnalysisCollection();

        List<SpectrumIdentification> specIdentList = analysisCollection
                .getSpectrumIdentification();
        specIdentList.add(specIdent);

        return analysisCollection;
    }

    @Override
    public AnalysisProtocolCollection getAnalysisProtocolCollection() {

        AnalysisProtocolCollection analysisProtocolCollection
                = new AnalysisProtocolCollection();
        List<SpectrumIdentificationProtocol> sipList
                = analysisProtocolCollection.getSpectrumIdentificationProtocol();

        sipList.add(siProtocol);

        return analysisProtocolCollection;
    }

    @Override
    public AnalysisSampleCollection getAnalysisSampleCollection() {
        return null;
    }

    @Override
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();

        analysisSoftwareList.getAnalysisSoftware().add(analysisSoftwareXtandem);

        /*
         * TO DO - need to work out how to use Param CvParam cvParam = new
         * CvParam(); cvParam.setName("xtandem"); cvParam.setCvRef(psiCvID);
         * cvParam.setAccession("MS:1001476"); ParamAlternative paramAlt = new
         * ParamAlternative(); paramAlt.setCvParam(cvParam);
         *
         * analysisSoftware.setSoftwareName(makeCvParam("MS:1001476","xtandem",psiCV));
         * analysisSoftware.setSoftwareName(paramAlt);
         */
        return analysisSoftwareList;
    }

    @Override
    public AuditCollection getAuditCollection() {
        Organization org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName("myworkplace");
        org.getCvParam().add(MzidLibUtils.makeCvParam("MS:1000586",
                                                      "contact name",
                                                      CvConstants.PSI_CV,
                                                      "address"));

        //org.setAddress(address);
        List<Affiliation> affList = docOwner.getAffiliation();
        Affiliation aff = new Affiliation();
        aff.setOrganization(org);
        affList.add(aff);
        AuditCollection auditCollection = new AuditCollection();
        List<AbstractContact> contactList = auditCollection
                .getPersonOrOrganization();
        contactList.add(docOwner);
        contactList.add(org);

        return auditCollection;
    }

    @Override
    public BibliographicReference getBibliographicReference() {
        return null;
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
    public Inputs getInputs() {
        Inputs inputs = new Inputs();
        List<SearchDatabase> searchDbList = inputs.getSearchDatabase();

        searchDbList.add(searchDb);

        SourceFile sourceFile = new SourceFile();
        sourceFile.setLocation(inputFileName);
        sourceFile.setId(SOURCE_FILE_ID);

        FileFormat ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam("MS:1001401",
                                               "X\\!Tandem xml file",
                                               CvConstants.PSI_CV));
        sourceFile.setFileFormat(ff);
        List<SourceFile> sourceFileList = inputs.getSourceFile();
        sourceFileList.add(sourceFile);

        SpectraData spectraData = new SpectraData();

        SpectrumIDFormat sif = new SpectrumIDFormat();
        sif.setCvParam(MzidLibUtils.makeCvParam("MS:1000774",
                                                "multiple peak list nativeID format",
                                                CvConstants.PSI_CV));
        spectraData.setSpectrumIDFormat(sif);

        ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam(this.massSpecFileFormatId,
                                               this.massSpecFileFormatName,
                                               CvConstants.PSI_CV));
        spectraData.setFileFormat(ff);
        spectraData.setId(SPECTRA_DATA_ID);
        spectraData.setLocation(inputParams.getSpectrumPath());
        List<SpectraData> spectraDataList = inputs.getSpectraData();
        spectraDataList.add(spectraData);

        return inputs;
    }

    @Override
    public ProteinDetectionList getProteinDetectionList() {
        return null;
    }

    @Override
    public Provider getProvider() {
        Provider provider = new Provider();
        provider.setId("PROVIDER");

        ContactRole contactRole = new ContactRole();
        contactRole.setContact(docOwner);

        Role role = new Role();
        role.setCvParam(MzidLibUtils.makeCvParam("MS:1001271", "researcher",
                                                 CvConstants.PSI_CV));
        contactRole.setRole(role);

        provider.setContactRole(contactRole);

        return provider;
    }

    @Override
    public SequenceCollection getSequenceCollection() {
        return this.sequenceCollection;
    }

    @Override
    public SpectrumIdentificationList getSpectrumIdentificationList() {
        return this.siList;
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

        this.tandemParams = xfile.getPerformParameters();
        this.tandemVersion = tandemParams.getProcVersion();
        this.dbName = tandemParams.getSequenceSourceDescription_1();
        this.dbLocation = tandemParams.getSequenceSource_1();
        // databaseFileFormatId
        if (dbFileFormatId != null) {
            this.databaseFileFormatId = dbFileFormatId;
            this.databaseFileFormatName = getCvName(databaseFileFormatId);
        } else {
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
            this.isMs2SpectrumIdStartingAtZero
                    = this.massSpecFileFormatId.equalsIgnoreCase("MS:1000584");
        } else {
            this.isMs2SpectrumIdStartingAtZero = isMs2SpecIdStartingAtZero;
        }

        // MzIdentMLVersion
        if (null == ver) {
            this.version = MzIdentMLVersion.Version_1_1;
        } else {
            this.version = ver;
        }

        // docOwner
        docOwner = createDocOwner("firstName", "secondName", "address");

        // analysisSoftwareXtandem
        analysisSoftwareXtandem = createAnalysisSoftware("xtandem",
                                                         ANALYSIS_SOFT_ID,
                                                         CvConstants.XTANDEM,
                                                         this.tandemVersion);
    }

    private Person createDocOwner(String firstName, String secondName,
                                  String address) {
        Person person = new Person();
        person.setId("PERSON_DOC_OWNER");
        person.setFirstName(firstName);
        person.setLastName(secondName);
        person.getCvParam().add(MzidLibUtils.makeCvParam("MS:1000587",
                                                         "contact address",
                                                         CvConstants.PSI_CV,
                                                         address));

        return person;
    }

    private AnalysisSoftware createAnalysisSoftware(String name, String id,
                                                    CvParam cp,
                                                    String softVersion) {
        AnalysisSoftware analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName(name);
        analysisSoftware.setId(id);
        Param tempParam = new Param();
        tempParam.setParam(cp);
        analysisSoftware.setSoftwareName(tempParam);
        analysisSoftware.setVersion(softVersion);
        return analysisSoftware;
    }

    private boolean isMonoFragment(InputParams inputParams) {
        return "average".equals(inputParams.getSpectrumFragMassType());
    }

    private SearchModification translateToSearchModification(String reportedMod,
                                                             boolean fragmentIsMono,
                                                             boolean isFixedMod) {
        List<String> residues = new ArrayList<>();
        String[] temp = reportedMod.split("@");
        double monoMass = Double.parseDouble(temp[0]);

        residues.add(temp[1]);
        ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError,
                                             fragmentIsMono, residues);

        SearchModification searchMod = new SearchModification();
        searchMod.setFixedMod(isFixedMod);

        if (unimod == null) {
            searchMod.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001460",
                                                                "unknown modification",
                                                                CvConstants.PSI_CV));
        } else {

            searchMod.getCvParam().add(MzidLibUtils.makeCvParam("UNIMOD:"
                    + unimod.getRecordId(), unimod.getTitle(),
                                                                CvConstants.UNIMOD_CV));
        }
        searchMod.setMassDelta(Float.parseFloat(String.valueOf(monoMass)));

        residues.stream()
                .map((residue) -> {
                    if (residue.equals("[")) {
                        searchMod.getCvParam().add(MzidLibUtils
                                .makeCvParam("MS:1001189",
                                             "modification specificity N-term",
                                             CvConstants.PSI_CV));
                        residue = ".";
                    } else if (residue.equals("]")) {
                        searchMod.getCvParam().add(MzidLibUtils
                                .makeCvParam("MS:1001190",
                                             "modification specificity C-term",
                                             CvConstants.PSI_CV));
                        residue = ".";     //The any char must be inserted into mzid
                    }
                    return residue;
                })
                .forEach((residue) -> {
                    searchMod.getResidues().add(residue);
                });
        return searchMod;
    }

    private void createSpectrumIdentificationListAndSequenceCollection() {
        siList = new SpectrumIdentificationList();
        siList.setId(SI_LIST_ID);
        if (outputFragmentation) {
            FragmentationTable fragTbl = createFragmentationTable();
            siList.setFragmentationTable(fragTbl);
        }

        List<SpectrumIdentificationResult> specIdentResults = siList
                .getSpectrumIdentificationResult();
        sequenceCollection = new SequenceCollection();
        Map<String, Peptide> uniquePeps = new HashMap<>();
        Map<String, DBSequence> foundProts = new HashMap<>();
        Map<String, PeptideEvidence> pepEvidLookup = new HashMap<>();
        List<Peptide> peptideList = sequenceCollection.getPeptide();
        SpectraData spectraData = new SpectraData();
        spectraData.setId(SPECTRA_DATA_ID);

        int sirCounter = 1; //Counter used to create unique ID for SpectrumIdentificationResult

        // Iterate over all the spectra
        @SuppressWarnings("unchecked")
        Iterator<Spectrum> iter = this.xfile.getSpectraIterator();
        while (iter.hasNext()) {
            // Get the next spectrum.
            Spectrum spectrum = iter.next();
            int spectrumNumber = spectrum.getSpectrumNumber();
            //note: spectrum number seems to be a sequential index. 
            //For the spectrum number as found in xtandem file use spectrumId

            /*
             * ***********************************************
             *  *** Setup SpectrumIdentificationResult ****
             * *********************************************
             */
            SpectrumIdentificationResult specIdentRes
                    = new SpectrumIdentificationResult();

            // Get the peptide hits.
            List<de.proteinms.xtandemparser.xtandem.Peptide> pepList
                    = this.xfile.getPeptideMap().getAllPeptides(spectrumNumber);

            //int pepEvidCounter = 1;
            int siiCounter = 1; //Counter used to create unique ID for SpectrumIdentificationItem

            SpectrumIdentificationItem sii = null;
            Peptide mzidPep = null;

            List<IonType> ionTypeList = null;

            Map<String, SpectrumIdentificationItem> siiMap = new HashMap<>();

            //This list can contain both different pep2Protein maps and
            for (de.proteinms.xtandemparser.xtandem.Peptide peptide : pepList) {

                for (Domain domain : peptide.getDomains()) {

                    //In mzIdentML we have 1 SepctrumIdentificationItem (SII) 
                    //linked to 1 Peptide item via the peptide_ref attribute. 
                    //Each Peptide item is a unique combination of:
                    //peptidesequence + modifications + substitutionModifications.
                    String uniquePepKey = getPeptideKey(domain, this.xfile);
                    //If it is a new global peptide, 
                    //then initialize a new mzIdentML Peptide object: 
                    if (!uniquePeps.containsKey(uniquePepKey)) {
                        mzidPep = new uk.ac.ebi.jmzidml.model.mzidml.Peptide();
                        mzidPep.setPeptideSequence(domain.getDomainSequence());
                        mzidPep.setId(uniquePepKey);

                        //Parse the modifications and add them to the mzidPep object:
                        parseModificationsAndSubstitutions(mzidPep, domain,
                                                           this.xfile,
                                                           fragmentIsMono);

                        if (!uniquePeps.containsKey(uniquePepKey)) {
                            peptideList.add(mzidPep);
                            uniquePeps.put(uniquePepKey, mzidPep);
                        }
                    } else {
                        //otherwise fetch it from the map:
                        mzidPep = uniquePeps.get(uniquePepKey);
                    }

                    //Here we have to decide whether to create a new 
                    //SepctrumIdentificationItem (SII) item OR if this new 
                    //domain item is just what we call in mzIdentML a 
                    //new PeptideEvidence (i.e. the same Peptide item found in
                    //another part of the protein sequence or even in another protein).
                    String siiKey = getSiiKey(domain, this.xfile);
                    if (isNewSii(siiKey, siiMap)) {
                        /*
                         ****************************************************
                         ****** Create new SpectrumIdentificationItem *******
                         * ***************************************************
                         */
                        sii = new SpectrumIdentificationItem();

                        //From XSD:
                        //Set to true if the producers of the file has deemed 
                        //that the identification has passed a given threshold 
                        //or been validated as correct. If no such threshold 
                        //has been set, value of true should be given for all results.
                        //Since this is not something related to X!Tandem 
                        //functionality (but rather expected in a post-X!Tandem 
                        //processing QC type of step), leave it as true always:
                        //TODO - check again if X!Tandem does not have any 
                        //post-processing built-in which results in an output 
                        //file that has items above and below a certain threshold... 
                        sii.setPassThreshold(true);

                        //add sii to sir:
                        specIdentRes.getSpectrumIdentificationItem().add(sii);
                        siiMap.put(siiKey, sii);

                        if (outputFragmentation) {
                            //ionTypeList = sii.getFragmentation();
                            Fragmentation frag = new Fragmentation();
                            sii.setFragmentation(frag);
                            ionTypeList = frag.getIonType();
                            parseFragmentationData(ionTypeList, domain, peptide,
                                                   this.xfile,
                                                   mzMeasure, intMeasure,
                                                   errorMeasure);

                            //if ionTypeList is still empty after parsing above, throw error:
                            if (ionTypeList.isEmpty()) {
                                throw new RuntimeException(
                                        "Error while parsing peptide ["
                                        + peptide.getDomains().get(0)
                                        .getDomainSequence()
                                        + "] identification in spectrum ["
                                        + spectrum.getSpectrumId()
                                        + "]: no fragmentation data found");
                            }
                        }

                        sii.setPeptide(mzidPep);

                        parseScoresAndOtherSiiAttributes(sii, domain, spectrum);
                        sii.setId("SII_" + sirCounter + "_" + siiCounter);
                        siiCounter++;
                    } else {
                        sii = siiMap.get(siiKey);
                    }

                    //Parse protein details into DBSequence objects:
                    String[] protAccessionAndDescription = parseProteinDetails(
                            foundProts, domain, peptide, this.xfile);
                    String protAccession = protAccessionAndDescription[0];
                    String protDescription = protAccessionAndDescription[1];

                    String pepEvidKey = getGlobalPeptideEvidenceKey(domain,
                                                                    protAccession,
                                                                    uniquePepKey);
                    PeptideEvidence pepEvid = null;
                    if (pepEvidLookup.get(pepEvidKey) == null) {
                        //Is new evidence, so initialize a new PeptideEvidence object:
                        pepEvid = parseNewPeptideEvidence(domain);
                        //store in map:
                        pepEvidLookup.put(pepEvidKey, pepEvid);
                        //extra details:
                        pepEvid.setId("PE" + sirCounter + "_" + siiCounter + "_"
                                + pepEvidCounter);
                        pepEvidCounter++;

                        //Decoy check part:
                        pepEvid.setIsDecoy(Boolean.FALSE);
                        if (this.decoyRegularExpression != null) {
                            if (protAccession.matches(
                                    this.decoyRegularExpression)
                                    || protAccession.split(
                                            this.decoyRegularExpression).length
                                    > 1) {
                                pepEvid.setIsDecoy(Boolean.TRUE);
                            }
                        }
                        //if protein description ends with ":reversed", then 
                        //peptide is decoy. Example:
                        //<protein expect="-1.4" id="5736.1" uid="4657" 
                        //label="tr|F1RRZ6|F1RRZ6_PIG Uncharacterized protein 
                        //(Fragment) OS=Sus scrofa GN=SH3PXD2B..." sumI="3.78" >
                        //  <note label="description">tr|F1RRZ6|F1RRZ6_PIG 
                        //Uncharacterized protein (Fragment) OS=Sus scrofa 
                        //GN=SH3PXD2B PE=4 SV=1:reversed</note>
                        if (protDescription.endsWith(XTANDEM_REVERSED_FLAG)) {
                            pepEvid.setIsDecoy(Boolean.TRUE);
                        }

                        //link it to our sequenceCollection list:
                        sequenceCollection.getPeptideEvidence().add(pepEvid);
                        //link to mzidPep:
                        pepEvid.setPeptide(mzidPep);
                        //link to DBSequence:
                        pepEvid.setDBSequence(foundProts.get(protAccession));

                    } else {
                        //otherwise fetch it from the map:
                        pepEvid = pepEvidLookup.get(pepEvidKey);
                    }

                    PeptideEvidenceRef peptideEvidenceRef
                            = new PeptideEvidenceRef();
                    peptideEvidenceRef.setPeptideEvidence(pepEvid);
                    //link to sii:
                    sii.getPeptideEvidenceRef().add(peptideEvidenceRef);
                }

            }

            /*
             * ***********************************************
             *  *** Complete SpectrumIdentificationResult ****
             * *********************************************
             */
            //setting the sectrumID so that it conforms to the specifications of 
            //the controlled vocabulary item MS:1000774 where spectrumID should start from 0. 
            //The problem is that we don't know for sure from the X!Tandem file alone what is the 
            //type of the spectrum file...below we have now a very simple check for whether
            //the spectrum file submitted to X!Tandem had ids starting from 0 or from 1.
            //For example,  int MZML format the spectrum id in xtandem output already starts from 0.
            int spectrumId;
            if (this.isMs2SpectrumIdStartingAtZero) {
                spectrumId = (spectrum.getSpectrumId());
            } else {
                //then we assume it was starting from 1:
                spectrumId = (spectrum.getSpectrumId() - 1);
            }

            //TODO add checks to find out what is the data format submitted to 
            //X!Tandem and to correct the spectrumID value accordingly            
            specIdentRes.setSpectrumID("index=" + spectrumId);
            specIdentRes.setSpectraData(spectraData);
            specIdentRes.setId("SIR_" + sirCounter);
            // Get the support data for each spectrum.
            SupportData supportData = this.xfile.getSupportData(spectrumNumber);
            String label = supportData.getFragIonSpectrumDescription();
            if (label != null && !label.equals("")) {
                List<CvParam> sirCvParamList = specIdentRes.getCvParam();
                CvParam cvp = MzidLibUtils.makeCvParam("MS:1000796",
                                                       "spectrum title",
                                                       CvConstants.PSI_CV,
                                                       label);
                sirCvParamList.add(cvp);

            }

            specIdentResults.add(specIdentRes);
            sirCounter++;

            //TO - currently only implements the case where the same peptide 
            //matches to different proteins
            // Initialize the array lists
            //ArrayList<double> mzValues;
            //ArrayList<double> intensityValues;
            // Get the spectrum fragment mz and intensity values
            //mzValues = supportData.getXValuesFragIonMass2Charge();
            //intensityValues = supportData.getYValuesFragIonMass2Charge();
            // Fill the maps
            //allMzValues.put(new Integer(spectrumNumber), mzValues);
            // allIntensityValues.put(new Integer(spectrumNumber), intensityValues);
        }

    }

    private FragmentationTable createFragmentationTable() {
        /*
         * <Measure id="m_mz"> <cvParam cvRef="PSI-MS" accession="MS:1001225"
         * name="product ion m/z"/> </Measure> <Measure id="m_intensity">
         * <cvParam cvRef="PSI-MS" accession="MS:1001226" name="product ion
         * intensity"/> </Measure> <Measure id="m_error"> <cvParam
         * cvRef="PSI-MS" accession="MS:1001227" name="product ion m/z error"
         * unitAccession="MS:1000040" unitName="m/z" unitCvRef="PSI-MS"/>
         * </Measure>
         */

        FragmentationTable fragTable = new FragmentationTable();
        List<Measure> measureList = fragTable.getMeasure();
        measureList.add(mzMeasure);
        measureList.add(intMeasure);
        measureList.add(errorMeasure);

        return fragTable;
    }

    private boolean isNewSii(String siiKey,
                             Map<String, SpectrumIdentificationItem> siiMap) {
        return siiMap.get(siiKey) == null;
    }

    private String getPeptideKey(Domain domain, XTandemFile ixtandemFile) {
        //is really the same as in getSIIKey, but the context of both maps is 
        //different (peptide map is global and siimap is local within 
        //a SepctrumIdentificationResult ):
        return getSiiKey(domain, ixtandemFile);
    }

    private String getSiiKey(Domain domain, XTandemFile ixtandemFile) {
        List<de.proteinms.xtandemparser.interfaces.Modification> fixModList
                = ixtandemFile.getModificationMap().getFixedModifications(
                        domain.getDomainKey());
        List<de.proteinms.xtandemparser.interfaces.Modification> varModList
                = ixtandemFile.getModificationMap().getVariableModifications(
                        domain.getDomainKey());

        StringBuilder fixMods = new StringBuilder();

        for (de.proteinms.xtandemparser.interfaces.Modification fixMod
                : fixModList) {
            String name = fixMod.getName();
            if (fixMod.isSubstitution()) {
                name += "_subs_" + fixMod.getSubstitutedAminoAcid();
            }
            int loc = Integer.parseInt(fixMod.getLocation());
            int pepLoc = loc - domain.getDomainStart() + 1;
            fixMods.append(name).append("$").append(pepLoc).append(";");
        }

        StringBuilder varMods = new StringBuilder();
        for (de.proteinms.xtandemparser.interfaces.Modification varMod
                : varModList) {
            String name = varMod.getName();
            if (varMod.isSubstitution()) {
                name += "_subs_" + varMod.getSubstitutedAminoAcid();
            }

            int loc = Integer.parseInt(varMod.getLocation());
            int pepLoc = loc - domain.getDomainStart() + 1;
            varMods.append(name).append("$").append(pepLoc).append(";");
        }
        String siiKey = domain.getDomainSequence() + "_" + varMods.toString()
                + "_" + fixMods.toString() + "_";
        return siiKey;
    }

    private void parseModificationsAndSubstitutions(Peptide mzidPep,
                                                    Domain domain,
                                                    XTandemFile ixtandemFile,
                                                    boolean fragmentIsMono) {
        //Parse the modifications
        List<de.proteinms.xtandemparser.interfaces.Modification> fixModList
                = ixtandemFile.getModificationMap().getFixedModifications(
                        domain.getDomainKey());
        List<de.proteinms.xtandemparser.interfaces.Modification> varModList
                = ixtandemFile.getModificationMap().getVariableModifications(
                        domain.getDomainKey());

        fixModList.stream()
                .forEach((reportedMod) -> {
                    if (!reportedMod.isSubstitution()) {
                        Modification mzidMod = MzidLibUtils
                                .translateToMzidModification(
                                        reportedMod,
                                        domain,
                                        fragmentIsMono, unimodMassError);
                        mzidPep.getModification().add(mzidMod);
                    } else {
                        SubstitutionModification mzidSubs = MzidLibUtils
                                .translateToMzidSubstitution(
                                        reportedMod, domain, fragmentIsMono);
                        mzidPep.getSubstitutionModification().add(mzidSubs);
                    }
                });

        varModList.stream()
                .forEach((reportedMod) -> {
                    if (!reportedMod.isSubstitution()) {
                        Modification mzidMod = MzidLibUtils
                                .translateToMzidModification(
                                        reportedMod,
                                        domain,
                                        fragmentIsMono, unimodMassError);
                        mzidPep.getModification().add(mzidMod);
                    } else {
                        SubstitutionModification mzidSubs = MzidLibUtils
                                .translateToMzidSubstitution(
                                        reportedMod, domain, fragmentIsMono);
                        mzidPep.getSubstitutionModification().add(mzidSubs);
                    }
                });

    }

    private void parseFragmentationData(List<IonType> ionTypeList, Domain domain,
                                        de.proteinms.xtandemparser.xtandem.Peptide peptide,
                                        XTandemFile ixtandemFile,
                                        Measure mzMeasure, Measure intMeasure,
                                        Measure errorMeasure) {
        // Get the fragment ions
        if (outputFragmentation) {
            try {

                //Gets the list with a number of arrays: 
                //b ions, one of y ions, b ions - NH3, x ions, a ions, b ions - H2O etc
                //NB: The underlying de.proteinms.xtandemparser simulates the peptide 
                //fragmentation (in silico) again as well as the spectrum matching in order 
                //to (re)generate the list of ions (which is apparently then not present in 
                //the X!Tandem output file). The in silico fragmentation of this underlying 
                //library makes indeed theoretical fragments of different charges, 
                //more specifically when the precursor ion has charge X it generates 
                //fragment ions of charges 1 to X (the rationale of the authors here 
                //must have been that fragment ions can lose charge during fragmentation, 
                //but never gain...not sure if this is always the case). 
                System.out.println(
                        "WARN: Fragment ion annotation data is inferred based on new calculations "
                        + "triggered by this (Tandem2mzid) conversion library. "
                        + "This is done because fragment ion annotation information is not present "
                        + "in X!Tandem file (X!Tandem SLEDGEHAMMER (2013.09.01) and previous) . ");
                //TODO add real logger option
                @SuppressWarnings("unchecked")
                List<FragmentIon[]> ionsForPeptide = ixtandemFile
                        .getFragmentIonsForPeptide(peptide, domain, ixtandemFile
                                                   .getInputParameters()
                                                   .getSpectrumMonoIsoMassError());
                //TODO improvement for item above: a good trade-off would be that the X!Tandem 
                //parser lib reads the search parameters (present in the X!Tandem output file) 
                //and use in it's in silico fragmentation only the same ion types (i.e. b and y) 
                //used in the search.

                /*
                 * <IonType index="2 4 4 9 7 10 8 11 8 13" charge="1"> <cvParam
                 * cvRef="PSI-MS" accession="MS:1001366" name="frag: internal ya
                 * ion"/> <FragmentArray values="286 644.2 329.9 329.9 514.2 "
                 * Measure_ref="m_mz"/> <FragmentArray values="32 194 2053 2053
                 * 125" Measure_ref="m_intensity"/> <FragmentArray
                 * values="-0.2125 -0.1151 -0.2772 -0.2772 -0.0620"
                 * Measure_ref="m_error"/> </IonType>
                 *
                 */
                //The part below creates a structure like the one commented above, 
                //one structure <IonType> per charge found. The fragments reported
                //in each structure are only the ones reported on the respective charge.
                Map<Integer, Map<Integer, List<FragmentIon>>> map
                        = orderFragmentIons(ionsForPeptide);
                map.entrySet().stream()
                        .forEach((entryCharge) -> {
                            entryCharge.getValue().entrySet().stream()
                                    .filter((entryType)
                                            -> !(entryType.getValue().isEmpty()))
                                    .forEach((entryType) -> {
                                        IonType ion = new IonType();
                                        ion.setCharge(entryCharge.getKey());
                                        ion.setCvParam(MzidLibUtils
                                                .getFragmentCVParam(entryType
                                                        .getKey()));
                                        List<Integer> ionIndexList = ion
                                                .getIndex();

                                        FragmentArray mzArray
                                                = new FragmentArray();
                                        FragmentArray intArray
                                                = new FragmentArray();
                                        FragmentArray errorArray
                                                = new FragmentArray();
                                        mzArray.setMeasure(mzMeasure);
                                        intArray.setMeasure(intMeasure);
                                        errorArray.setMeasure(errorMeasure);
                                        List<Float> mzValues = mzArray
                                                .getValues();
                                        List<Float> intValues = intArray
                                                .getValues();
                                        List<Float> errorValues = errorArray
                                                .getValues();
                                        List<FragmentArray> fragmentList = ion
                                                .getFragmentArray();
                                        fragmentList.add(mzArray);
                                        fragmentList.add(intArray);
                                        fragmentList.add(errorArray);

                                        ionTypeList.add(ion);

                                        entryType.getValue().stream()
                                                .map((fragIon) -> {
                                                    //Reported MZ is the theoretical value
                                                    mzValues.add(
                                                            (float) (fragIon
                                                            .getMZ() + fragIon
                                                            .getTheoreticalExperimentalMassError()));
                                                    return fragIon;
                                                })
                                                .map((fragIon) -> {
                                                    intValues.add(
                                                            //Note intensity values in Tandem 
                                                            //do not match the source spectrum, 
                                                            //appears that some processing happens                                                    
                                                            (float) fragIon
                                                            .getIntensity());
                                                    return fragIon;
                                                })
                                                .map((fragIon) -> {
                                                    errorValues.add(
                                                            (float) fragIon
                                                            .getTheoreticalExperimentalMassError());
                                                    return fragIon;
                                                })
                                                .forEach((fragIon) -> {
                                                    ionIndexList.add(fragIon
                                                            .getNumber());  //index position
                                                });
                                    });
                        });

            } catch (Exception e) {
                throw new RuntimeException(
                        "Error while parsing the MS/MS fragmentation data. Please check if "
                        + " your input files indeed contain fragmentation data "
                        + "(this is optional and depending "
                        + "on your X!Tandem settings it could be absent). Error details: ",
                        e);
            }
        }
    }

    private Map<Integer, Map<Integer, List<FragmentIon>>> orderFragmentIons(
            List<FragmentIon[]> ionsForPeptide) {
        Map<Integer, Map<Integer, List<FragmentIon>>> map = new HashMap<>();
        ionsForPeptide.stream()
                .forEach((ions) -> {
                    for (FragmentIon ion : ions) {
                        Map<Integer, List<FragmentIon>> mapCharge = map.get(
                                (int) ion
                                .getCharge());
                        if (mapCharge == null) {
                            mapCharge = new HashMap<>();
                            map.put((int) ion.getCharge(), mapCharge);
                        }
                        List<FragmentIon> list = mapCharge.get(ion.getType());
                        if (list == null) {
                            list = new ArrayList<>();
                            mapCharge.put(ion.getType(), list);
                        }
                        list.add(ion);
                    }
                });
        return map;
    }

    private void parseScoresAndOtherSiiAttributes(SpectrumIdentificationItem sii,
                                                  Domain domain,
                                                  Spectrum spectrum) {

        //Get precursorMh reported by X!Tandem and convert it back
        //to m/z using:
        //     m/z=(mh + z*1.007276 - 1*1.007276)/z
        //where mh is M+H    and H=1.007276466812 (proton mass rounded from 1.007276466812)
        //  X!Tandem is using H=1.007276 (proton mass rounded from 1.007276466812) 
        //so we will use that
        double protonMass = 1.007276;
        int precursorCharge = spectrum.getPrecursorCharge();
        double expMz = (spectrum.getPrecursorMh() + precursorCharge * protonMass
                - protonMass) / precursorCharge;
        //round at 6 decimals:
        expMz = Utils.round(expMz, 6);
        sii.setExperimentalMassToCharge(expMz);

        double calcMz = (domain.getDomainMh() + precursorCharge * protonMass
                - protonMass) / precursorCharge;
        calcMz = Utils.round(calcMz, 6);
        sii.setCalculatedMassToCharge(calcMz);
        sii.setChargeState(precursorCharge);
        sii.setRank(1);

        List<CvParam> cvParamList = sii.getCvParam();
        //<cvParam accession="MS:1001330" name="xtandem:expect" cvRef="PSI-MS"  value="1.1e-003" />
        //<cvParam accession="MS:1001331" name="xtandem:hyperscore" cvRef="PSI-MS"  value="60.4" />
        double evalue = domain.getDomainExpect();
        double hyperscore = domain.getDomainHyperScore();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001330",
                                                 "X\\!Tandem:expect",
                                                 CvConstants.PSI_CV, ""
                                                 + evalue));
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001331",
                                                 "X\\!Tandem:hyperscore",
                                                 CvConstants.PSI_CV,
                                                 "" + hyperscore));
    }

    private String[] parseProteinDetails(Map<String, DBSequence> foundProts,
                                         Domain domain,
                                         de.proteinms.xtandemparser.xtandem.Peptide peptide,
                                         XTandemFile ixtandemFile) {
        Protein protein = ixtandemFile.getProteinMap().getProtein(domain
                .getProteinKey());

        String protAccession = "";
        String protDescription = "";
        String protSeq = "";
        //int protLength;

        if (protein != null) {
            protAccession = protein.getLabel();
            //If the proteinCodeRegexPattern is set, then extract the proteinAccession 
            //according to regex:
            if (this.proteinCodeRegexPattern != null) {
                Matcher matcher = this.proteinCodeRegexPattern.matcher(
                        protAccession);
                //If there is a match, then extract this as the new proteinAccession value:
                if (matcher.find()) {
                    protAccession = matcher.group();
                }
            }
            protDescription = protein.getDescription();
            //System.out.println("prot: " + protAccession);
            protSeq = peptide.getSequence();        //getSequence returns the protein sequence
            protSeq = protSeq.replaceAll("\\s+", "");
        } else {
            throw new RuntimeException(
                    "Unexpected problem: protein not found in parsed protMap");
        }

        //Use Hash map to test if Protein sequence has been added to DBSeq before
        if (!foundProts.containsKey(protAccession)) {
            DBSequence dbSeq = new DBSequence();
            foundProts.put(protAccession, dbSeq);
            dbSeq.setAccession(protAccession);
            dbSeq.setName(protDescription);
            dbSeq.setSeq(protSeq);
            dbSeq.setLength(protSeq.length());
            dbSeq.setId("dbseq_" + protAccession);
            dbSeq.setSearchDatabase(searchDb);
            //dbSeq.setSearchDatabase(searchDB);
            sequenceCollection.getDBSequence().add(dbSeq);
        }
        String[] result = {protAccession, protDescription};
        return result;
    }

    /**
     * The method returns the key to check whether the domain represents a new
     * peptide evidence, which means in mzIdentML terms: the combination of
     * (peptidesequence + modifications + substitutionModifications) +
     * proteinAccession + proteinLocation(i.e start,end)is unique so far.
     * This will mean that we have to make a new mzIdentML PeptideEvidence
     * object.
     *
     * @param domain        domain
     *
     * @param protAccession protein accession
     * @param peptideKey    : as returned by getPeptideKey method
     *
     *
     * @return key for global peptide evidence
     */
    private String getGlobalPeptideEvidenceKey(Domain domain,
                                               String protAccession,
                                               String peptideKey) {
        int end = domain.getDomainEnd();
        int start = domain.getDomainStart();
        String uniqueProtLocation = protAccession + "_" + start + "_" + end;

        String testPepMods = uniqueProtLocation + "_" + peptideKey;
        return testPepMods;
    }

    /**
     * Parses the X!Tandem details in to a new mzIdentML PeptideEvidence object.
     *
     * @param domain domain
     *
     * @return PeptideEvidence object
     */
    private PeptideEvidence parseNewPeptideEvidence(Domain domain) {
        PeptideEvidence pepEvid = new PeptideEvidence();
        pepEvid.setEnd(domain.getDomainEnd());
        pepEvid.setStart(domain.getDomainStart());

        //pepEvid.setMissedCleavages(domain.getMissedCleavages());
        char post = domain.getDownFlankSequence().charAt(0);
        if (post == ']') {
            post = '-';

        }
        pepEvid.setPost("" + post); //Reports 4 chars, we need the first only

        char pre = domain.getUpFlankSequence().charAt(domain
                .getUpFlankSequence().length() - 1);
        if (pre == '[') {
            pre = '-';
        }
        pepEvid.setPre("" + pre);    //Reports 4 chars, we need last only
        return pepEvid;
    }

    private SpectrumIdentificationProtocol createSpectrumIdentificationProtocol() {
        SpectrumIdentificationProtocol siProtocolRet
                = new SpectrumIdentificationProtocol();
        siProtocolRet.setId(SI_PROTOCOL_ID);
        siProtocolRet.setAnalysisSoftware(analysisSoftwareXtandem);

        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
        //siProtocolRet.setSearchType(makeCvParam("MS:1001083","ms-ms search",psiCV));
        Param tempParam = new Param();
        tempParam.setParam(CvConstants.MS_MS_SEARCH);
        siProtocolRet.setSearchType(tempParam);

        ParamList paramList = siProtocolRet.getAdditionalSearchParams();
        if (paramList == null) {
            paramList = new ParamList();
            siProtocolRet.setAdditionalSearchParams(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();

        //does not appear to be a way in Tandem of specifying parent mass is average
        boolean parentIsMono = true;
        if (parentIsMono) {
            cvParamList.add(CvConstants.PARENT_MASS_TYPE_MONO);
        } else {
            cvParamList.add(CvConstants.PARENT_MASS_TYPE_AVERAGE);
        }

        if (fragmentIsMono) {
            cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_MONO);

        } else {
            cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_AVERAGE);
        }

        // Add "no special processing" cv term if this is mzid 1.2 version
        if (this.version.equals(MzIdentMLVersion.Version_1_2)) {
            cvParamList.add(CvConstants.NO_SPECIAL_PROCESSING);
        }

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams
                .getSearchModification();

        //residue, potential modification mass
        String varMods = inputParams.getResiduePotModMass();
        String fixedMods = inputParams.getResidueModMass();

        if (varMods != null) {
            String[] allVarMods = varMods.split(",");

            //<note type="input" label="residue, modification mass">57.021469@C</note>
            // <note type="input" label="residue, potential modification mass">15.99492@M</note>
            for (String varMod : allVarMods) {

                //If there are no modification parameters, 
                //this is identified with the value 'None'. 
                //In this case, skip the block below:
                if (!varMod.equalsIgnoreCase("None")) {
                    SearchModification searchMod
                            = translateToSearchModification(varMod,
                                                            fragmentIsMono,
                                                            false);
                    searchModList.add(searchMod);
                }
            }
        }
        if (fixedMods != null) {
            String[] allFixedMods = fixedMods.split(",");

            for (String fixedMod : allFixedMods) {
                //If there are no modification parameters, 
                //this is identified with the value 'None'. 
                //In this case, skip the block below:
                if (!fixedMod.equalsIgnoreCase("None")) {
                    SearchModification searchMod
                            = translateToSearchModification(fixedMod,
                                                            fragmentIsMono, true);
                    searchModList.add(searchMod);
                }
            }
        }

        /*
         * <ModificationParams> <SearchModification fixedMod="false"> <ModParam
         * massDelta="15.994919" residues="M"> <cvParam accession="UNIMOD:35"
         * name="Oxidation" cvRef="UNIMOD" /> </ModParam> </SearchModification>
         */
        //Only add this group if there are any modifications in the list:
        if (searchModList.size() > 0) {
            siProtocolRet.setModificationParams(modParams);
        }

        /*
         * <Enzymes independent="0"> <Enzyme id="ENZ_1" CTermGain="OH"
         * NTermGain="H" missedCleavages="1" semiSpecific="0"> <EnzymeName>
         * <cvParam accession="MS:1001251" name="Trypsin" cvRef="PSI-MS" />
         * </EnzymeName> </Enzyme> </Enzymes>
         */
        Enzymes enzymes = siProtocolRet.getEnzymes();

        if (enzymes == null) {
            enzymes = new Enzymes();
            siProtocolRet.setEnzymes(enzymes);
        }

        enzymes.setIndependent(false);

        List<Enzyme> enzymeList = enzymes.getEnzyme();

        Enzyme enzyme = Utils.getXtandemEnzyme(inputParams
                .getProteinCleavageSite(), inputParams
                                               .getScoringMissCleavageSites());

        enzymeList.add(enzyme);

        boolean isDaltons;

        if (inputParams.getSpectrumMonoIsoMassErrorUnits().equalsIgnoreCase(
                "ppm")) {
            isDaltons = false;
        } else {
            isDaltons = true;
            //Dynamically set Unimod lookup mass error
            unimodMassError = inputParams.getSpectrumMonoIsoMassError();
        }

        /*
         * <FragmentTolerance> <cvParam accession="MS:1001412" name="search
         * tolerance plus value" value="0.5" cvRef="PSI-MS"
         * unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
         * <cvParam accession="MS:1001413" name="search tolerance minus value"
         * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
         * unitName="dalton" unitCvRef="UO" /> </FragmentTolerance>
         * <ParentTolerance> <cvParam accession="MS:1001412" name="search
         * tolerance plus value" value="0.5" cvRef="PSI-MS"
         * unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
         * <cvParam accession="MS:1001413" name="search tolerance minus value"
         * value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
         * unitName="dalton" unitCvRef="UO" /> </ParentTolerance>
         */
        CvParam fragCvPlus = MzidLibUtils.getCvParamWithMassUnits(isDaltons);
        CvParam fragCvMinus = MzidLibUtils.getCvParamWithMassUnits(isDaltons);
        fragCvPlus.setAccession("MS:1001412");
        fragCvPlus.setName("search tolerance plus value");
        fragCvMinus.setAccession("MS:1001413");
        fragCvMinus.setName("search tolerance minus value");
        fragCvPlus.setValue("" + inputParams.getSpectrumMonoIsoMassError());
        fragCvMinus.setValue("" + inputParams.getSpectrumMonoIsoMassError());
        Tolerance fragTol = new Tolerance();
        List<CvParam> fragCvList = fragTol.getCvParam();
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        boolean isDaltonsParent = !inputParams
                .getSpectrumParentMonoIsoMassErrorUnits()
                .equalsIgnoreCase("ppm");

        CvParam parCvPlus = MzidLibUtils
                .getCvParamWithMassUnits(isDaltonsParent);
        CvParam parCvMinus = MzidLibUtils.getCvParamWithMassUnits(
                isDaltonsParent);
        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue("" + inputParams
                .getSpectrumParentMonoIsoMassErrorPlus());
        parCvMinus.setValue("" + inputParams
                .getSpectrumParentMonoIsoMassErrorMinus());

        Tolerance parTol = new Tolerance();
        List<CvParam> parCvList = parTol.getCvParam();
        parCvList.add(parCvPlus);
        parCvList.add(parCvMinus);

        siProtocolRet.setFragmentTolerance(fragTol);
        siProtocolRet.setParentTolerance(parTol);

        ParamList sipParamList = siProtocolRet.getThreshold();
        if (sipParamList == null) {
            sipParamList = new ParamList();
            siProtocolRet.setThreshold(sipParamList);
        }
        cvParamList = sipParamList.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001494", "no threshold",
                                                 CvConstants.PSI_CV));
        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />

        return siProtocolRet;
    }

    private SearchDatabase createSearchDatabase() {
        SearchDatabase searchDbRet = new SearchDatabase();
        searchDbRet.setId(SEARCH_DB_ID);
        long numProts = this.tandemParams.getTotalProteinsUsed();
        searchDbRet.setNumDatabaseSequences(numProts);
        //<cvParam accession="MS:1001401" name="xtandem xml file" cvRef="PSI-MS"/>

        UserParam param = new UserParam();
        param.setName(this.dbName);
        Param tempParam = new Param();
        tempParam.setParam(param);
        searchDbRet.setDatabaseName(tempParam);

        searchDbRet.setLocation(this.dbLocation);

        FileFormat ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam(this.databaseFileFormatId,
                                               this.databaseFileFormatName,
                                               CvConstants.PSI_CV));
        searchDbRet.setFileFormat(ff);

        return searchDbRet;
    }

}
