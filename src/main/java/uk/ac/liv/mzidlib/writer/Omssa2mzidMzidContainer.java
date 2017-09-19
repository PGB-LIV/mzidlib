/*
 * Date: 02-Sep-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.Omssa2mzidMzidContainer.java
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import de.proteinms.omxparser.OmssaOmxFile;
import de.proteinms.omxparser.util.MSHitSet;
import de.proteinms.omxparser.util.MSHitSet_hits;
import de.proteinms.omxparser.util.MSHits;
import de.proteinms.omxparser.util.MSHits_mods;
import de.proteinms.omxparser.util.MSHits_mzhits;
import de.proteinms.omxparser.util.MSHits_pephits;
import de.proteinms.omxparser.util.MSMZHit;
import de.proteinms.omxparser.util.MSMZHit_ion;
import de.proteinms.omxparser.util.MSModHit;
import de.proteinms.omxparser.util.MSModHit_modtype;
import de.proteinms.omxparser.util.MSPepHit;
import de.proteinms.omxparser.util.MSRequest;
import de.proteinms.omxparser.util.MSRequest_settings;
import de.proteinms.omxparser.util.MSResponse;
import de.proteinms.omxparser.util.MSSearch;
import de.proteinms.omxparser.util.MSSearchSettings;
import de.proteinms.omxparser.util.MSSearch_request;
import de.proteinms.omxparser.util.MSSearch_response;
import de.proteinms.omxparser.util.MSSpectrum;
import de.proteinms.omxparser.util.MSSpectrum_abundance;
import de.proteinms.omxparser.util.MSSpectrum_ids;
import de.proteinms.omxparser.util.MSSpectrum_mz;
import de.proteinms.omxparser.util.OmssaEnumerators;
import de.proteinms.omxparser.util.OmssaModification;
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
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.ReadUnimod;
import uk.ac.liv.mzidlib.util.MzidLibUtils;
import uk.ac.liv.unimod.ModT;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 02-Sep-2017 00:12:41
 */
public class Omssa2mzidMzidContainer implements MzidContainer {

    private String inputOmssaFileName;
    private String modsFile;
    private String userModsFile;
    private MzIdentMLVersion version;
    private OmssaOmxFile omxFile;
    private Boolean outputFragmentation;
    private String decoyRegularExpression;
    private AnalysisSoftware analysisSoftwareOmssa;
    private SequenceCollection sequenceCollection;
    private SpectrumIdentificationList spectrumIdentList;
    private SearchDatabase searchDatabase;
    private FragmentationTable fragmentationTable;
    private SpectraData spectraData;
    private SpectrumIdentificationProtocol spectrumIdentificationProtocol;

    private Person docOwner;

    private static final String ANALYSIS_SOFT_ID = "ID_software";
    private static final String SI_PROTOCOL_ID = "SearchProtocol_1";
    private static final String SOURCE_FILE_ID = "SourceFile_1";
    private static final String SEARCH_DB_ID = "SearchDB_1";
    private static final String SPECTRA_DATA_ID = "SID_1";
    private static final String SPEC_IDENT_ID = "SpecIdent_1";
    private static final String SI_LIST_ID = "SI_List_1";
    private static final String MEASURE_MZ_ID = "Measure_MZ";
    private static final String MEASURE_INT_ID = "Measure_Int";
    private static final String MEASURE_ERROR_ID = "Measure_Error";
    private static final double UNIMOD_MASS_ERROR = 0.1;
    //TODO  - current grabs protein accessions from the defline, 
    //by space - need to implement other options
    private final String defline_regex = " ";

    private ReadUnimod unimodDoc;
    private MSSearchSettings settings;
    private MSSearch mssearch;
    private Measure mzMeasure;
    private Measure intMeasure;
    private Measure errorMeasure;

    //This is the scale to get correct MZ values out - get reset from the omx
    private int responseScale = 100;
    //Counter used to create unique ID for SpectrumIdentificationResult
    private int sirCounter = 1;

    private Map<String, DBSequence> foundProts;
    private Map<String, String> pepProtMap;

    //lookup to get a peptide by peptideseq_varmods_fixedmods
    private Map<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide> peptideLookup;

    //lookup to get a peptide evidence object by peptideID_proteinacc_start_end
    private Map<String, PeptideEvidence> pepEvidLookup;
    private Map<String, Boolean> uniquePeps;

    //Mapping from Omssa Mod integers to objects
    private Map<Integer, OmssaModification> intToModMap;
    private Map<MSSpectrum, MSHitSet> specToHitMap;
    private Map<MSSpectrum, MSHitSet> results;

    public Omssa2mzidMzidContainer(String input,
                                   Boolean outputFrags,
                                   String decoyRegex,
                                   String omssaModsFile,
                                   String omssaUserModsFile,
                                   MzIdentMLVersion ver) {
        this.inputOmssaFileName = input;
        boolean cleanupOmssaMods = false;
        if (omssaModsFile == null || omssaModsFile.equals("")) {
            System.out.println("Using the default mods.xml file.\n");
            this.modsFile = "mods.xml";
            InputStream inMods = ClassLoader.getSystemClassLoader()
                    .getResourceAsStream(modsFile);
            extractFileFromJar(inMods, modsFile);
            cleanupOmssaMods = true;
        } else {
            this.modsFile = omssaModsFile;
        }

        boolean cleanupUserMods = false;
        if (userModsFile == null || userModsFile.equals("")) {
            System.out.println("Using the default usermods.xml file\n");
            userModsFile = "usermods.xml";
            InputStream inUserMods = ClassLoader.getSystemClassLoader()
                    .getResourceAsStream("usermods.xml");
            extractFileFromJar(inUserMods, userModsFile);
            cleanupUserMods = true;
        } else {
            this.userModsFile = omssaUserModsFile;
        }

        this.omxFile = new OmssaOmxFile(input, modsFile, userModsFile);

        if (decoyRegex != null) {
            this.decoyRegularExpression = decoyRegex;
        }

        this.outputFragmentation = outputFrags;

        //determin mzid file version
        if (null == ver) {
            this.version = MzIdentMLVersion.Version_1_1;
        } else {
            this.version = ver;
        }

        parseFile(this.omxFile);
        init();

        if (cleanupOmssaMods) {
            File modFile = new File(modsFile);
            modFile.delete();
        }
        if (cleanupUserMods) {
            File userModFile = new File(userModsFile);
            userModFile.delete();
        }

    }

    public Omssa2mzidMzidContainer(String input,
                                   String decoyRegex,
                                   String omssaModsFile,
                                   String omssaUserModsFile,
                                   MzIdentMLVersion ver) {
        this(input, false, decoyRegex, omssaModsFile, omssaUserModsFile, ver);
    }

    /*
     * Helper method to get the mod file out of the jar, since local file is
     * needed for OMXParser
     */
    private void extractFileFromJar(InputStream in, String filename) {

        StringBuilder builder = new StringBuilder();
        try (BufferedReader br
                = new BufferedReader(new InputStreamReader(in,
                                                           StandardCharsets.UTF_8))) {
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                builder.append(line);
            }
            String text = builder.toString();
            try (FileWriter writer = new FileWriter(filename)) {
                writer.write(text);
            }
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1]
                    .getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + e
                    .getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information "
                    + "on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    private void parseFile(OmssaOmxFile omxFile) {
        this.results = omxFile.getSpectrumToHitSetMap();
        this.mssearch = omxFile.getParserResult();
        this.specToHitMap = omxFile
                .getSpectrumToHitSetMap();
        this.intToModMap = omxFile.getModifications();
    }

    @Override
    public AnalysisCollection getAnalysisCollection() {

        SpectrumIdentification specIdent = new SpectrumIdentification();
        specIdent.setId(SPEC_IDENT_ID);
        specIdent.setSpectrumIdentificationList(spectrumIdentList);
        specIdent.setSpectrumIdentificationProtocol(
                this.spectrumIdentificationProtocol);
        List<SearchDatabaseRef> searchDbRefList = specIdent
                .getSearchDatabaseRef();
        SearchDatabaseRef searchDbRef = new SearchDatabaseRef();
        searchDbRef.setSearchDatabase(this.searchDatabase);
        searchDbRefList.add(searchDbRef);

        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
        InputSpectra inputSpec = new InputSpectra();
        inputSpec.setSpectraData(spectraData);
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
        sipList.add(this.spectrumIdentificationProtocol);

        return analysisProtocolCollection;
    }

    @Override
    public AnalysisSampleCollection getAnalysisSampleCollection() {
        return null;
    }

    @Override
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();
        analysisSoftwareList.getAnalysisSoftware().add(analysisSoftwareOmssa);
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
        searchDbList.add(this.searchDatabase);

        SourceFile sourceFile = new SourceFile();
        File inputFile = new File(this.inputOmssaFileName);
        sourceFile.setLocation(inputFile.getAbsolutePath());
        sourceFile.setId(SOURCE_FILE_ID);
        FileFormat ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam("MS:1001400", "OMSSA xml file",
                                               CvConstants.PSI_CV));

        List<SourceFile> sourceFileList = inputs.getSourceFile();
        sourceFile.setFileFormat(ff);
        sourceFileList.add(sourceFile);

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
        return sequenceCollection;
    }

    @Override
    public SpectrumIdentificationList getSpectrumIdentificationList() {
        return this.spectrumIdentList;
    }

    private void init() {

        //analysisSoftwareXtandem
        this.analysisSoftwareOmssa = MzidLibUtils.createAnalysisSoftware(
                "OMSSA",
                ANALYSIS_SOFT_ID,
                CvConstants.OMSSA,
                "");

        //docOwner
        docOwner = createDocOwner("firstName", "secondName", "address");

        unimodDoc = new ReadUnimod();

        //searchDatabase
        this.searchDatabase = createSearchDatabase();

        //Must initialise measures before createFragmentationTable
        /*
         * <Measure id="m_mz"> <cvParam cvRef="PSI-MS" accession="MS:1001225"
         * name="product ion m/z"/> </Measure> <Measure id="m_intensity">
         * <cvParam cvRef="PSI-MS" accession="MS:1001226" name="product ion
         * intensity"/> </Measure> <Measure id="m_error"> <cvParam
         * cvRef="PSI-MS" accession="MS:1001227" name="product ion m/z error"
         * unitAccession="MS:1000040" unitName="m/z" unitCvRef="PSI-MS"/>
         * </Measure>
         *
         */
        mzMeasure = new Measure();
        mzMeasure.setId(MEASURE_MZ_ID);
        List<CvParam> cvParamList = mzMeasure.getCvParam();
        cvParamList.add(MzidLibUtils
                .makeCvParam("MS:1001225", "product ion m/z", CvConstants.PSI_CV,
                             "MS:1000040",
                             "m/z", CvConstants.UNIT_CV));
        intMeasure = new Measure();
        intMeasure.setId(MEASURE_INT_ID);
        cvParamList = intMeasure.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001226",
                                                 "product ion intensity",
                                                 CvConstants.PSI_CV,
                                                 "MS:1000131",
                                                 "number of counts",
                                                 CvConstants.UNIT_CV));
        errorMeasure = new Measure();
        errorMeasure.setId(MEASURE_ERROR_ID);
        cvParamList = errorMeasure.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001227",
                                                 "product ion m/z error",
                                                 CvConstants.PSI_CV,
                                                 "MS:1000040", "m/z",
                                                 CvConstants.UNIT_CV));

        //fragmentationTable
        this.fragmentationTable = createFragmentationTable();

        //spectraData
        this.spectraData = createSpectraData();

        //spectrumIdentificationProtocol
        this.spectrumIdentificationProtocol
                = createSpectrumIdentificationProtocol();

        //spectrumIdList and sequenceCollection
        createSpectrumIdentificationListAndSequenceCollection();

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

    private SearchDatabase createSearchDatabase() {
        MSSearch_response responses = mssearch.MSSearch_response;
        MSSearch_request requests = mssearch.MSSearch_request;

        int reqCounter = 0;

        for (MSRequest request : requests.MSRequest) {

            MSRequest_settings reqSettings = request.MSRequest_settings;
            settings = reqSettings.MSSearchSettings;
            reqCounter++;

            if (reqCounter > 1) {
                System.out.println(
                        "Error: multiple requests in the OMX file - "
                        + "this is not currently supported");
            }
        }

        MSResponse response = new MSResponse();
        int respCounter = 0;
        for (MSResponse resp : responses.MSResponse) {
            response = resp;
            responseScale = response.MSResponse_scale;

            respCounter++;

            if (respCounter > 1) {
                System.out.println(
                        "Error: multiple responses in the OMX file - "
                        + "this is not currently supported");
            }

            // String version = response.MSResponse_version;
            //TO DO get other running params etc.
            //String bioseq = response.MSResponse_bioseqs;
            //System.out.println("Version: " + version + " bioseqs: " + bioseq);
        }

        SearchDatabase searchDb = new SearchDatabase();
        searchDb.setId(SEARCH_DB_ID);
        searchDb.setNumDatabaseSequences((long) response.MSResponse_dbversion);

        UserParam param = new UserParam();
        param.setName(settings.MSSearchSettings_db);
        Param tempParam = new Param();
        tempParam.setParam(param);
        searchDb.setDatabaseName(tempParam);
        searchDb.setLocation(settings.MSSearchSettings_db);

        FileFormat ff = new FileFormat();
        ff.setCvParam(CvConstants.FASTA_FORMAT);
        //TODO - this should not be hard coded 
        //<cvParam accession="MS:1001348" name="FASTA format" cvRef="PSI-MS"/>
        searchDb.setFileFormat(ff);
        return searchDb;
    }

    private FragmentationTable createFragmentationTable() {

        FragmentationTable fragTable = new FragmentationTable();
        List<Measure> measureList = fragTable.getMeasure();
        measureList.add(mzMeasure);
        measureList.add(intMeasure);
        measureList.add(errorMeasure);

        return fragTable;
    }

    private String convertVarMods(Peptide mzidPep, MSHits_mods mods) {

        String modString = "";
        List<Modification> allMods = mzidPep.getModification();
        String pepSeq = mzidPep.getPeptideSequence();
        for (MSModHit mod : mods.MSModHit) {
            MSModHit_modtype modType = mod.MSModHit_modtype;
            int modNum = modType.MSMod;
            int modSite = mod.MSModHit_site;
            modString += "_" + modNum + "@" + modSite;
            //System.out.println("\t" + modString);
            Modification mzidmod = new Modification();

            CvParam modParam = new CvParam();
            OmssaModification omod = intToModMap.get(modNum);
            double monoMass = omod.getModMonoMass();

            mzidmod.setMonoisotopicMassDelta(monoMass);
            //+1 since Omssa counts from zero, mzid counts first position in peptide as 1
            mzidmod.setLocation(modSite + 1);

            boolean isMono = true;
            ModT unimod = unimodDoc.getModByMass(monoMass, UNIMOD_MASS_ERROR,
                                                 isMono, pepSeq.charAt(modSite));

            if (unimod != null) {
                mzidmod.getResidues().add("" + pepSeq.charAt(modSite));
            }

            if (unimod == null && modSite == 0) {
                //See if this is a possible N-terminal mod
                System.out.println(
                        "\tNot found, so look to see if it is N-terminal\n");
                unimod = unimodDoc.getModByMass(monoMass, UNIMOD_MASS_ERROR,
                                                isMono, '[');
                mzidmod.setLocation(0);
            }
            if (unimod == null && modSite == pepSeq.length()) {
                //See if this is a possible C-terminal mod
                System.out.println(
                        "\tNot found, so look to see if it is C-terminal\n");
                unimod = unimodDoc.getModByMass(monoMass, UNIMOD_MASS_ERROR,
                                                isMono, ']');
                mzidmod.setLocation(pepSeq.length() + 1);
            }

            if (unimod != null) {
                modParam.setAccession("UNIMOD:" + unimod.getRecordId());
                modParam.setCv(CvConstants.UNIMOD_CV);
                modParam.setName(unimod.getTitle());

            } else {
                System.out.println(
                        "Error: modification with mass not recognized");
                modParam.setName("unknown modification");
                modParam.setCv(CvConstants.PSI_CV);
                modParam.setAccession("MS:1001460");
            }

            List<CvParam> paramList = mzidmod.getCvParam();
            paramList.add(modParam);
            allMods.add(mzidmod);
        }
        return modString;
    }

    private void convertFixedMods(Peptide mzidPep) {

        // get the list of fixed modifications
        List<Integer> fixedModifications
                = settings.MSSearchSettings_fixed.MSMod;

        String pepSeq = mzidPep.getPeptideSequence();

        List<Modification> allMods = mzidPep.getModification();

        // fixed modifications
        if (fixedModifications.size() > 0) {

            for (Integer fixedModification : fixedModifications) {

                OmssaModification omod = intToModMap.get(fixedModification);
                List<String> modifiedResidues = omod.getModResidues();

                for (String modifiedResidue : modifiedResidues) {
                    int index = pepSeq.indexOf(modifiedResidue);

                    while (index != -1) {
                        Modification mzidmod = new Modification();
                        mzidmod.getResidues().add(modifiedResidue);
                        CvParam modParam = new CvParam();
                        double monoMass = omod.getModMonoMass();
                        mzidmod.setMonoisotopicMassDelta(monoMass);

                        //If second res is modified, index would return 1, 
                        //but mzid position should be 2
                        mzidmod.setLocation(index + 1);
                        boolean isMono = true;
                        ModT unimod = unimodDoc.getModByMass(monoMass,
                                                             UNIMOD_MASS_ERROR,
                                                             isMono, pepSeq
                                                             .charAt(index));
                        if (unimod != null) {
                            modParam.setAccession("UNIMOD:" + unimod
                                    .getRecordId());
                            modParam.setCv(CvConstants.UNIMOD_CV);
                            modParam.setName(unimod.getTitle());
                        } else {
                            System.out.println(
                                    "Error: modification with mass not recognized");
                            modParam.setName("unknown modification");
                            modParam.setCv(CvConstants.PSI_CV);
                            modParam.setAccession("MS:1001460");
                        }
                        List<CvParam> paramList = mzidmod.getCvParam();
                        paramList.add(modParam);
                        allMods.add(mzidmod);
                        index = pepSeq.indexOf(modifiedResidue, index + 1);
                    }
                }

                //Candidate N or C terminal mods
                if (modifiedResidues == null || modifiedResidues.isEmpty()) {

                    boolean isPepNTerminalMod = false;
                    boolean isPepCTerminalMod = false;
                    boolean isProtNTerminalMod = false;
                    boolean isProtCTerminalMod = false;

                    int modType = omod.getModType();
                    switch (modType) {
                        case 1:
                        case 2:
                            isProtNTerminalMod = true;
                            break;
                        case 5:
                        case 6:
                            isPepNTerminalMod = true;
                            break;
                        case 3:
                        case 4:
                            isProtCTerminalMod = true;
                            break;
                        case 7:
                        case 8:
                            isPepCTerminalMod = true;
                            break;
                        default:
                            break;
                    }
                    boolean isMono = true;
                    ModT unimod = null;

                    if (isPepNTerminalMod || isPepCTerminalMod) {

                        Modification mzidmod = new Modification();

                        CvParam modParam = new CvParam();
                        double monoMass = intToModMap.get(fixedModification)
                                .getModMonoMass();
                        mzidmod.setMonoisotopicMassDelta(monoMass);

                        if (isPepNTerminalMod) {
                            unimod = unimodDoc.getModByMass(monoMass,
                                                            UNIMOD_MASS_ERROR,
                                                            isMono, '[');
                            mzidmod.setLocation(0);
                        } else if (isPepCTerminalMod) {
                            unimod = unimodDoc.getModByMass(monoMass,
                                                            UNIMOD_MASS_ERROR,
                                                            isMono, ']');
                            mzidmod.setLocation(pepSeq.length() + 1);
                        }

                        if (unimod != null) {
                            modParam.setAccession("UNIMOD:" + unimod
                                    .getRecordId());
                            modParam.setCv(CvConstants.UNIMOD_CV);
                            modParam.setName(unimod.getTitle());
                        } else {
                            System.out.println(
                                    "Error: modification with mass not recognized");
                            modParam.setName("unknown modification");
                            modParam.setCv(CvConstants.PSI_CV);
                            modParam.setAccession("MS:1001460");
                        }

                        List<CvParam> paramList = mzidmod.getCvParam();
                        paramList.add(modParam);
                        allMods.add(mzidmod);
                    } else {
                        //TODO Fixed protein N or C terminal mods not yet supported
                        /*
                         * else if(isProtNTerminalMod){
                         * paramList.add(makeCvParam("MS:1002057","modification
                         * specificity protein N-term",psiCV)); unimod =
                         * unimodDoc.getModByMass(monoMass, unimodMassError,
                         * isMono, '['); mzidmod.setLocation(0); } else
                         * if(isProtCTerminalMod){
                         * paramList.add(makeCvParam("MS:1002058","modification
                         * specificity protein C-term",psiCV)); unimod =
                         * unimodDoc.getModByMass(monoMass, unimodMassError,
                         * isMono, ']'); mzidmod.setLocation(pepSeq.length()+1);
                         * }
                         */
                    }
                }
            }
        }
    }

    private double[] getMatchedIon(MSSpectrum spectrum, double ionMz) {

        double[] matchedPeak = new double[2];
        MSSpectrum_mz specMz = spectrum.MSSpectrum_mz;
        List<Integer> mzVals = specMz.MSSpectrum_mz_E;
        MSSpectrum_abundance specInt = spectrum.MSSpectrum_abundance;
        List<Integer> abundValues = specInt.MSSpectrum_abundance_E;

        double error = 100;

        int count = 0;
        int foundPos = -1;

        double foundMz = -1;
        for (int mzVal : mzVals) {
            double mz = (double) mzVal / responseScale;

            double tempError = mz - ionMz;

            if (Math.abs(tempError) < Math.abs(error)) {
                error = tempError;
                foundPos = count;
                foundMz = mz;
            }

            count++;
            //System.out.print(mz + "\t");
        }

        //System.out.println("");
        double foundInt = 0.0;
        if (foundPos != -1) {
            /*
             * for(Integer intVal: abundValues){ double intensity = (double)
             * intVal/spectrum.MSSpectrum_iscale; //System.out.print(intensity +
             * "\t"); }
             */

            //System.out.println("\n");
            double intensityScale = spectrum.MSSpectrum_iscale;
            foundInt = (double) (abundValues.get(foundPos));
            foundInt = foundInt / intensityScale;

        }

        matchedPeak[0] = foundMz;
        matchedPeak[1] = foundInt;
        return matchedPeak;
    }

    private SearchModification handleOmssaMod(OmssaModification omod,
                                              boolean isFixed) {
        SearchModification searchMod = new SearchModification();
        double monoMass = omod.getModMonoMass();
        List<String> residues = omod.getModResidues();

        boolean isPepNTerminalMod = false;
        boolean isPepCTerminalMod = false;
        boolean isProtNTerminalMod = false;
        boolean isProtCTerminalMod = false;

        int modType = omod.getModType();

        switch (modType) {
            case 1:
            case 2:
                isProtNTerminalMod = true;
                break;
            case 5:
            case 6:
                isPepNTerminalMod = true;
                break;
            case 3:
            case 4:
                isProtCTerminalMod = true;
                break;
            case 7:
            case 8:
                isPepCTerminalMod = true;
                break;
            default:
                break;
        }
        boolean isMono = true;

        if (isPepNTerminalMod || isProtNTerminalMod) {
            residues.add("[");
            residues.remove(residues.size() - 1); //Remove the ] or [ char
        } else if (isPepCTerminalMod || isProtCTerminalMod) {
            residues.add("]");
            residues.remove(residues.size() - 1); //Remove the ] or [ char
        }
        ModT unimod = unimodDoc
                .getModByMass(monoMass, UNIMOD_MASS_ERROR, isMono,
                              residues);

        searchMod.setFixedMod(isFixed);

        List<CvParam> modCvParamList = searchMod.getCvParam();
        CvParam modParam = new CvParam();

        if (unimod != null) {
            modParam = MzidLibUtils
                    .makeCvParam("UNIMOD:" + unimod.getRecordId(), unimod.
                                 getTitle(), CvConstants.UNIMOD_CV);

        } else {
            modParam.setName("unknown modification");
            modParam.setCv(CvConstants.PSI_CV);
            modParam.setAccession("MS:1001460");
        }

        modCvParamList.add(modParam);
        searchMod.setMassDelta(new Float(monoMass));
        List<String> residueList = searchMod.getResidues();

        if (isPepNTerminalMod) {
            modCvParamList.add(MzidLibUtils.makeCvParam("MS:1001189",
                                                        "modification specificity peptide N-term",
                                                        CvConstants.PSI_CV));
        } else if (isPepCTerminalMod) {
            modCvParamList.add(MzidLibUtils.makeCvParam("MS:1001190",
                                                        "modification specificity peptide C-term",
                                                        CvConstants.PSI_CV));
        } else if (isProtNTerminalMod) {
            modCvParamList.add(MzidLibUtils.makeCvParam("MS:1002057",
                                                        "modification specificity protein N-term",
                                                        CvConstants.PSI_CV));
        } else if (isProtCTerminalMod) {
            modCvParamList.add(MzidLibUtils.makeCvParam("MS:1002058",
                                                        "modification specificity protein C-term",
                                                        CvConstants.PSI_CV));
        }

        if (!residues.isEmpty()) {
            for (String residue : residues) {
                residueList.add(residue);
            }
        } else {
            residueList.add(".");               //any chars character
        }

        return searchMod;

    }

    private Enzyme getEnzyme(String omssaEnzyme, int missedCleavage) {

        Enzyme enzyme = new Enzyme();
        //[KR]|{P}

        //TODO only a few enzymes implemented
        enzyme.setId("Enz1");
        enzyme.setCTermGain("OH");
        enzyme.setNTermGain("H");
        enzyme.setMissedCleavages(missedCleavage);
        enzyme.setSemiSpecific(false);
        ParamList paramList = enzyme.getEnzymeName();
        if (paramList == null) {
            paramList = new ParamList();
            enzyme.setEnzymeName(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();

        if (omssaEnzyme.equalsIgnoreCase("trypsin")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001251", "Trypsin",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Arg-C")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001303", "Arg-C",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("CNBr")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001307", "CNBr",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Chymotrypsin")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001306",
                                                     "Chymotrypsin",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Formic Acid")) {
            cvParamList.add(MzidLibUtils.makeCvParam("TODO", "TODO",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Lys-C")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001309", "Lys-C",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Lys-C, no P rule")) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001310", "Lys-C/P",
                                                     CvConstants.PSI_CV));
        } else if (omssaEnzyme.equalsIgnoreCase("Pepsin A")) {
            cvParamList.add(MzidLibUtils.makeCvParam("TODO", "TODO",
                                                     CvConstants.PSI_CV));
        } else {
            cvParamList.add(MzidLibUtils.makeCvParam("TODO", "TODO",
                                                     CvConstants.PSI_CV));
            /*
             * 0: Trypsin 1: Arg-C 2: CNBr 3: Chymotrypsin 4: Formic Acid 5:
             * Lys-C 6: Lys-C, no P rule 7: Pepsin A 8: Trypsin+CNBr 9:
             * Trypsin+Chymotrypsin 10: Trypsin, no P rule 11: Whole protein 12:
             * Asp-N 13: Glu-C 14: Asp-N+Glu-C 15: Top-Down 16: Semi-Tryptic 17:
             * No Enzyme 18: Chymotrypsin, no P rule 19: Asp-N (DE) 20: Glu-C
             * (DE)
             */

            //TODO
            /*
             *
             * [Term] id: MS:1001303 name: Arg-C is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001272 ! (?<=R)(?!P)
             *
             * [Term] id: MS:1001304 name: Asp-N is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001273 ! (?=[BD])
             *
             * [Term] id: MS:1001305 name: Asp-N_ambic is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001274 !
             * (?=[DE])
             *
             * [Term] id: MS:1001306 name: Chymotrypsin is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001332 !
             * (?<=[FYWL])(?!P)
             *
             * [Term] id: MS:1001307 name: CNBr is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001333 ! (?<=M)
             *
             * [Term] id: MS:1001308 name: Formic_acid is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001334 !
             * ((?<=D))|((?=D))
             *
             * [Term] id: MS:1001309 name: Lys-C is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001335 ! (?<=K)(?!P)
             *
             * [Term] id: MS:1001310 name: Lys-C/P is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001336 ! (?<=K)
             *
             * [Term] id: MS:1001311 name: PepsinA is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001337 ! (?<=[FL])
             *
             * [Term] id: MS:1001312 name: TrypChymo is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001338 !
             * (?<=[FYWLKR])(?!P)
             *
             * [Term] id: MS:1001313 name: Trypsin/P is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001339 ! (?<=[KR])
             *
             * [Term] id: MS:1001314 name: V8-DE is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001340 ! (?<=[BDEZ])(?!P)
             *
             * [Term] id: MS:1001315 name: V8-E is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001341 ! (?<=[EZ])(?!P)
             */
        }

        return enzyme;

    }

    private CvParam getFileFormatCVParam(int omssaFileType) {

        /*
         * MSSpectrumFileType.put(new Integer(index++), "dta");
         * MSSpectrumFileType.put(new Integer(index++), "dtablank");
         * MSSpectrumFileType.put(new Integer(index++), "dtaxml");
         * MSSpectrumFileType.put(new Integer(index++), "asc");
         * MSSpectrumFileType.put(new Integer(index++), "pkl");
         * MSSpectrumFileType.put(new Integer(index++), "pks");
         * MSSpectrumFileType.put(new Integer(index++), "sciex");
         * MSSpectrumFileType.put(new Integer(index++), "mgf");
         * MSSpectrumFileType.put(new Integer(index++), "unknown");
         * MSSpectrumFileType.put(new Integer(index++), "oms - asn.1 binary for
         * iterative search"); MSSpectrumFileType.put(new Integer(index++), "omx
         * - xml for iterative search"); MSSpectrumFileType.put(new
         * Integer(index++), "xml - xml MSRequest"); MSSpectrumFileType.put(new
         * Integer(index++), "omxbz - bzip2 omx file2");
         *
         * "MS:1001062","Mascot MGF file"
         */
        CvParam cvParam = new CvParam();
        cvParam.setCv(CvConstants.PSI_CV);

        switch (omssaFileType) {
            case 0:
                cvParam.setAccession("MS:1000613");
                cvParam.setName("DTA file");
                break;
            case 1:
                cvParam.setAccession("MS:00000NOID");
                cvParam.setName("dtablank");
                break;
            case 2:
                cvParam.setAccession("MS:00000NOID");
                cvParam.setName("dtaxml");
                break;
            case 3:
                cvParam.setAccession("MS:00000NOID");
                cvParam.setName("asc");
                break;
            case 4:
                cvParam.setAccession("MS:1000565");
                cvParam.setName("Micromass PKL file");
                break;
            case 5:
                cvParam.setAccession("MS:1001245");
                cvParam.setName("PerSeptive PKS file");
                break;
            case 6:
                cvParam.setAccession("MS:1001246");
                cvParam.setName("Sciex API III file");
                break;
            case 7:
                cvParam.setAccession("MS:1001062");
                cvParam.setName("Mascot MGF file");
                break;
            case 8:
                cvParam.setAccession("MS:00000NOID");
                cvParam.setName("oms - asn.1 binary for iterative search");
                break;
            case 9:
                cvParam.setAccession(" MS:1001400");
                cvParam.setName("OMSSA xml file");
                break;
            case 10:
                cvParam.setAccession("MS:00000NOID");
                cvParam.setName("oms - asn.1 binary for iterative search");
                break;
            default:
                cvParam.setAccession("MS:Unknown");
                cvParam.setName("unrecognized file format");

        }

        return cvParam;
    }

    private SpectraData createSpectraData() {
        SpectraData specData = new SpectraData();
        SpectrumIDFormat sif = new SpectrumIDFormat();
        sif.setCvParam(MzidLibUtils.makeCvParam("MS:1000774",
                                                "multiple peak list nativeID format",
                                                CvConstants.PSI_CV));
        specData.setSpectrumIDFormat(sif);

        int fileType
                = settings.MSSearchSettings_infiles.MSInFile.MSInFile_infiletype.MSSpectrumFileType;
        FileFormat ff = new FileFormat();
        ff.setCvParam(getFileFormatCVParam(fileType));
        specData.setFileFormat(ff);

        specData.setId(SPECTRA_DATA_ID);
        specData.setLocation(
                settings.MSSearchSettings_infiles.MSInFile.MSInFile_infile);

        return specData;
    }

    private SpectrumIdentificationProtocol createSpectrumIdentificationProtocol() {
        SpectrumIdentificationProtocol specIdProtocol
                = new SpectrumIdentificationProtocol();
        specIdProtocol.setId(SI_PROTOCOL_ID);
        specIdProtocol.setAnalysisSoftware(analysisSoftwareOmssa);

        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
        Param tempParam = new Param();
        tempParam.setParam(CvConstants.MS_MS_SEARCH);
        specIdProtocol.setSearchType(tempParam);

        ParamList paramList = specIdProtocol.getAdditionalSearchParams();
        if (paramList == null) {
            paramList = new ParamList();
            specIdProtocol.setAdditionalSearchParams(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();

        //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact
        int msSearchType
                = settings.MSSearchSettings_precursorsearchtype.MSSearchType;
        switch (msSearchType) {
            case 0:
                cvParamList.add(CvConstants.PARENT_MASS_TYPE_MONO);
                break;
            case 1:
                cvParamList.add(CvConstants.PARENT_MASS_TYPE_AVERAGE);
                break;
            case 2:
                cvParamList.add(MzidLibUtils.makeCvParam("MS:No_acc",
                                                         "monoisotopic N15",
                                                         CvConstants.PSI_CV));
                System.out.println(
                        "Warning: No CV term for monoisotopic N15 search type");
                break;
            case 3:
                cvParamList.add(MzidLibUtils
                        .makeCvParam("MS:No_acc", "exact", CvConstants.PSI_CV));
                System.out.println(
                        "Warning: No CV term for exact mass search type");
                break;
            default:
                System.out.println("Error search type not recognised");
                break;
        }

        //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact
        int prodSearchType
                = settings.MSSearchSettings_productsearchtype.MSSearchType;

        switch (prodSearchType) {
            case 0:
                cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_MONO);
                break;
            case 1:
                cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_AVERAGE);
                break;
            case 2:
                cvParamList.add(MzidLibUtils.makeCvParam("MS:No_acc",
                                                         "monoisotopic N15",
                                                         CvConstants.PSI_CV));
                System.out.println(
                        "Warning: No CV term for monoisotopic N15 search type");
                break;
            case 3:
                cvParamList.add(MzidLibUtils
                        .makeCvParam("MS:No_acc", "exact", CvConstants.PSI_CV));
                System.out.println(
                        "Warning: No CV term for exact mass search type");
                break;
            default:
                System.out.println("Error search type not recognised");
                break;
        }

        // Add "no special processing" cv term if this is mzid 1.2 version
        if (version.equals(MzIdentMLVersion.Version_1_2)) {
            cvParamList.add(CvConstants.NO_SPECIAL_PROCESSING);
        }

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams
                .getSearchModification();

        settings.MSSearchSettings_fixed.MSMod.stream()
                .map((fixedMod) -> intToModMap.get(fixedMod))
                .map((omod) -> handleOmssaMod(omod, true))
                .forEach((searchMod) -> {
                    searchModList.add(searchMod);
                });

        settings.MSSearchSettings_variable.MSMod.stream()
                .map((varMod) -> intToModMap.get(varMod))
                .map((omod) -> handleOmssaMod(omod, false))
                .forEach((searchMod) -> {
                    searchModList.add(searchMod);
                });

        specIdProtocol.setModificationParams(modParams);

        Enzymes enzymes = specIdProtocol.getEnzymes();
        if (enzymes == null) {
            enzymes = new Enzymes();
            specIdProtocol.setEnzymes(enzymes);
        }
        enzymes.setIndependent(false);
        List<Enzyme> enzymeList = enzymes.getEnzyme();
        List<Integer> msEnzymeList = settings.MSSearchSettings_enzyme.MSEnzymes;
        for (Integer msEnzyme : msEnzymeList) {
            //OmssaEnumerators omssaEnums = new OmssaEnumerators();
            Enzyme enzyme
                    = getEnzyme(OmssaEnumerators.getEnzymeAsText(msEnzyme),
                                settings.MSSearchSettings_missedcleave);
            enzymeList.add(enzyme);
        }

        CvParam fragCvPlus = MzidLibUtils.getCvParamWithMassUnits(true);
        CvParam fragCvMinus = MzidLibUtils.getCvParamWithMassUnits(true);

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
        fragCvPlus.setAccession("MS:1001412");
        fragCvPlus.setName("search tolerance plus value");
        fragCvMinus.setAccession("MS:1001413");
        fragCvMinus.setName("search tolerance minus value");
        fragCvPlus.setValue("" + settings.MSSearchSettings_msmstol);
        fragCvMinus.setValue("" + settings.MSSearchSettings_msmstol);
        Tolerance fragTol = new Tolerance();
        List<CvParam> fragCvList = fragTol.getCvParam();
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        CvParam parCvPlus = MzidLibUtils.getCvParamWithMassUnits(true);
        CvParam parCvMinus = MzidLibUtils.getCvParamWithMassUnits(true);

        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue("" + settings.MSSearchSettings_peptol);
        parCvMinus.setValue("" + settings.MSSearchSettings_peptol);
        Tolerance parTol = new Tolerance();
        List<CvParam> parCvList = parTol.getCvParam();
        parCvList.add(parCvPlus);
        parCvList.add(parCvMinus);

        specIdProtocol.setFragmentTolerance(fragTol);
        specIdProtocol.setParentTolerance(parTol);

        ParamList sipParamList = specIdProtocol.getThreshold();
        if (sipParamList == null) {
            sipParamList = new ParamList();
            specIdProtocol.setThreshold(sipParamList);
        }
        cvParamList = sipParamList.getCvParam();

        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001494", "no threshold",
                                                 CvConstants.PSI_CV));
        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />

        return specIdProtocol;

    }

    private void createSpectrumIdentificationListAndSequenceCollection() {
        double maxEvalue = 0.0;
        foundProts = new HashMap<>();
        pepProtMap = new HashMap<>();

        //lookup to get a peptide by peptideseq_varmods_fixedmods
        peptideLookup = new HashMap<>();
        pepEvidLookup = new HashMap<>();
        uniquePeps = new HashMap<>();

        spectrumIdentList = new SpectrumIdentificationList();
        spectrumIdentList.setId(SI_LIST_ID);

        spectrumIdentList.setFragmentationTable(this.fragmentationTable);

        sequenceCollection = new SequenceCollection();
        List<SpectrumIdentificationResult> specIdentResults = spectrumIdentList
                .getSpectrumIdentificationResult();
        List<PeptideEvidence> peptideEvidenceList = sequenceCollection
                .getPeptideEvidence();
        List<DBSequence> dbSeqList = sequenceCollection.getDBSequence();
        List<Peptide> peptideList = sequenceCollection.getPeptide();

        // Iterate over all the spectra
        Iterator<MSSpectrum> iter = results.keySet().iterator();
        while (iter.hasNext()) {

            // Get the next spectrum.
            MSSpectrum spectrum = iter.next();

            MSHitSet msHitSet = specToHitMap.get(spectrum);
            SpectrumIdentificationResult specIdentRes = null;
            List<SpectrumIdentificationItem> siiList = null;

            MSHitSet_hits msHits = msHitSet.MSHitSet_hits;

            /**
             * ***********************************************
             *  *** Setup SpectrumIdentificationResult ****
             * *********************************************
             */
            if (!msHits.MSHits.isEmpty()) {
                specIdentRes = new SpectrumIdentificationResult();
                siiList = specIdentRes.getSpectrumIdentificationItem();
                specIdentResults.add(specIdentRes);

                specIdentRes.setSpectraData(spectraData);
                specIdentRes
                        .setSpectrumID("index=" + spectrum.MSSpectrum_number);

                MSSpectrum_ids specIDs = spectrum.MSSpectrum_ids;

                if (!specIDs.MSSpectrum_ids_E.isEmpty()) {
                    List<CvParam> sirCvParamList = specIdentRes.getCvParam();
                    for (String id : specIDs.MSSpectrum_ids_E) {
                        CvParam cvp = MzidLibUtils.makeCvParam("MS:1000796",
                                                               "spectrum title",
                                                               CvConstants.PSI_CV,
                                                               id);
                        sirCvParamList.add(cvp);
                    }
                }

                int rank = 1;

                //Counter used to create unique ID for SpectrumIdentificationItem
                int siiCounter = 1;
                int pepCounter = 0;
                int pepEvidCounter = 0;
                for (MSHits hits : msHits.MSHits) {

                    SpectrumIdentificationItem sii
                            = new SpectrumIdentificationItem();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = sii
                            .getPeptideEvidenceRef();
                    siiList.add(sii);

                    int charge = hits.MSHits_charge;
                    String pre = hits.MSHits_pepstart;
                    String post = hits.MSHits_pepstop;

                    sii.setId("SII_" + sirCounter + "_" + siiCounter);
                    siiCounter++;

                    sii.setChargeState(charge);
                    sii.setRank(rank);
                    //by default these are supposed to be set to true
                    sii.setPassThreshold(true);
                    rank++;

                    long expMz = java.lang.Math.round(hits.MSHits_mass
                            / hits.MSHits_charge);
                    long theoMass = java.lang.Math.round(hits.MSHits_theomass
                            / hits.MSHits_charge);

                    double evalue = hits.MSHits_evalue;

                    if (evalue > maxEvalue) {
                        maxEvalue = evalue;
                    }

                    sii.setExperimentalMassToCharge((double) expMz / 1000);
                    sii.setCalculatedMassToCharge((double) theoMass / 1000);

                    List<CvParam> cvParamList = sii.getCvParam();

                    cvParamList.add(MzidLibUtils.makeCvParam("MS:1001328",
                                                             "OMSSA:evalue",
                                                             CvConstants.PSI_CV,
                                                             "" + evalue));

                    double pvalue = hits.MSHits_pvalue;
                    cvParamList.add(MzidLibUtils.makeCvParam("MS:1001329",
                                                             "OMSSA:pvalue",
                                                             CvConstants.PSI_CV,
                                                             "" + pvalue));

                    MSHits_mods mods = hits.MSHits_mods;

                    List<IonType> ionTypeList = null;
                    if (outputFragmentation) {
                        MSHits_mzhits mzhits = hits.MSHits_mzhits;
                        Map<String, List<MSMZHit>> mapIonsToHits
                                = new HashMap<>();

                        Fragmentation frag = new Fragmentation();
                        sii.setFragmentation(frag);
                        ionTypeList = frag.getIonType();

                        for (MSMZHit mzhit : mzhits.MSMZHit) {
                            MSMZHit_ion ion = mzhit.MSMZHit_ion;
                            int ionType = ion.MSIonType;

                            //Logic is slightly different from Tandem parser, 
                            //since all ion types appear to be mixed up in one structure
                            //1. Work out all the ion types we have, put in a temporary HashMap?
                            String testString = ionType + "_" + charge;

                            List<MSMZHit> mzHitList = null;
                            if (mapIonsToHits.containsKey(testString)) {
                                mzHitList = mapIonsToHits.get(testString);
                            } else {
                                mzHitList = new ArrayList<>();
                                mapIonsToHits.put(testString, mzHitList);
                            }
                            mzHitList.add(mzhit);

                        }

                        for (String ionKey : mapIonsToHits.keySet()) {

                            IonType mzidIon = new IonType();
                            List<Integer> ionIndexList = mzidIon.getIndex();

                            FragmentArray mzArray = new FragmentArray();
                            FragmentArray intArray = new FragmentArray();
                            FragmentArray errorArray = new FragmentArray();
                            mzArray.setMeasure(mzMeasure);
                            intArray.setMeasure(intMeasure);
                            errorArray.setMeasure(errorMeasure);

                            List<Float> mzValues = mzArray.getValues();
                            List<Float> intValues = intArray.getValues();
                            List<Float> errorValues = errorArray.getValues();

                            int j = 0;
                            for (MSMZHit mzhit : mapIonsToHits.get(ionKey)) {

                                MSMZHit_ion ion = mzhit.MSMZHit_ion;

                                int ionType = ion.MSIonType;

                                if (j == 0) {
                                    int ionCharge = mzhit.MSMZHit_charge;
                                    mzidIon.setCharge(ionCharge);
                                    //System.out.println("Lookup:" + ionType);

                                    CvParam cvParam
                                            = MzidLibUtils.getFragmentCvParam(
                                                    ionType);
                                    mzidIon.setCvParam(cvParam);
                                }
                                j++;

                                double theoMz = (double) (mzhit.MSMZHit_mz)
                                        / responseScale;
                                double[] matchedPeak = getMatchedIon(spectrum,
                                                                     theoMz);
                                mzValues.add((float) matchedPeak[0]);
                                intValues.add((float) matchedPeak[1]);
                                errorValues.add(
                                        (float) (matchedPeak[0] - theoMz));
                                ionIndexList.add(mzhit.MSMZHit_number);  //index position

                            }
                            List<FragmentArray> fragmentList = mzidIon
                                    .getFragmentArray();
                            fragmentList.add(mzArray);
                            fragmentList.add(intArray);
                            fragmentList.add(errorArray);

                            ionTypeList.add(mzidIon);
                        }
                    }
                    String pepSeq = hits.MSHits_pepstring;
                    Peptide mzidPep = new Peptide();
                    mzidPep.setPeptideSequence(pepSeq);

                    String modString = convertVarMods(mzidPep, mods);
                    //Uses the fixed mods in the search settings to guess the fixed modifications
                    convertFixedMods(mzidPep);
                    String uniquePep = pepSeq + modString;
                    if (!uniquePeps.containsKey(uniquePep)) {
                        peptideList.add(mzidPep);
                    }
                    uniquePeps.put(uniquePep, true);

                    /*
                     *********************************************************
                     ******** These are peptide evidences ****************
                     * ********************************************************
                     */
                    MSHits_pephits pephits = hits.MSHits_pephits;
                    for (MSPepHit pephit : pephits.MSPepHit) {
                        //start and end positions are mapped from zero in Omssa 
                        //but from one in the rest of the known universe
                        int start = pephit.MSPepHit_start + 1;

                        //Note this is not an issue in the CSV version
                        int end = pephit.MSPepHit_stop + 1;

                        String defline = pephit.MSPepHit_defline;
                        //String protObjID = pephit.MSPepHit_accession;
                        // Added By FG
                        String protAcc;
                        if (defline.contains(defline_regex)) {
                            protAcc = defline.substring(0, defline.indexOf(
                                                        defline_regex));
                        } else {
                            protAcc = defline;
                        }

                        //String protSeq = "";
                        int protLength = pephit.MSPepHit_protlength;

                        //TODO unclear how to retrieve protein sequences
                        //Use Hash map to test if Protein sequence has been added to DBSeq before
                        DBSequence dbSeq = null;
                        if (!foundProts.containsKey(protAcc)) {
                            dbSeq = new DBSequence();
                            foundProts.put(protAcc, dbSeq);
                            dbSeq.setAccession(protAcc);
                            //dbSeq.setSeq(protSeq);
                            dbSeq.setLength(protLength);
                            dbSeq.setId("dbseq_" + protAcc);
                            dbSeq.setSearchDatabase(this.searchDatabase);
                            List<CvParam> dbCvParamList = dbSeq.getCvParam();
                            String desc = defline;

                            dbCvParamList.add(MzidLibUtils.makeCvParam(
                                    "MS:1001088", "protein description",
                                    CvConstants.PSI_CV,
                                    desc));
                            dbSeqList.add(dbSeq);
                        } else {
                            dbSeq = foundProts.get(protAcc);
                        }

                        String testPepMods = pepSeq + "_" + modString + "_"
                                + start + "_" + end;
                        String testProt = protAcc + "_" + start + "_" + end;

                        //Check this is a unique peptide
                        //TODO - Wasteful of memory, since mzidPep doesn't need 
                        //to be created if this is not a unique peptide
                        boolean newPepEvid = false;
                        String pepId = null;

                        if (!pepProtMap.containsKey(testPepMods)) {

                            //Add peptide sequence to mzid Peptide object
                            pepId = "Peptide" + pepCounter;
                            pepCounter++;
                            mzidPep.setId(uniquePep);

                            if (!uniquePeps.containsKey(uniquePep)) {
                                peptideList.add(mzidPep);
                            }
                            uniquePeps.put(uniquePep, true);

                            peptideLookup.put(testPepMods, mzidPep);

                            pepProtMap.put(testPepMods, testProt);
                            newPepEvid = true;
                        } else {
                            String mappedProts = pepProtMap.get(testPepMods);
                            //test if testProt is within mapped prots - 
                            //if yes, do nothing, if not, create a new peptideevidence
                            if (!mappedProts.contains(testProt)) {
                                pepProtMap.put(testPepMods, mappedProts + ";"
                                               + testProt);
                                newPepEvid = true;
                            }
                            mzidPep = peptideLookup.get(testPepMods);
                        }

                        PeptideEvidence pepEvid = null;
                        if (newPepEvid) {
                            pepEvid = new PeptideEvidence();
                            pepEvidLookup.put(mzidPep.getId() + "_" + testProt,
                                              pepEvid);
                            pepEvid.setEnd(end);
                            pepEvid.setStart(start);
                            pepEvid.setPeptide(mzidPep);

                            //pepEvid.setMissedCleavages(0);  //TODO not sure we can get this
                            pepEvid.setDBSequence(foundProts.get(protAcc));

                            if (post == null || post.equals("")) {
                                post = "-";
                            }
                            pepEvid.setPost(post);

                            if (pre == null || pre.equals("")) {
                                pre = "-";
                            }
                            pepEvid.setPre(pre);
                            pepEvid.setId("PE" + sirCounter + "_" + siiCounter
                                    + "_" + pepEvidCounter);
                            pepEvidCounter++;
                            pepEvid.setIsDecoy(Boolean.FALSE);
                            if (this.decoyRegularExpression != null) {
                                if (protAcc
                                        .contains(this.decoyRegularExpression)) {
                                    pepEvid.setIsDecoy(Boolean.TRUE);
                                }
                            }
                            peptideEvidenceList.add(pepEvid);
                        } else {
                            pepEvid = pepEvidLookup.get(mzidPep.getId() + "_"
                                    + testProt);

                        }

                        PeptideEvidenceRef peptideEvidenceRef
                                = new PeptideEvidenceRef();
                        peptideEvidenceRef.setPeptideEvidence(pepEvid);
                        peptideEvidenceRefList.add(peptideEvidenceRef);
                        sii.setPeptide(mzidPep);
                    }
                }

                specIdentRes.setId("SIR_" + sirCounter);
                sirCounter++;
            }

        }

        System.out.println("Max evalue: " + maxEvalue);

    }

}
