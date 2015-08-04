package uk.ac.liv.mzidlib.converters;

import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import au.com.bytecode.opencsv.CSVReader;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.model.mzidml.*;

/**
 *
 * @author jonesar
 */
public class Csv2mzid {

    private MzIdentMLMarshaller marshaller;
    private String outFile = "temp.mzid";
    //private static String inputCsvFile = "example_files/Toxo_1D_Slice43_omssa.csv";
    private String parser_params = "example_files/toxo_omssa_params.csv";
    private String inputCsvFile = "example_files/Toxo_1D_Slice43_omssa_fiddle_ranks_for_testing.csv";
    //private URL configFile =  this.getClass().getResource("/resources/csv_config_file.csv");
    private InputStream inputConfigFile = ClassLoader.getSystemClassLoader().getResourceAsStream("csv_config_file.csv");
    private Map<String, String> dataTypeMap = new HashMap<>();         //mapping for config file from internal datatypes to datatypes specified by each Software
    private String softwareName = "";
    //These will be set after reading the config file - 
    private String cvAccForEvalue = "";
    private String cvNameForEvalue = "";
    private String cvAccForPvalue = "";
    private String cvNameForPvalue = "";
    private String cvScoreToOrderBy = "MS:1001328";     //This is a default only
    private Map<String, Peptide> peptideIDToPeptideMap = new HashMap<>();
    private Map<String, PeptideEvidence> idToPeptideEvidenceMap = new HashMap<>();
    private Map<String, SpectrumIdentificationResult> spectrumIDToSIRMap = new HashMap<>();
    //Param	cvTerm	Accession	Value  (structure of the parameter file containing search metadata)
    private Map<String, String> paramToCvParamName = new HashMap<>();
    private Map<String, String> paramToCvParamValue = new HashMap<>();
    private Map<String, String> paramToCvParamAcc = new HashMap<>();
    //Mods results name	Unimod name	Unimod ID	Residue	Fixed	Mass Delta   (structure of the parameter file containing the mods searched)
    private Map<String, String> modNameToUnimodName = new HashMap<>();
    private Map<String, String> modNameToUnimodID = new HashMap<>();
    private Map<String, String> modNameToResidue = new HashMap<>();
    private Map<String, String> modNameToPosition = new HashMap<>();
    private Map<String, Boolean> modNameToIsFixed = new HashMap<>();
    private Map<String, Float> modNameToMassDelta = new HashMap<>();
    private Map<String, DBSequence> accessionToDBSequenceHashMap = new HashMap<>();
    static String decoyRegex = null;
    ReadUnimod unimodDoc;
    static String siiListID = "SII_LIST_1";
    static String spectraDataID = "SID_1";
    static String psiCvID = "PSI-MS";
    static String siProtocolID = "SearchProtocol_1";
    static String searchDBID = "SearchDB_1";
    static String pepEvidenceListID = "PepEvidList_1";
    static String analysisSoftID = "ID_software";
    static String specIdentID = "SpecIdent_1";
    static String unimodID = "UNIMOD";
    static String unitCvID = "UO";
    static String measureMzID = "Measure_MZ";
    static String measureIntID = "Measure_Int";
    static String measureErrorID = "Measure_Error";
    static String sourceFileID = "SourceFile_1";
    //Some objects we will need globally
    private Cv unimodCV;
    private Cv psiCV;
    private Cv unitCV;
    private SpectrumIdentificationProtocol siProtocol;
    private SearchDatabase searchDB;
    private SpectraData spectraData;
    private Person docOwner;
    private SequenceCollection sequenceCollection;
    private SpectrumIdentificationList siList;
    private CvList cvList;
    private AnalysisSoftwareList analysisSoftwareList;
    private Provider provider;
    private AuditCollection auditCollection;
    private AnalysisSampleCollection analysisSampleCollection;
    private AnalysisCollection analysisCollection;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private Inputs inputs;
    private SpectrumIdentification specIdent;
    private List<SpectrumIdentification> specIdentList;
    private AnalysisSoftware analysisSoftware;
    Map<String, SpectraData> locationToSpectraDataMap = new HashMap<>();
    private boolean addFixedModsBasedOnSearchParams = true;
    public String regexForExtractingAccsFromDefline = " ";     //OMSSA often fails to extract the accessions from the defline, so we use a regex to grab the accession e.g. prior to the first space
    public boolean orderSIIRanksLowToHigh = true;               // The score we are using is ordered low to high. This should be changed for scores going high to low

    public static void main(String[] args) {

        //For internal testing only
        Csv2mzid csv2mzid = new Csv2mzid();
    }

    /*
     * Constructor for testing with in-built example files.
     */
    public Csv2mzid() {
        System.out.println("Running conversion for the example files...");

        this.init();
    }

    public Csv2mzid(String inputfile, String outputfile, String paramsFile, String cvParamAccForRankingPSMs, String decoyRegularExpression, boolean addFixedModificationsBasedOnSearchParams) {
        decoyRegex = decoyRegularExpression;
        this.inputCsvFile = inputfile;
        this.outFile = outputfile;
        this.parser_params = paramsFile;
        this.cvScoreToOrderBy = cvParamAccForRankingPSMs;
        this.addFixedModsBasedOnSearchParams = addFixedModificationsBasedOnSearchParams;

        this.init();
    }

    public Csv2mzid(String inputfile, String outputfile, String paramsFile, String cvParamAccForRankingPSMs, boolean addFixedModificationsBasedOnSearchParams) {
        this.inputCsvFile = inputfile;
        this.outFile = outputfile;
        this.parser_params = paramsFile;

        this.cvScoreToOrderBy = cvParamAccForRankingPSMs;
        this.addFixedModsBasedOnSearchParams = addFixedModificationsBasedOnSearchParams;
        this.init();

    }

    /*
     * Method to set going to conversion once all the global variables have been
     * set
     */
    private void init() {
        unimodDoc = new ReadUnimod();
        marshaller = new MzIdentMLMarshaller();

        buildParameters();
        readConfigFile();
        buildMetaDeta();
        validateMetadata();
        validateModifications();
        handleAnalysisCollection();
        buildPSMs();
        writeMzidFile();
    }

    private void readConfigFile() {
        try {
            CSVReader reader = new CSVReader(new InputStreamReader(inputConfigFile));
            String[] nextLine;

            int lineCounter = 0;
            boolean mappingCVTerms = false;

            int softwareColumn = -1;        //Needs to match header to know how data type mappings will work
            while ((nextLine = reader.readNext()) != null) {

                String firstCell = nextLine[0].trim();

                if (firstCell != null) {
                    if (lineCounter == 0) {
                        for (int i = 1; i < nextLine.length; i++) {
                            if (nextLine[i].equals(softwareName)) {
                                softwareColumn = i;
                                break;
                            }
                        }
                        if (softwareColumn == -1) {
                            System.out.println("No valid software mapping found from params file \"Software Name\":" + softwareName + " value to config file column header, exiting...");
                        }
                    } else if (firstCell.equals("CV mappings")) {     //Reached the part of the file that deals with the mods
                        mappingCVTerms = true;
                    } else {
                        if (!mappingCVTerms) {
                            //Just grab the column with the correct software name
                            dataTypeMap.put(nextLine[0], nextLine[softwareColumn]);

                        } else {
                            if (!firstCell.equals("Data type")) {

                                if (nextLine[0].equals("evalue") && nextLine[3].equals(softwareName)) {
                                    cvAccForEvalue = nextLine[2];
                                    cvNameForEvalue = nextLine[1];
                                } else if (nextLine[0].equals("pvalue") && nextLine[3].equals(softwareName)) {
                                    cvAccForPvalue = nextLine[2];
                                    cvNameForPvalue = nextLine[1];
                                }
                            }
                        }
                    }
                    lineCounter++;
                }
            }
            reader.close();

            if (cvAccForEvalue.equals("") || cvAccForPvalue.equals("") || cvNameForPvalue.equals("") || cvNameForEvalue.equals("")) {
                System.out.println("Error - not recognized e or p-value equivalent from config file");
            }

            System.out.println("From params file and config file, two PSM data types set as: " + cvAccForEvalue + " and " + cvAccForPvalue + "\nIf these are not present in your data file, an error will result");
        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }


    }

    /*
     * Method for reading the params file
     */
    private void buildParameters() {

        try {
            CSVReader reader = new CSVReader(new FileReader(parser_params));
            String[] nextLine;

            int lineCounter = 0;
            boolean parsingMods = false;
            while ((nextLine = reader.readNext()) != null) {

                String firstCell = nextLine[0];
                if (firstCell != null) {
                    if (lineCounter == 0) {
                        //ignore headers
                    } else if (firstCell.equals("Mods results name")) {     //Reached the part of the file that deals with the mods
                        parsingMods = true;
                    } else {
                        if (!parsingMods) {
                            paramToCvParamName.put(nextLine[0], nextLine[1]);
                            paramToCvParamAcc.put(nextLine[0], nextLine[2]);
                            paramToCvParamValue.put(nextLine[0], nextLine[3]);
                            //System.out.println("Read param: " + nextLine[0] + " " + nextLine[1] + " " + nextLine[2] + " " + nextLine[3]);
                        } else {
                            if (nextLine.length > 6) {
                                modNameToUnimodName.put(nextLine[0], nextLine[1]);
                                modNameToUnimodID.put(nextLine[0], nextLine[2]);
                                modNameToResidue.put(nextLine[0], nextLine[3]);
                                modNameToPosition.put(nextLine[0], nextLine[4]);
                                modNameToIsFixed.put(nextLine[0], Boolean.parseBoolean(nextLine[5]));
                                modNameToMassDelta.put(nextLine[0], Float.parseFloat(nextLine[6]));
                            } else {
                                System.out.println("Ignoring mod line since there are not the correct number of values: " + Arrays.toString(nextLine));
                            }
                        }
                    }
                    lineCounter++;
                }

            }
            reader.close();

            softwareName = paramToCvParamName.get("Software name");
            if (softwareName == null) {
                System.out.println("Error, params file did not contain a software name - required for configuring converter. Exiting..");
            }

        } catch (Exception e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 03 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    /*
     * The list of params MUST contain exactly these CV terms with a valid MS:
     * accession (only structure checked): mandatoryCVParams = {"Software name",
     * "Parent mass type","Fragment mass type","Enzyme","Fragment search
     * tolerance plus","Fragment search tolerance minus", "Parent search
     * tolerance plus","Parent search tolerance minus","PSM threshold","Input
     * file format","Database file format","Spectra data file format","Spectrum
     * ID format"};
     *
     * Also these params with a value MUST be present: otherMandatoryParams =
     * {"Local database path","Database name","Searched spectrum","Software
     * version","Missed cleavages","File contact first name", "File contact last
     * name","File contact organization name","File contact address"};
     *
     */
    private void validateMetadata() {

        System.out.print("Doing basic metadata validation...");
        String[] mandatoryCVParams = {"Software name",
            "Parent mass type", "Fragment mass type", "Enzyme", "Fragment search tolerance plus", "Fragment search tolerance minus",
            "Parent search tolerance plus", "Parent search tolerance minus", "PSM threshold", "Input file format", "Database file format", "Spectra data file format", "Spectrum ID format"};

        String[] otherMandatoryParams = {"Local database path", "Database name", "Searched spectrum", "Software version", "Missed cleavages", "File contact first name",
            "File contact last name", "File contact organization name", "File contact address"};

        boolean passed = true;
        for (String mandatoryCVParam : mandatoryCVParams) {
            if (!paramToCvParamName.containsKey(mandatoryCVParam)) {
                System.out.println("\nFATAL ERROR: Missing CV param in params file: " + mandatoryCVParam);
                passed = false;
            } else {
                String cvAcc = paramToCvParamAcc.get(mandatoryCVParam);
                if (cvAcc.equals("") || cvAcc == null || !cvAcc.contains("MS:")) {
                    System.out.println("\nFATAL ERROR: Missing CV accession or not a real PSI-MS term for " + mandatoryCVParam);
                    passed = false;
                }
            }
        }

        for (String otherMandatoryParam : otherMandatoryParams) {
            if (!paramToCvParamName.containsKey(otherMandatoryParam)) {
                System.out.println("\nFATAL ERROR: Missing  param in params file: " + otherMandatoryParam);
                passed = false;
            } else {
                String value = paramToCvParamValue.get(otherMandatoryParam);
                if (value.equals("") || value == null) {
                    System.out.println("\nFATAL ERROR: Missing value for mandatory param: " + otherMandatoryParam);
                    passed = false;
                }
            }
        }

        if (!passed) {
            System.out.println("Exiting, due to error in params specified in the param file, please consult the help guide");
        }

        System.out.print("...passed\n");
    }

    private void validateModifications() {


        System.out.print("Doing basic validation of modifications entered...");
        boolean passed = true;
        String[] allowedResidues = {"G", "P", "A", "V", "L", "I", "M", "C", "F", "Y", "W", "H", "K", "R", "Q", "N", "E", "D", "S", "T", "."};
        String[] allowedPositions = {"Peptide N-term", "Protein N-term", "Peptide C-term", "Protein N-term", "Any"};

        for (String modName : modNameToUnimodName.keySet()) {

            String unimodID = modNameToUnimodID.get(modName);
            String residue = modNameToResidue.get(modName);
            String position = modNameToPosition.get(modName);

            if (unimodID.equals("") || unimodID == null || !(unimodID.contains("UNIMOD:") || unimodID.equals("MS:1001460"))) {
                System.out.println("\nFATAL ERROR: Bad modification ID entered in the params file: " + unimodID
                        + ". Modification MUST be UNIMOD:[int] or \"MS:1001460\" - the term for an unknown modification");
                passed = false;
            }
            boolean residueValidated = false;
            for (String res : allowedResidues) {
                if (res.equals(residue)) {
                    residueValidated = true;
                }
            }
            if (!residueValidated) {
                System.out.println("\nFATAL ERROR in search modifications, illegal residue entered: " + residue
                        + " allowed residues:" + Arrays.toString(allowedResidues));
                passed = false;
            }

            boolean positionValidated = false;
            for (String pos : allowedPositions) {
                if (pos.equals(position)) {
                    positionValidated = true;
                }
            }
            if (!positionValidated) {
                System.out.println("\nFATAL ERROR in search modifications, illegal position entered: " + position
                        + " allowed positions:" + Arrays.toString(allowedResidues));
                passed = false;
            }
        }

        if (!passed) {
            System.out.println("Exiting, due to error in modifications specified in the param file, please consult the help guide");
        }

        System.out.print("...passed\n");
    }

    private void buildMetaDeta() {

        handleCVs();
        handleAnalysisSoftware();
        handleAuditCollection();
        handleProvider();
        handleAnalysisProtocolCollection();
        handleInputs();

    }

    /**
     *
     * Aim is to write out set up the analysisSoftwareList following this
     * structure: <AnalysisSoftware id="ID_software" name="xtandem"
     * version="2008.12.1.1" > <SoftwareName> <cvParam accession="MS:1001476"
     * name="xtandem" cvRef="PSI-MS" /> </SoftwareName>
     *
     */
    public void handleAnalysisSoftware() {
        analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        analysisSoftware = new AnalysisSoftware();

        String softwareVersion = paramToCvParamValue.get("Software version");
        analysisSoftware.setName(softwareName);
        Param tempParam = new Param();
        //tempParam.setParamGroup(makeCvParam("MS:1001475","OMSSA",psiCV));
        tempParam.setParam(makeCvParam(paramToCvParamAcc.get("Software name"), paramToCvParamName.get("Software name"), psiCV));
        analysisSoftware.setSoftwareName(tempParam);
        analysisSoftware.setId(analysisSoftID);
        analysisSoftware.setVersion(softwareVersion);

        analysisSoftwares.add(analysisSoftware);

    }

    public void handleCVs() {


        //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
        //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
        //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>

        cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        psiCV = new Cv();
        psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
        psiCV.setId(psiCvID);
        psiCV.setVersion("2.25.0");
        psiCV.setFullName("PSI-MS");

        unimodCV = new Cv();
        unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
        unimodCV.setId(unimodID);
        unimodCV.setFullName("UNIMOD");

        unitCV = new Cv();
        unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
        unitCV.setId(unitCvID);
        unitCV.setFullName("UNIT-ONTOLOGY");

        localCvList.add(psiCV);
        localCvList.add(unimodCV);
        localCvList.add(unitCV);
    }

    /**
     * Setup Provider element as follows <Provider id="PROVIDER"> <ContactRole
     * Contact_ref="PERSON_DOC_OWNER"> <role> <cvParam accession="MS:1001271"
     * name="researcher" cvRef="PSI-MS"/> </role> </ContactRole> </Provider>
     *
     */
    public void handleProvider() {
        provider = new Provider();
        provider.setId("PROVIDER");

        ContactRole contactRole = new ContactRole();
        contactRole.setContact(docOwner);


        Role role = new Role();
        role.setCvParam(makeCvParam("MS:1001271", "researcher", psiCV));
        contactRole.setRole(role);

        provider.setContactRole(contactRole);

    }

    /**
     * TO DO Capture name and email of the user <AuditCollection> <Person
     * id="PERSON_DOC_OWNER" firstName="Andy" lastName="Jones"
     * email="someone@someuniversity.com"> <affiliations
     * Organization_ref="ORG_DOC_OWNER"/> </Person> <Organization
     * id="ORG_DOC_OWNER" address="Some address" name="Some place" />
     * </AuditCollection>
     *
     *
     */
    public void handleAuditCollection() {
        auditCollection = new AuditCollection();
        //List<Contact> contactList = auditCollection.getContactGroup();
        List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
        docOwner = new Person();
        docOwner.setId("PERSON_DOC_OWNER");

        docOwner.setFirstName(paramToCvParamValue.get("File contact first name"));
        docOwner.setLastName(paramToCvParamValue.get("File contact last name"));
        docOwner.getCvParam().add(makeCvParam("MS:1000587", "contact address", psiCV, paramToCvParamValue.get("File contact address")));

        //docOwner.setEmail(email);

        Organization org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName(paramToCvParamValue.get("File contact organization name"));
        org.getCvParam().add(makeCvParam("MS:1000586", "contact name", psiCV, paramToCvParamValue.get("File contact organization name")));
        //org.setAddress(address);


        List<Affiliation> affList = docOwner.getAffiliation();
        Affiliation aff = new Affiliation();
        aff.setOrganization(org);
        affList.add(aff);
        contactList.add(docOwner);
        contactList.add(org);

    }

    /**
     * TODO This part is optional in the file - not yet completed
     *
     *
     *
     */
    public void handleAnalysisSampleCollection() {
        analysisSampleCollection = new AnalysisSampleCollection();

    }

    /**
     * <AnalysisCollection> <SpectrumIdentification id="SI_1"
     * SpectrumIdentificationProtocol_ref="SearchProtocol"
     * SpectrumIdentificationList_ref="siiListID"
     * activityDate="2008-02-27T08:22:12"> <InputSpectra
     * SpectraData_ref="SD_1"/> <SearchDatabase
     * SearchDatabase_ref="search_database"/> </SpectrumIdentification>
     * </AnalysisCollection>
     *
     */
    public void handleAnalysisCollection() {
        analysisCollection = new AnalysisCollection();
        specIdentList = analysisCollection.getSpectrumIdentification();
        specIdent = new SpectrumIdentification();
        specIdent.setId(specIdentID);
        specIdent.setSpectrumIdentificationProtocol(siProtocol);
        specIdentList.add(specIdent);
        List<SearchDatabaseRef> searchDBRefList = specIdent.getSearchDatabaseRef();
        SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
        searchDBRef.setSearchDatabase(searchDB);
        searchDBRefList.add(searchDBRef);



        if (paramToCvParamValue.get("Searched spectrum") != null) {
            List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
            InputSpectra inputSpec = new InputSpectra();
            inputSpec.setSpectraData(spectraData);
            inputSpecList.add(inputSpec);
            //specIdentList.add(specIdent);     //bug caused by this line - removed by ARJ 27/03/2013
        }

    }

    /**
     * <AnalysisProtocolCollection> <SpectrumIdentificationProtocol
     * id="SearchProtocol" AnalysisSoftware_ref="ID_software"> <SearchType>
     * <cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
     * </SearchType> <AdditionalSearchParams> <cvParam accession="MS:1001211"
     * name="parent mass type mono" cvRef="PSI-MS"/> <cvParam
     * accession="MS:1001256" name="fragment mass type mono" cvRef="PSI-MS"/>
     * </AdditionalSearchParams> <ModificationParams> <SearchModification
     * fixedMod="true"> <ModParam massDelta="57.021464" residues="C"> <cvParam
     * accession="UNIMOD:4" name="Carbamidomethyl" cvRef="UNIMOD" /> </ModParam>
     * </SearchModification> <SearchModification fixedMod="false"> <ModParam
     * massDelta="15.994919" residues="M"> <cvParam accession="UNIMOD:35"
     * name="Oxidation" cvRef="UNIMOD" /> </ModParam> </SearchModification>
     * </ModificationParams> <Enzymes independent="0"> <Enzyme id="ENZ_1"
     * CTermGain="OH" NTermGain="H" missedCleavages="1" semiSpecific="0">
     * <EnzymeName> <cvParam accession="MS:1001251" name="Trypsin"
     * cvRef="PSI-MS" /> </EnzymeName> </Enzyme> </Enzymes> <MassTable id="0"
     * msLevel="2"> </MassTable> <FragmentTolerance> <cvParam
     * accession="MS:1001412" name="search tolerance plus value" value="0.5"
     * cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
     * unitCvRef="UO" /> <cvParam accession="MS:1001413" name="search tolerance
     * minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
     * unitName="dalton" unitCvRef="UO" /> </FragmentTolerance>
     * <ParentTolerance> <cvParam accession="MS:1001412" name="search tolerance
     * plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
     * unitName="dalton" unitCvRef="UO" /> <cvParam accession="MS:1001413"
     * name="search tolerance minus value" value="0.5" cvRef="PSI-MS"
     * unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
     * </ParentTolerance> <Threshold> <cvParam accession="MS:1001494" name="no
     * threshold" cvRef="PSI-MS" /> </Threshold>
     * </SpectrumIdentificationProtocol> </AnalysisProtocolCollection>
     *
     *
     */
    public void handleAnalysisProtocolCollection() {
        analysisProtocolCollection = new AnalysisProtocolCollection();
        List<SpectrumIdentificationProtocol> sipList = analysisProtocolCollection.getSpectrumIdentificationProtocol();

        siProtocol = new SpectrumIdentificationProtocol();
        siProtocol.setId(siProtocolID);
        siProtocol.setAnalysisSoftware(analysisSoftware);

        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
        Param tempParam = new Param();
        tempParam.setParam(makeCvParam("MS:1001083", "ms-ms search", psiCV));
        siProtocol.setSearchType(tempParam);

        //List<CvParam> cvParamList = siProtocol.getAdditionalSearchCvParams();
        ParamList paramList = siProtocol.getAdditionalSearchParams();
        if (paramList == null) {
            paramList = new ParamList();
            siProtocol.setAdditionalSearchParams(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();
        cvParamList.add(makeCvParam(paramToCvParamName.get("Parent mass type"), paramToCvParamAcc.get("Parent mass type"), psiCV));
        cvParamList.add(makeCvParam(paramToCvParamName.get("Fragment mass type"), paramToCvParamAcc.get("Fragment mass type"), psiCV));

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams.getSearchModification();

        for (String modName : modNameToUnimodName.keySet()) {
            SearchModification searchMod = new SearchModification();
            searchMod.setFixedMod(modNameToIsFixed.get(modName));
            List<CvParam> modCvParamList = searchMod.getCvParam();
            modCvParamList.add(makeCvParam(modNameToUnimodID.get(modName), modNameToUnimodName.get(modName), unimodCV));
            searchMod.setMassDelta(modNameToMassDelta.get(modName));
            List<String> residueList = searchMod.getResidues();
            String modResidue = modNameToResidue.get(modName);
            residueList.add(modResidue);

            if (modNameToPosition.get(modName).equals("Peptide N-term")) {
                modCvParamList.add(makeCvParam("MS:1001189", "modification specificity peptide N-term", psiCV));
            } else if (modNameToPosition.get(modName).equals("Peptide C-term")) {
                modCvParamList.add(makeCvParam("MS:1001190", "modification specificity peptide C-term", psiCV));
            } else if (modNameToPosition.get(modName).equals("Protein N-term")) {
                modCvParamList.add(makeCvParam("MS:1002057", "modification specificity protein N-term", psiCV));

                if (addFixedModsBasedOnSearchParams) {
                    System.out.println("Please note - adding fixed mods on protein N-term not yet supported");
                }
            } else if (modNameToPosition.get(modName).equals("Protein C-term")) {
                modCvParamList.add(makeCvParam("MS:1002058", "modification specificity protein C-term", psiCV));
                if (addFixedModsBasedOnSearchParams) {
                    System.out.println("Please note - adding fixed mods on protein N-term not yet supported");
                }
            }

            searchModList.add(searchMod);
        }


        siProtocol.setModificationParams(modParams);

        Enzymes enzymes = siProtocol.getEnzymes();

        if (enzymes == null) {
            enzymes = new Enzymes();
            siProtocol.setEnzymes(enzymes);
        }
        enzymes.setIndependent(false);

        List<Enzyme> enzymeList = enzymes.getEnzyme();


        Enzyme enzyme = new Enzyme();
        //[KR]|{P}

        enzyme.setId("Enz1");
        enzyme.setCTermGain("OH");
        enzyme.setNTermGain("H");
        enzyme.setMissedCleavages(Integer.parseInt(paramToCvParamValue.get("Missed cleavages")));
        enzyme.setSemiSpecific(false);
        ParamList eParamList = enzyme.getEnzymeName();
        if (eParamList == null) {
            eParamList = new ParamList();
            enzyme.setEnzymeName(eParamList);
        }
        List<CvParam> eCvParamList = eParamList.getCvParam();
        eCvParamList.add(makeCvParam(paramToCvParamName.get("Enzyme"), paramToCvParamAcc.get("Enzyme"), psiCV));
        enzymeList.add(enzyme);
        Tolerance fragTol = new Tolerance();
        Tolerance parTol = new Tolerance();

        List<CvParam> fragCvList = fragTol.getCvParam();
        CvParam fragCvPlus = getCvParamWithMassUnits(true);
        CvParam fragCvMinus = getCvParamWithMassUnits(true);


        fragCvPlus.setAccession("MS:1001412");
        fragCvPlus.setName("search tolerance plus value");
        fragCvMinus.setAccession("MS:1001413");
        fragCvMinus.setName("search tolerance minus value");
        fragCvPlus.setValue(paramToCvParamValue.get("Fragment search tolerance plus"));
        fragCvMinus.setValue(paramToCvParamValue.get("Fragment search tolerance minus"));
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        List<CvParam> parCvList = parTol.getCvParam();
        CvParam parCvPlus = getCvParamWithMassUnits(true);
        CvParam parCvMinus = getCvParamWithMassUnits(true);

        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue(paramToCvParamValue.get("Parent search tolerance plus"));
        parCvMinus.setValue(paramToCvParamValue.get("Parent search tolerance plus"));
        parCvList.add(parCvPlus);
        parCvList.add(parCvMinus);

        siProtocol.setFragmentTolerance(fragTol);
        siProtocol.setParentTolerance(parTol);

        // siProtocol.getThresholdCvParams();
        ParamList sip_paramList = siProtocol.getThreshold();
        if (sip_paramList == null) {
            sip_paramList = new ParamList();
            siProtocol.setThreshold(sip_paramList);
        }
        cvParamList = sip_paramList.getCvParam();

        cvParamList.add(makeCvParam(paramToCvParamAcc.get("PSM threshold"), paramToCvParamName.get("PSM threshold"), psiCV));
        sipList.add(siProtocol);

    }

    public void handleInputs() {

        inputs = new Inputs();
        List<SearchDatabase> searchDBList = inputs.getSearchDatabase();

        searchDB = new SearchDatabase();
        searchDB.setId(searchDBID);

        UserParam param = new UserParam();
        //param.setName(settings.MSSearchSettings_db);
        Param tempParam = new Param();

        param.setName(paramToCvParamValue.get("Database name"));
        tempParam.setParam(param);
        searchDB.setDatabaseName(tempParam);
        searchDB.setLocation(paramToCvParamValue.get("Local database path"));

        List<CvParam> searchDBCvParamList = searchDB.getCvParam();

        if (paramToCvParamName.get("Decoy database composition") != null) {
            searchDBCvParamList.add(makeCvParam(paramToCvParamAcc.get("Decoy database composition"), paramToCvParamName.get("Decoy database composition"), psiCV));
        }
        if (paramToCvParamName.get("Decoy database regex") != null) {
            searchDBCvParamList.add(makeCvParam(paramToCvParamAcc.get("Decoy database regex"), paramToCvParamName.get("Decoy database regex"), psiCV, paramToCvParamValue.get("Decoy database regex")));
        }
        if (paramToCvParamName.get("Decoy database type") != null) {
            searchDBCvParamList.add(makeCvParam(paramToCvParamAcc.get("Decoy database type"), paramToCvParamName.get("Decoy database type"), psiCV));
        }

        FileFormat ff = new FileFormat();
        ff.setCvParam(makeCvParam(paramToCvParamAcc.get("Database file format"), paramToCvParamName.get("Database file format"), psiCV));
        searchDB.setFileFormat(ff);
        searchDBList.add(searchDB);

        List<SourceFile> sourceFileList = inputs.getSourceFile();
        SourceFile sourceFile = new SourceFile();
        sourceFile.setLocation(inputCsvFile);
        sourceFile.setId(sourceFileID);
        ff = new FileFormat();
        ff.setCvParam(makeCvParam(paramToCvParamAcc.get("Input file format"), paramToCvParamName.get("Input file format"), psiCV));

        sourceFile.setFileFormat(ff);
        sourceFileList.add(sourceFile);

        if (paramToCvParamValue.get("Searched spectrum") != null) {
            List<SpectraData> spectraDataList = inputs.getSpectraData();
            spectraData = new SpectraData();
            SpectrumIDFormat sif = new SpectrumIDFormat();
            sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", psiCV));
            spectraData.setSpectrumIDFormat(sif);

            ff = new FileFormat();
            ff.setCvParam(makeCvParam(paramToCvParamAcc.get("Spectra data file format"), paramToCvParamName.get("Spectra data file format"), psiCV));
            spectraData.setFileFormat(ff);

            spectraData.setId(spectraDataID);
            spectraData.setLocation(paramToCvParamValue.get("Searched spectrum"));
            spectraDataList.add(spectraData);
        }


    }

    private void buildPSMs() {

        try {
            CSVReader reader = new CSVReader(new FileReader(inputCsvFile));
            System.out.println("Processing..." + inputCsvFile);
            String[] nextLine;


            sequenceCollection = new SequenceCollection();
            List<Peptide> peptideList = sequenceCollection.getPeptide();

            siList = new SpectrumIdentificationList();
            siList.setId(siiListID);
            specIdent.setSpectrumIdentificationList(siList);
            List<DBSequence> dbSequenceList = sequenceCollection.getDBSequence();
            List<PeptideEvidence> peptideEvidenceList = sequenceCollection.getPeptideEvidence();
            AnalysisData analysisData = new AnalysisData();

            List<SpectrumIdentificationResult> sirList = siList.getSpectrumIdentificationResult();

            Map<String, Integer> headerToColumnMap = new HashMap<>();
            Map<Integer, String> columnToHeaderMap = new HashMap<>();
            int lineCounter = 0;

            int sirCounter = 1;
            int specDataCounter = 1;  //Only required if spectrum location is in the searched file (not in a typical Omssa file)

            while ((nextLine = reader.readNext()) != null) {

                if (lineCounter == 0) {
                    for (int i = 0; i < nextLine.length; i++) {

                        headerToColumnMap.put(nextLine[i].trim(), i);
                        columnToHeaderMap.put(i, nextLine[i]);
                    }
                } else {

                    String specUniqueID = null;
                    String spectrumID = "index=" + nextLine[headerToColumnMap.get("Spectrum number")];
                    String spectraLocation = null;

                    if (dataTypeMap.get("spectrum_location") != null) {
                        if (headerToColumnMap.get(dataTypeMap.get("spectrum_location")) != null) {
                            spectraLocation = nextLine[headerToColumnMap.get(dataTypeMap.get("spectrum_location"))];
                        }
                    }


                    if (spectraLocation != null && !spectraLocation.equals("")) {
                        specUniqueID = spectrumID + "_" + spectraLocation;
                    } else {
                        specUniqueID = spectrumID;
                    }

                    SpectrumIdentificationResult sir;
                    if (spectrumIDToSIRMap.containsKey(specUniqueID)) {
                        sir = spectrumIDToSIRMap.get(specUniqueID);
                    } else {
                        sir = new SpectrumIdentificationResult();
                        sirList.add(sir);
                        sir.setId("SIR_" + sirCounter);
                        sir.setSpectrumID(spectrumID);
                        sirCounter++;


                        //If spectra Location is found in the file (not Omssa)
                        if (spectraLocation != null) {
                            spectraData = locationToSpectraDataMap.get(spectraLocation);

                            //Create a new spectra data object and add input spectra to the SpectrumIdentification object
                            if (spectraData == null) {
                                List<SpectraData> spectraDataList = inputs.getSpectraData();
                                spectraData = new SpectraData();
                                SpectrumIDFormat sif = new SpectrumIDFormat();
                                sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", psiCV));
                                spectraData.setSpectrumIDFormat(sif);

                                FileFormat ff = new FileFormat();
                                ff.setCvParam(makeCvParam(paramToCvParamAcc.get("Spectra data file format"), paramToCvParamName.get("Spectra data file format"), psiCV));
                                spectraData.setFileFormat(ff);

                                spectraDataID = "SD_" + specDataCounter;
                                specDataCounter++;
                                spectraData.setId(spectraDataID);
                                spectraData.setLocation(spectraLocation);
                                spectraDataList.add(spectraData);
                                List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
                                InputSpectra inputSpec = new InputSpectra();
                                inputSpec.setSpectraData(spectraData);
                                inputSpecList.add(inputSpec);
                                //specIdentList.add(specIdent);
                                locationToSpectraDataMap.put(spectraLocation, spectraData);
                            }
                        }
                        sir.setSpectraData(spectraData);

                        spectrumIDToSIRMap.put(specUniqueID, sir);

                    }

                    String defline = nextLine[headerToColumnMap.get("Defline")];

                    if (defline == null) {
                        System.out.println("Error - Unable to extract protein accessions from the Defline - quitting");
                    }
                    String protAcc = null;

                    if (defline.indexOf(regexForExtractingAccsFromDefline) != -1) { //Set this is a regex to grab the correct accession. Omssa often fails to grab the accession but we can usually get it from the defline
                        protAcc = defline.substring(0, defline.indexOf(regexForExtractingAccsFromDefline));
                    } else {
                        protAcc = defline;      //Try grabbing the whole defline
                    }

                    if (protAcc == null) {
                        System.out.println("Error - Unable to extract protein accessions from the Defline - quitting");
                    }


                    DBSequence dbSequence = null;
                    if (accessionToDBSequenceHashMap.containsKey(protAcc)) {
                        dbSequence = accessionToDBSequenceHashMap.get(protAcc);
                    } else {
                        dbSequence = new DBSequence();
                        dbSequence.setId("dbseq_" + protAcc);
                        dbSequence.setAccession(protAcc);
                        dbSequenceList.add(dbSequence);
                        dbSequence.getCvParam().add(makeCvParam("MS:1001088", "protein description", psiCV, defline));
                        dbSequence.setSearchDatabase(searchDB);
                        accessionToDBSequenceHashMap.put(protAcc, dbSequence);
                    }


                    String pepSeq = nextLine[headerToColumnMap.get("Peptide")].toUpperCase();
                    String modString = nextLine[headerToColumnMap.get("Mods")];

                    Peptide pep;
                    String pepID = pepSeq + "_" + modString;

                    if (peptideIDToPeptideMap.containsKey(pepID)) {
                        pep = peptideIDToPeptideMap.get(pepID);
                    } else {
                        pep = new Peptide();
                        peptideList.add(pep);
                        pep.setPeptideSequence(pepSeq);
                        pep.setId(pepID);
                        peptideIDToPeptideMap.put(pepID, pep);

                        convertVarMods(pep, modString);

                        if (addFixedModsBasedOnSearchParams) {
                            convertFixedMods(pep);
                        }


                    }

                    int start = Integer.parseInt(nextLine[headerToColumnMap.get("Start")]);
                    int end = Integer.parseInt(nextLine[headerToColumnMap.get("Stop")]);
                    String pepEvidID = pepSeq + "_" + protAcc + "_" + start + "_" + end;

                    PeptideEvidence peptideEvidence = null;
                    if (idToPeptideEvidenceMap.containsKey(pepEvidID)) {
                        peptideEvidence = idToPeptideEvidenceMap.get(pepEvidID);

                    } else {
                        peptideEvidence = new PeptideEvidence();
                        peptideEvidence.setId(pepEvidID);
                        peptideEvidence.setDBSequence(dbSequence);
                        peptideEvidence.setPeptide(pep);

                        peptideEvidence.setStart(start);
                        peptideEvidence.setEnd(end);

                        peptideEvidence.setIsDecoy(Boolean.FALSE);
                        if (decoyRegex != null) {
                            if (protAcc.contains(decoyRegex)) {
                                peptideEvidence.setIsDecoy(Boolean.TRUE);
                            }
                        }
                        peptideEvidenceList.add(peptideEvidence);

                        idToPeptideEvidenceMap.put(pepEvidID, peptideEvidence);

                    }



                    List<SpectrumIdentificationItem> siiList = sir.getSpectrumIdentificationItem();

                    /*
                     * Cases: If PepID is the same, then this is another
                     * PeptideEvidence, use same sii If PepID is different, new
                     * SII
                     */

                    SpectrumIdentificationItem sii = null;

                    for (SpectrumIdentificationItem currentSii : siiList) {
                        String currentPepID = currentSii.getPeptideRef();
                        if (currentPepID.equals(pepID)) {
                            sii = currentSii;
                            break;
                        }
                    }

                    //Now create a new SII
                    if (sii == null) {
                        sii = new SpectrumIdentificationItem();
                        sii.setChargeState(Integer.parseInt(nextLine[headerToColumnMap.get("Charge")]));
                        sii.setExperimentalMassToCharge(Double.parseDouble(nextLine[headerToColumnMap.get("Mass")]));
                        sii.setCalculatedMassToCharge(Double.parseDouble(nextLine[headerToColumnMap.get("Theo Mass")]));
                        List<CvParam> cvParamList = sii.getCvParam();
                        cvParamList.add(makeCvParam(cvAccForEvalue, cvNameForEvalue, psiCV, nextLine[headerToColumnMap.get(dataTypeMap.get("evalue"))]));
                        cvParamList.add(makeCvParam(cvAccForPvalue, cvNameForPvalue, psiCV, nextLine[headerToColumnMap.get(dataTypeMap.get("pvalue"))]));
                        sii.setPeptide(pep);
                        addSIIToListAndSetRank(siiList, sii, cvScoreToOrderBy, orderSIIRanksLowToHigh, sir.getId());
                    }
                    PeptideEvidenceRef pepEvidRef = new PeptideEvidenceRef();
                    pepEvidRef.setPeptideEvidence(peptideEvidence);
                    sii.getPeptideEvidenceRef().add(pepEvidRef);
                }
                lineCounter++;
            }
            reader.close();
        } catch (FileNotFoundException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 01 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);

        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    private void convertVarMods(Peptide pep, String modString) {
        List<Modification> modList = pep.getModification();
        String[] mods = modString.split(",");

        for (String mod : mods) {
            Modification mzidMod = new Modification();
            mod = mod.trim();
            String[] temp = mod.split(":");
            String oneModString = temp[0].trim();
            if (!oneModString.equals("")) {
                if (temp.length == 2) {
                    int location = Integer.parseInt(temp[1]);
                    Boolean foundOkay = true;
                    if (modNameToMassDelta.get(oneModString) != null) {
                        mzidMod.setMonoisotopicMassDelta((double) modNameToMassDelta.get(oneModString));
                    } else {
                        System.out.println("Unable to insert correct mass delta:" + oneModString);
                        foundOkay = false;
                    }

                    mzidMod.setLocation(location);
                    String residue = modNameToResidue.get(oneModString);
                    if (residue != null) {
                        if (residue.equals("Peptide N-term")) {
                            mzidMod.getCvParam().add(makeCvParam("MS:1001189", "modification specificity peptide N-term", psiCV));
                            //mzidMod.getResidues().add(residue); //don't set residues for terminal mods
                        } else if (residue.equals("Peptide C-term")) {
                            mzidMod.getCvParam().add(makeCvParam("MS:1001190", "modification specificity peptide C-term", psiCV));
                            //mzidMod.getResidues().add(residue);
                        } else {
                            mzidMod.getResidues().add(residue);
                        }
                    } else {
                        System.out.println("Unable to insert correct residue for:" + oneModString);
                        foundOkay = false;
                    }

                    if (foundOkay) {
                        CvParam modParam = new CvParam();
                        modParam.setAccession(modNameToUnimodID.get(oneModString));
                        modParam.setCv(unimodCV);
                        modParam.setName(modNameToUnimodName.get(oneModString));

                        mzidMod.getCvParam().add(modParam);
                        modList.add(mzidMod);
                    } else {
                        mzidMod.getCvParam().add(makeCvParam("MS:1001460", "unknown modification", psiCV, mod));
                    }
                } else {
                    System.out.println("Incorrectly formatted mod:" + mod);
                }
            }
        }

    }

    /*
     * Helper method to add fixed mods based on the search parameters i.e. for
     * OMSSA fixed mods are not reported, so this method must be used to insert
     * the correct mods
     *
     */
    private void convertFixedMods(Peptide mzidPep) {

        String pepSeq = mzidPep.getPeptideSequence();
        List<uk.ac.ebi.jmzidml.model.mzidml.Modification> allMods = mzidPep.getModification();

        for (String modName : modNameToUnimodName.keySet()) {

            if (modNameToIsFixed.get(modName)) {

                double monoMass = modNameToMassDelta.get(modName);

                String modifiedResidue = modNameToResidue.get(modName);
                int index = pepSeq.indexOf(modifiedResidue);

                if (!modifiedResidue.equals(".")) {   //This char is only allowed for N or C terminal mods
                    while (index != -1) {
                        uk.ac.ebi.jmzidml.model.mzidml.Modification mzidmod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
                        List<CvParam> paramList = mzidmod.getCvParam();
                        CvParam modParam = makeCvParam(modNameToUnimodID.get(modName), modNameToUnimodName.get(modName), unimodCV);
                        mzidmod.setMonoisotopicMassDelta(monoMass);
                        int mzidModLocation = index + 1;              //If second res is modified, index would return 1, but mzid position should be 2
                        mzidmod.setLocation(mzidModLocation);
                        mzidmod.getResidues().add(modifiedResidue);
                        boolean isMono = true;

                        paramList.add(modParam);
                        allMods.add(mzidmod);
                        index = pepSeq.indexOf(modifiedResidue, index + 1);
                    }
                } else {
                    if (modNameToPosition.get(modName).equals("Peptide N-term")) {

                        if (modifiedResidue.equals(".") || modifiedResidue.equals(pepSeq.substring(0, 1))) {
                            uk.ac.ebi.jmzidml.model.mzidml.Modification mzidmod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
                            List<CvParam> paramList = mzidmod.getCvParam();
                            CvParam modParam = makeCvParam(modNameToUnimodID.get(modName), modNameToUnimodName.get(modName), unimodCV);
                            mzidmod.setMonoisotopicMassDelta(monoMass);
                            mzidmod.setLocation(0);
                            paramList.add(modParam);
                            allMods.add(mzidmod);
                        }
                    } else if (modNameToPosition.get(modName).equals("Peptide C-term")) {
                        if (modifiedResidue.equals(".") || modifiedResidue.equals(pepSeq.substring(pepSeq.length() - 1, pepSeq.length()))) {
                            uk.ac.ebi.jmzidml.model.mzidml.Modification mzidmod = new uk.ac.ebi.jmzidml.model.mzidml.Modification();
                            List<CvParam> paramList = mzidmod.getCvParam();
                            CvParam modParam = makeCvParam(modNameToUnimodID.get(modName), modNameToUnimodName.get(modName), unimodCV);
                            mzidmod.setMonoisotopicMassDelta(monoMass);
                            mzidmod.setLocation(pepSeq.length() + 1);
                            paramList.add(modParam);
                            allMods.add(mzidmod);
                        }
                    }
                }


                //TODO Protein N or C term mods not yet supported
            }

        }

    }

    /*
     * Accepts the SIIs associated with any SIR and inserts the correct rank
     * values
     *
     */
    private void addSIIToListAndSetRank(List<SpectrumIdentificationItem> siiList, SpectrumIdentificationItem sii, String scoreCvParamNameToOrderBy, boolean orderLowToHigh, String sirID) {

        double scoreOfNewSII = getScoreFromSII(sii, scoreCvParamNameToOrderBy);
        final String cvParamScore = scoreCvParamNameToOrderBy;
        final boolean lowToHigh = orderLowToHigh;

        siiList.add(sii);
        Collections.sort(siiList, new Comparator<SpectrumIdentificationItem>() {

            @Override
            public int compare(SpectrumIdentificationItem sii1, SpectrumIdentificationItem sii2) {
                double sii1Score = getScoreFromSII(sii1, cvParamScore);
                double sii2Score = getScoreFromSII(sii2, cvParamScore);
                int i = 0;

                if (lowToHigh) {
                    if (sii1Score < sii2Score) {
                        i = -1;
                    } else if (sii1Score > sii2Score) {
                        i = +1;
                    } else {
                        i = 0;
                    }
                } else {
                    if (sii2Score < sii1Score) {
                        i = -1;
                    } else if (sii2Score > sii1Score) {
                        i = +1;
                    } else {
                        i = 0;
                    }
                }
                return i;
            }
        });

        int rank = 0;
        int innerRankCounter = 2;
        double lastScore = -999.0;

        for (SpectrumIdentificationItem spii : siiList) {
            double score = getScoreFromSII(spii, scoreCvParamNameToOrderBy);
            boolean sameScore = false;
            if (score != lastScore) { //Otherwise set the same rank as previous
                rank++;
                sameScore = true;
                spii.setId(sirID + "_SII_" + rank);
                spii.setRank(rank);
                innerRankCounter = 2;
            } else {
                spii.setId(sirID + "SII_" + rank + "_" + innerRankCounter);
                spii.setRank(rank);
                innerRankCounter++;
            }
        }
    }

    /*
     * Helper method to retrive a particular score type from an SII, based on
     * the CV accession. Will return 0 if the accession is not found.
     */
    public double getScoreFromSII(SpectrumIdentificationItem sii, String cvParamAccForScore) {

        double score = 0.0;

        for (CvParam cvParam : sii.getCvParam()) {
            if (cvParam.getAccession().equals(cvParamAccForScore)) {
                score = Double.parseDouble(cvParam.getValue());
            }
        }

        return score;

    }

    /**
     * Helper method to create and return a CvParam from accession, name and CV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name and CV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setValue(value);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name, CV,
     * unitAccession and unitName (unitCV is automatically provided)
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setUnitAccession(unitAccession);
        cvParam.setUnitCv(unitCV);
        cvParam.setUnitName(unitName);
        return cvParam;
    }

    /**
     * Helper method to create and return a CvParam from accession, name, CV,
     * unitAccession, unitName and unitCV
     *
     * @return CvParam
     */
    public CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName, Cv alternateUnitCV) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setUnitAccession(unitAccession);
        cvParam.setUnitCv(alternateUnitCV);
        cvParam.setUnitName(unitName);
        return cvParam;
    }

    public void writeMzidFile() {

        try {
            Writer writer = new FileWriter(outFile);


            // mzIdentML
            //     cvList
            //     AnalysisSoftwareList
            //     Provider
            //     AuditCollection
            //     AnalysisSampleCollection
            //     SequenceCollection
            //     AnalysisCollection
            //     AnalysisProtocolCollection
            //     DataCollection
            //         Inputs
            //         AnalysisData
            //             SpectrumIdentificationList
            //             ProteinDetectionList
            //         /AnalysisData
            //     /DataCollection
            //     BibliographicReference
            // /mzIdentML


            // Note: writing of '\n' characters is optional and only for readability of the produced XML document
            // Also note: since the XML is produced in individual parts, the overall formatting of the document
            //            is not as nice as it would be when marshalling the whole structure at once.

            // XML header
            writer.write(marshaller.createXmlHeader() + "\n");


            // mzIdentML start tag

            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");



            marshaller.marshall(cvList, writer);
            writer.write("\n");
            AnalysisSoftware analysisSoftware = new AnalysisSoftware();
            Date date = new Date() ;
            SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss") ;
            analysisSoftware.setName(this.getClass().getSimpleName()+"_"+dateFormat.format(date)); analysisSoftware.setId(this.getClass().getSimpleName()+"_"+dateFormat.format(date));

            marshaller.marshall(analysisSoftwareList, writer);
            writer.write("\n");


            marshaller.marshall(provider, writer);
            writer.write("\n");


            marshaller.marshall(auditCollection, writer);
            writer.write("\n");

            marshaller.marshall(sequenceCollection, writer);
            writer.write("\n");

            marshaller.marshall(analysisCollection, writer);
            writer.write("\n");


            marshaller.marshall(analysisProtocolCollection, writer);
            writer.write("\n");


            writer.write(marshaller.createDataCollectionStartTag() + "\n");
            marshaller.marshall(inputs, writer);
            writer.write("\n");


            //Inputs inputs = unmarshaller.unmarshal(MzIdentMLElement.Inputs.getXpath());
            //m.marshall(inputs, writer);
            //writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");



            // writer.write(m.createSpectrumIdentificationListStartTag("SIL_1", null, 71412L) + "\n");

            //FragmentationTable table = unmarshaller.unmarshal(MzIdentMLElement.FragmentationTable.getXpath());
            //m.marshall(table, writer);
            //writer.write("\n");


            //Iterator<SpectrumIdentificationResult> specResIter = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);

            /*
             * Iterator<SpectrumIdentificationResult> specResIter =
             * specIdentResults.iterator(); while (specResIter.hasNext()) {
             * SpectrumIdentificationResult specIdentRes = specResIter.next();
             * m.marshall(specIdentRes, writer); writer.write("\n"); }
             */

            marshaller.marshall(siList, writer);
            writer.write("\n");


            // writer.write(m.createSpectrumIdentificationListClosingTag() + "\n");

            writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");

            /*
             * Iterator<ProteinAmbiguityGroup> protAmbGroupIter =
             * unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
             * while (protAmbGroupIter.hasNext()) { ProteinAmbiguityGroup
             * protAmbGroup = protAmbGroupIter.next(); m.marshall(protAmbGroup,
             * writer); writer.write("\n"); }
             *
             */

            writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            //BibliographicReference ref = unmarshaller.unmarshal(MzIdentMLElement.BibliographicReference.getXpath());
            // m.marshall(ref, writer);
            // writer.write("\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    /**
     * Helper method to setup a CvParam with CVRef, with either Daltons or ppm
     * as units
     *
     */
    public CvParam getCvParamWithMassUnits(Boolean isDaltonUnit) {
        CvParam cvParam = new CvParam();

        //<cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
        cvParam.setCv(psiCV);
        cvParam.setUnitCv(unitCV);

        if (isDaltonUnit) {
            cvParam.setUnitAccession("UO:0000221");
            cvParam.setUnitName("dalton");
        } else {
            cvParam.setUnitAccession("UO:0000169");
            cvParam.setUnitName("parts per million");
        }
        return cvParam;
    }
}
