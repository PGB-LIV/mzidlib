
package uk.ac.liv.mzidlib.converters;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import de.proteinms.xtandemparser.xtandem.Domain;
import de.proteinms.xtandemparser.xtandem.FragmentIon;
import de.proteinms.xtandemparser.xtandem.InputParams;
import de.proteinms.xtandemparser.xtandem.PerformParams;
import de.proteinms.xtandemparser.xtandem.Protein;
import de.proteinms.xtandemparser.xtandem.Spectrum;
import de.proteinms.xtandemparser.xtandem.SupportData;
import de.proteinms.xtandemparser.xtandem.XTandemFile;
import uk.ac.ebi.jmzidml.model.mzidml.AbstractContact;
import uk.ac.ebi.jmzidml.model.mzidml.Affiliation;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
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
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.liv.mzidlib.util.CVUtils;
import uk.ac.liv.mzidlib.util.MzidLibUtils;
import uk.ac.liv.mzidlib.util.Utils;
import uk.ac.liv.unimod.ModT;

/**
 *
 * @author jonesar
 */
public class Tandem2mzid {

    //Constants:
    private static final String XTANDEM_REVERSED_FLAG = ":reversed";
    //These are the main structures to be output by the main writing method
    SequenceCollection sequenceCollection;
    SpectrumIdentificationList siList;
    CvList cvList;
    AnalysisSoftwareList analysisSoftwareList;
    Provider provider;
    AuditCollection auditCollection;
    AnalysisSampleCollection analysisSampleCollection;
    AnalysisCollection analysisCollection;
    AnalysisProtocolCollection analysisProtocolCollection;
    Inputs inputs;
    //Some IDs to be used throughout;
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
    String decoyRegex = null;    //Added by ARJ for setting is decoy
    //Some objects we will need globally
    Cv unimodCV;
    Cv psiCV;
    Cv unitCV;
    SpectrumIdentificationProtocol siProtocol;
    SearchDatabase searchDB;
    SpectraData spectraData;
    Person docOwner;
    AnalysisSoftware analysisSoftware;
    Map<String, String> cvMap = null;
    Map<String, Peptide> uniquePeps;
    //int pepCounter = 0;
    int pepEvidCounter = 0;
    //List<SpectrumIdentificationResult> specIdentResults;
    ReadUnimod unimodDoc;
    double unimodMassError = 0.01;              //TODO This parameter is hard-coded (ARJ changed from 0.001 to 0.01 - Aug2012; perhaps should be set dynamically from search params)
    //defaults:
    private boolean isMs2SpectrumIdStartingAtZero = false;
    private String databaseFileFormatID;
    private String databaseFileFormatName;
    private String massSpecFileFormatID;
    private String massSpecFileFormatName;
    private Pattern proteinCodeRegexPattern = null;
    boolean outputFragmentation = true;
    MzIdentMLVersion mzidVer;

    public Tandem2mzid(String inputfile, String outputfile, MzIdentMLVersion ver)
            throws Exception {
        //Calling the more detailed initialization method with all optionals set to null:
        XTandemFile xfile = new XTandemFile(inputfile);
        initializeVariables(xfile, null, null, null, null, null, true);
        convertFile(xfile, inputfile, outputfile, ver);
    }

    public Tandem2mzid(String inputfile, String outputfile,
                       String databaseFileFormatID, String massSpecFileFormatID,
                       boolean isMs2SpectrumIdStartingAtZero,
                       boolean outputFragmentation, MzIdentMLVersion ver)
            throws Exception {
        XTandemFile xfile = new XTandemFile(inputfile);
        initializeVariables(xfile, databaseFileFormatID, massSpecFileFormatID,
                            isMs2SpectrumIdStartingAtZero, null, null,
                            outputFragmentation);
        convertFile(xfile, inputfile, outputfile, ver);
    }

    public Tandem2mzid(String inputfile, String outputfile,
                       String databaseFileFormatID, String massSpecFileFormatID,
                       boolean isMs2SpectrumIdStartingAtZero,
                       String decoyRegularExpression,
                       boolean outputFragmentation, MzIdentMLVersion ver)
            throws Exception {
        XTandemFile xfile = new XTandemFile(inputfile);
        initializeVariables(xfile, databaseFileFormatID, massSpecFileFormatID,
                            isMs2SpectrumIdStartingAtZero,
                            decoyRegularExpression, null, outputFragmentation);
        convertFile(xfile, inputfile, outputfile, ver);
    }

    public Tandem2mzid(String inputfile, String outputfile,
                       String databaseFileFormatID, String massSpecFileFormatID,
                       boolean isMs2SpectrumIdStartingAtZero,
                       String decoyRegularExpression,
                       String proteinCodeRegex,
                       boolean outputFragmentation,
                       MzIdentMLVersion ver)
            throws Exception {
        XTandemFile xfile = new XTandemFile(inputfile);
        initializeVariables(xfile, databaseFileFormatID, massSpecFileFormatID,
                            isMs2SpectrumIdStartingAtZero,
                            decoyRegularExpression, proteinCodeRegex,
                            outputFragmentation);
        convertFile(xfile, inputfile, outputfile, ver);
    }

    private void initializeVariables(XTandemFile xfile,
                                     String databaseFileFormatID,
                                     String massSpecFileFormatID,
                                     Boolean isMs2SpectrumIdStartingAtZero,
                                     String decoyRegularExpression,
                                     String proteinCodeRegex,
                                     boolean outputFragmentation)
            throws IOException {
        if (decoyRegularExpression != null && !decoyRegularExpression.trim().
                equals("")) {
            this.decoyRegex = decoyRegularExpression;
        }

        if (proteinCodeRegex != null && !proteinCodeRegex.trim().equals("")) {
            //this regex should ensure the protein code is parsed from the longer string that X!Tandem is
            //currently making of this. See also https://code.google.com/p/mzidentml-lib/issues/detail?id=14
            try {
                this.proteinCodeRegexPattern = Pattern.compile(proteinCodeRegex);
            } catch (PatternSyntaxException pe) {
                throw new RuntimeException("Given proteinCodeRegex ["
                        + proteinCodeRegex + "] is not a "
                        + "valid regular expression. ", pe);
            }

        }

        this.outputFragmentation = outputFragmentation;

        if (databaseFileFormatID != null) {
            this.databaseFileFormatID = databaseFileFormatID;
            this.databaseFileFormatName = getCVName(databaseFileFormatID);
        } else {
            //Try to infer from the file itself:
            PerformParams tandemParams = xfile.getPerformParameters();
            String dbLocation = tandemParams.getSequenceSource_1();

            String[] cvIdAndName = CVUtils.getDatabaseFileFormat(dbLocation);
            this.databaseFileFormatID = cvIdAndName[0];
            this.databaseFileFormatName = cvIdAndName[1];
        }
        if (massSpecFileFormatID != null) {
            this.massSpecFileFormatID = massSpecFileFormatID;
            this.massSpecFileFormatName = getCVName(massSpecFileFormatID);
        } else {
            //Try to infer from the file itself:
            InputParams inputParams = xfile.getInputParameters();
            //Validate: if the spectrum path is null, then we can assume all input parameters are missing
            //as the spectrum path is a mandatory for X!tandem to run:
            String spectrumFile = inputParams.getSpectrumPath();
            if (spectrumFile == null) {
                throw new RuntimeException(
                        "Expected parameter 'spectrum, path' not found in X!Tandem file. Please run your X!Tandem search with option 'output, parameters=yes'. See http://thegpm.org/tandem/api/opara.html for more details.");
            }

            String[] cvIdAndName = CVUtils.getMassSpecFileFormatID(spectrumFile);
            this.massSpecFileFormatID = cvIdAndName[0];
            this.massSpecFileFormatName = cvIdAndName[1];

        }

        if (isMs2SpectrumIdStartingAtZero == null) {
            //if file format is mzML (MS:1000584), then spectrum starts at 0, otherwise it starts at 1
            if (this.massSpecFileFormatID.equalsIgnoreCase("MS:1000584")) {
                this.isMs2SpectrumIdStartingAtZero = true;
            } else {
                this.isMs2SpectrumIdStartingAtZero = false;
            }
        } else {
            this.isMs2SpectrumIdStartingAtZero = isMs2SpectrumIdStartingAtZero;
        }
    }

    private String getCVName(String cvItemID)
            throws IOException {
        //If CV map is not yet initialized, do it:
        if (this.cvMap == null) {
            this.cvMap = Utils.getInitializedCVMap();
        }

        //validate:
        if (this.cvMap.get(cvItemID) == null) {
            throw new RuntimeException(
                    "Given item not found in Controlled Vocabulary : "
                    + cvItemID);
        } else {
            return this.cvMap.get(cvItemID);
        }
    }

    private void convertFile(XTandemFile xfile, String inputfile,
                             String outputfile, MzIdentMLVersion ver)
            throws Exception {

        unimodDoc = new ReadUnimod();
        if (null == ver) {
            mzidVer = MzIdentMLVersion.Version_1_1;
        } else {
            mzidVer = ver;
        }
        parseFile(xfile, inputfile);
        writeMzidFile(outputfile);
    }

    public void parseFile(XTandemFile iXTandemFile, String inputfile) {

        // Setup the mzid objects
        handleCVs();

        PerformParams tandemParams = iXTandemFile.getPerformParameters();
        String version = tandemParams.getProcVersion();
        String dbLocation = tandemParams.getSequenceSource_1();
        String dbName = tandemParams.getSequenceSourceDescription_1();
        long totalProts = (long) tandemParams.getTotalProteinsUsed();

        InputParams inputParams = iXTandemFile.getInputParameters();
        //If the spectrum path is null, then we can assume all input parameters are missing
        //as the spectrum path is a mandatory for X!tandem to run:
        if (inputParams.getSpectrumPath() == null) {
            throw new RuntimeException(
                    "Expected parameter not found in X!Tandem file. Please run your X!Tandem search with option 'output, parameters=yes'. See http://thegpm.org/tandem/api/opara.html for more details.");
        }

        //TODO - This only works if the user specified to output the input params - need to document this clearly
        double massError = inputParams.getSpectrumParentMonoIsoMassErrorMinus();

        Map<String, DBSequence> foundProts = new HashMap<>();
        Map<String, PeptideEvidence> pepEvidLookup = new HashMap<>();
        //peptideLookup= new HashMap<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide>();   //lookup to get a peptide by peptideseq_varmods_fixedmods_start_stop
        uniquePeps = new HashMap<>();      //lookup to get a peptide by peptideseq_varmods_fixedmods (i.e. uniqueness check)
        sequenceCollection = new SequenceCollection();

        siList = new SpectrumIdentificationList();
        siList.setId(siiListID);

        handleAnalysisSoftware(version);
        handleAuditCollection("firstname", "secondName", "email@place.com",
                              "address", "myworkplace");
        handleProvider();                //Performed after auditcollection, since contact is needed for provider

        boolean fragmentIsMono = handleAnalysisProtocolCollection(inputParams);
        handleInputs(inputfile, dbName, dbLocation, totalProts, inputParams.
                     getSpectrumPath());
        handleAnalysisCollection(tandemParams.getProcStartTime());

        // Create fragmentation table
        /*
         * <Measure id="m_mz"> <cvParam cvRef="PSI-MS" accession="MS:1001225"
         * name="product ion m/z"/> </Measure> <Measure id="m_intensity">
         * <cvParam cvRef="PSI-MS" accession="MS:1001226" name="product ion
         * intensity"/> </Measure> <Measure id="m_error"> <cvParam
         * cvRef="PSI-MS" accession="MS:1001227" name="product ion m/z error"
         * unitAccession="MS:1000040" unitName="m/z" unitCvRef="PSI-MS"/>
         * </Measure>
         */
        Measure mzMeasure = null;
        Measure intMeasure = null;
        Measure errorMeasure = null;

        if (outputFragmentation) {
            FragmentationTable fragTable = new FragmentationTable();
            List<Measure> measureList = fragTable.getMeasure();
            mzMeasure = new Measure();
            mzMeasure.setId(measureMzID);
            List<CvParam> cvParamList = mzMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001225",
                                                     "product ion m/z", psiCV,
                                                     "MS:1000040", "m/z", psiCV));
            intMeasure = new Measure();
            intMeasure.setId(measureIntID);
            cvParamList = intMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001226",
                                                     "product ion intensity",
                                                     psiCV, "MS:1000131",
                                                     "number of counts", psiCV));
            errorMeasure = new Measure();
            errorMeasure.setId(measureErrorID);
            cvParamList = errorMeasure.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001227",
                                                     "product ion m/z error",
                                                     psiCV, "MS:1000040", "m/z",
                                                     psiCV));
            measureList.add(mzMeasure);
            measureList.add(intMeasure);
            measureList.add(errorMeasure);

            siList.setFragmentationTable(fragTable);
        }

        List<SpectrumIdentificationResult> specIdentResults = siList.
                getSpectrumIdentificationResult();
        List<Peptide> peptideList = sequenceCollection.getPeptide();
        spectraData = new SpectraData();
        spectraData.setId(spectraDataID);

        int sirCounter = 1; //Counter used to create unique ID for SpectrumIdentificationResult

        // Iterate over all the spectra
        @SuppressWarnings("unchecked")
        Iterator<Spectrum> iter = iXTandemFile.getSpectraIterator();
        while (iter.hasNext()) {
            // Get the next spectrum.
            Spectrum spectrum = iter.next();
            int spectrumNumber = spectrum.getSpectrumNumber();//note: spectrum number seems to be a sequential index. For the spectrum number as found in xtandem file use spectrumId

            /**
             * ***********************************************
             *  *** Setup SpectrumIdentificationResult ****
             * *********************************************
             */
            SpectrumIdentificationResult specIdentRes
                    = new SpectrumIdentificationResult();

            // Get the peptide hits.
            List<de.proteinms.xtandemparser.xtandem.Peptide> pepList
                    = iXTandemFile.getPeptideMap().
                    getAllPeptides(spectrumNumber);

            //int pepEvidCounter = 1;
            int siiCounter = 1; //Counter used to create unique ID for SpectrumIdentificationItem

            SpectrumIdentificationItem sii = null;
            Peptide mzidPep = null;

            List<IonType> ionTypeList = null;

            Map<String, SpectrumIdentificationItem> sIIMap
                    = new HashMap<String, SpectrumIdentificationItem>();

            for (de.proteinms.xtandemparser.xtandem.Peptide peptide : pepList) {    //This list can contain both different pep2Protein maps and 

                for (Domain domain : peptide.getDomains()) {

                    //In mzIdentML we have 1 SepctrumIdentificationItem (SII) linked to 1 Peptide item
                    //via the peptide_ref attribute. Each Peptide item is a unique combination of:
                    // peptidesequence + modifications + substitutionModifications.
                    String uniquePepKey = getPeptideKey(domain, iXTandemFile);
                    //If it is a new global peptide, then initialize a new mzIdentML Peptide object: 
                    if (!uniquePeps.containsKey(uniquePepKey)) {
                        mzidPep = new uk.ac.ebi.jmzidml.model.mzidml.Peptide();
                        mzidPep.setPeptideSequence(domain.getDomainSequence());
                        mzidPep.setId(uniquePepKey);

                        //Parse the modifications and add them to the mzidPep object:
                        parseModificationsAndSubstitutions(mzidPep, domain,
                                                           iXTandemFile,
                                                           fragmentIsMono);

                        if (!uniquePeps.containsKey(uniquePepKey)) {
                            peptideList.add(mzidPep);
                            uniquePeps.put(uniquePepKey, mzidPep);
                        }
                    } else {
                        //otherwise fetch it from the map:
                        mzidPep = uniquePeps.get(uniquePepKey);
                    }

                    //Here we have to decide whether to create a new SepctrumIdentificationItem (SII) 
                    //item OR if this new domain item is just what we call in mzIdentML a 
                    //new PeptideEvidence (i.e. the same Peptide item found in
                    //another part of the protein sequence or even in another protein).
                    String sIIKey = getSIIKey(domain, iXTandemFile);
                    if (isNewSII(sIIKey, sIIMap)) {
                        /*
                         ****************************************************
                         ****** Create new SpectrumIdentificationItem *******
                         * ***************************************************
                         */
                        sii = new SpectrumIdentificationItem();

                        //From XSD:
                        //Set to true if the producers of the file has deemed that the identification has passed a 
                        //given threshold or been validated as correct. If no such threshold has been set, value 
                        //of true should be given for all results.
                        //Since this is not something related to X!Tandem functionality (but rather expected in 
                        //a post-X!Tandem processing QC type of step), leave it as true always:
                        //TODO - check again if X!Tandem does not have any post-processing built-in which results
                        //in an output file that has items above and below a certain threshold... 
                        sii.setPassThreshold(true);

                        //add sii to sir:
                        specIdentRes.getSpectrumIdentificationItem().add(sii);
                        sIIMap.put(sIIKey, sii);

                        if (outputFragmentation) {
                            //ionTypeList = sii.getFragmentation();
                            Fragmentation frag = new Fragmentation();
                            sii.setFragmentation(frag);
                            ionTypeList = frag.getIonType();
                            parseFragmentationData(ionTypeList, domain, peptide,
                                                   iXTandemFile,
                                                   mzMeasure, intMeasure,
                                                   errorMeasure);

                            //if ionTypeList is still empty after parsing above, throw error:
                            if (ionTypeList.size() == 0) {
                                throw new RuntimeException(
                                        "Error while parsing peptide ["
                                        + peptide.getDomains().get(0).
                                        getDomainSequence()
                                        + "] identification in spectrum ["
                                        + spectrum.getSpectrumId()
                                        + "]: no fragmentation data found");
                            }
                        }

                        sii.setPeptide(mzidPep);

                        parseScoresAndOtherSIIAttributes(sii, domain, spectrum);
                        sii.setId("SII_" + sirCounter + "_" + siiCounter);
                        siiCounter++;
                    } else {
                        sii = sIIMap.get(sIIKey);
                    }

                    //Parse protein details into DBSequence objects:
                    String[] protAccessionAndDescription = parseProteinDetails(
                            foundProts, domain, peptide, iXTandemFile);
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
                        if (decoyRegex != null) {
                            if (protAccession.matches(decoyRegex)
                                    || protAccession.split(decoyRegex).length
                                    > 1) {
                                pepEvid.setIsDecoy(Boolean.TRUE);
                            }
                        }
                        //if protein description ends with ":reversed", then peptide is decoy. Example:
                        //<protein expect="-1.4" id="5736.1" uid="4657" label="tr|F1RRZ6|F1RRZ6_PIG Uncharacterized protein (Fragment) OS=Sus scrofa GN=SH3PXD2B..." sumI="3.78" >
                        //  <note label="description">tr|F1RRZ6|F1RRZ6_PIG Uncharacterized protein (Fragment) OS=Sus scrofa GN=SH3PXD2B PE=4 SV=1:reversed</note>
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

            // Get the support data for each spectrum.
            SupportData supportData = iXTandemFile.
                    getSupportData(spectrumNumber);

            // Fill the peptide map: for each spectrum get the corressponding peptide list.
            //peptideMap.put(spectrumNumber, pepList);
            /**
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
            int spectrumID;
            if (this.isMs2SpectrumIdStartingAtZero) {
                spectrumID = (spectrum.getSpectrumId());
            } else {
                //then we assume it was starting from 1:
                spectrumID = (spectrum.getSpectrumId() - 1);
            }

            //TODO add checks to find out what is the data format submitted to X!Tandem and to correct the spectrumID value accordingly
            String label = supportData.getFragIonSpectrumDescription();
            specIdentRes.setSpectrumID("index=" + spectrumID);
            specIdentRes.setSpectraData(spectraData);
            specIdentRes.setId("SIR_" + sirCounter);
            if (label != null && !label.equals("")) {
                List<CvParam> sir_cvParamList = specIdentRes.getCvParam();
                CvParam cvp = MzidLibUtils.makeCvParam("MS:1000796",
                                                       "spectrum title", psiCV,
                                                       label);
                sir_cvParamList.add(cvp);

            }

            specIdentResults.add(specIdentRes);
            sirCounter++;

            //TO - currently only implements the case where the same peptide matches to different proteins
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

    /**
     * Parses the X!Tandem details in to a new mzIdentML PeptideEvidence object
     *
     * @param foundProts
     * @param domain
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

        char pre = domain.getUpFlankSequence().charAt(domain.
                getUpFlankSequence().length() - 1);
        if (pre == '[') {
            pre = '-';
        }
        pepEvid.setPre("" + pre);    //Reports 4 chars, we need last only
        return pepEvid;
    }

    /**
     * This method returns the key to check whether the domain represents a new
     * peptide evidence, which means in mzIdentML terms: the combination of
     *
     * (peptidesequence + modifications + substitutionModifications) +
     * proteinAccession + proteinLocation(i.e. start,end)
     *
     * is unique so far. This will mean that we have to make a new mzIdentML
     * PeptideEvidence object.
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

    private String[] parseProteinDetails(Map<String, DBSequence> foundProts,
                                         Domain domain,
                                         de.proteinms.xtandemparser.xtandem.Peptide peptide,
                                         XTandemFile iXTandemFile) {
        Protein protein = iXTandemFile.getProteinMap().getProtein(domain.
                getProteinKey());

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
            dbSeq.setSearchDatabase(searchDB);
            //dbSeq.setSearchDatabase(searchDB);
            sequenceCollection.getDBSequence().add(dbSeq);
        }
        String[] result = {protAccession, protDescription};
        return result;
    }

    private void parseScoresAndOtherSIIAttributes(SpectrumIdentificationItem sii,
                                                  Domain domain,
                                                  Spectrum spectrum) {
        double evalue = domain.getDomainExpect();
        double hyperscore = domain.getDomainHyperScore();

        int precursorCharge = spectrum.getPrecursorCharge();

        //Get precursorMh reported by X!Tandem and convert it back
        //to m/z using:
        //     m/z=(mh + z*1.007276 - 1*1.007276)/z
        //where mh is M+H    and H=1.007276466812 (proton mass rounded from 1.007276466812)
        //  X!Tandem is using H=1.007276 (proton mass rounded from 1.007276466812) so we will use that
        double protonMass = 1.007276;
        double expMZ = (spectrum.getPrecursorMh() + precursorCharge * protonMass
                - protonMass) / precursorCharge;
        //round at 6 decimals:
        expMZ = Utils.round(expMZ, 6);
        sii.setExperimentalMassToCharge(expMZ);

        double calcMZ = (domain.getDomainMh() + precursorCharge * protonMass
                - protonMass) / precursorCharge;
        calcMZ = Utils.round(calcMZ, 6);
        sii.setCalculatedMassToCharge(calcMZ);
        sii.setChargeState(precursorCharge);
        sii.setRank(1);

        List<CvParam> cvParamList = sii.getCvParam();
        //<cvParam accession="MS:1001330" name="xtandem:expect" cvRef="PSI-MS"  value="1.1e-003" />
        //<cvParam accession="MS:1001331" name="xtandem:hyperscore" cvRef="PSI-MS"  value="60.4" />
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001330",
                                                 "X\\!Tandem:expect", psiCV, ""
                                                 + evalue));
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001331",
                                                 "X\\!Tandem:hyperscore", psiCV,
                                                 "" + hyperscore));
    }

    private void parseFragmentationData(List<IonType> ionTypeList, Domain domain,
                                        de.proteinms.xtandemparser.xtandem.Peptide peptide,
                                        XTandemFile iXTandemFile,
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
                //library makes indeed theoretical fragments of different charges, more specifically 
                //when the precursor ion has charge X it generates fragment ions of 
                //charges 1 to X (the rationale of the authors here must have been that fragment 
                //ions can lose charge during fragmentation, but never gain...not sure if this is 
                //always the case). 
                System.out.println(
                        "WARN: Fragment ion annotation data is inferred based on new calculations "
                        + "triggered by this (Tandem2mzid) conversion library. "
                        + "This is done because fragment ion annotation information is not present in X!Tandem file (X!Tandem SLEDGEHAMMER (2013.09.01) and previous) . "); //TODO add real logger option
                @SuppressWarnings("unchecked")
                List<FragmentIon[]> ionsForPeptide = iXTandemFile.
                        getFragmentIonsForPeptide(peptide, domain, iXTandemFile.
                                                  getInputParameters().
                                                  getSpectrumMonoIsoMassError());
                //TODO improvement for item above: a good trade-off would be that the X!Tandem parser lib reads the 
                //                     search parameters (present in the X!Tandem output file) and use in it's in silico 
                //                     fragmentation only the same ion types (i.e. b and y) used in the search.

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
                for (Entry<Integer, Map<Integer, List<FragmentIon>>> entryCharge
                        : map.entrySet()) {
                    for (Entry<Integer, List<FragmentIon>> entryType
                            : entryCharge.getValue().entrySet()) {
                        if (entryType.getValue().isEmpty()) {
                            continue;
                        }

                        IonType ion = new IonType();
                        ion.setCharge(entryCharge.getKey());
                        ion.setCvParam(getFragmentCVParam(entryType.getKey()));
                        List<FragmentArray> fragmentList = ion.
                                getFragmentArray();
                        List<Integer> ionIndexList = ion.getIndex();

                        FragmentArray mzArray = new FragmentArray();
                        FragmentArray intArray = new FragmentArray();
                        FragmentArray errorArray = new FragmentArray();
                        mzArray.setMeasure(mzMeasure);
                        intArray.setMeasure(intMeasure);
                        errorArray.setMeasure(errorMeasure);
                        List<Float> mzValues = mzArray.getValues();
                        List<Float> intValues = intArray.getValues();
                        List<Float> errorValues = errorArray.getValues();
                        fragmentList.add(mzArray);
                        fragmentList.add(intArray);
                        fragmentList.add(errorArray);

                        ionTypeList.add(ion);

                        for (FragmentIon fragIon : entryType.getValue()) {
                            //Reported MZ is the theoretical value
                            mzValues.add((float) (fragIon.getMZ() + fragIon.
                                    getTheoreticalExperimentalMassError()));
                            intValues.add((float) fragIon.getIntensity());       //Note intensity values in Tandem do not match the source spectrum, appears that some processing happens
                            errorValues.add((float) fragIon.
                                    getTheoreticalExperimentalMassError());
                            ionIndexList.add(fragIon.getNumber());  //index position
                        }
                    }
                }

            } catch (Exception e) {
                throw new RuntimeException(
                        "Error while parsing the MS/MS fragmentation data. Please check if "
                        + " your input files indeed contain fragmentation data (this is optional and depending "
                        + "on your X!Tandem settings it could be absent). Error details: ",
                        e);
            }
        }
    }

    private Map<Integer, Map<Integer, List<FragmentIon>>> orderFragmentIons(
            List<FragmentIon[]> ionsForPeptide) {
        Map<Integer, Map<Integer, List<FragmentIon>>> map = new HashMap<>();
        for (FragmentIon[] ions : ionsForPeptide) {
            for (FragmentIon ion : ions) {
                Map<Integer, List<FragmentIon>> mapCharge = map.get((int) ion.
                        getCharge());
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
        }
        return map;
    }

    private void parseModificationsAndSubstitutions(Peptide mzidPep,
                                                    Domain domain,
                                                    XTandemFile iXTandemFile,
                                                    boolean fragmentIsMono) {
        //Parse the modifications
        List<de.proteinms.xtandemparser.interfaces.Modification> fixModList
                = iXTandemFile.getModificationMap().getFixedModifications(
                        domain.getDomainKey());
        List<de.proteinms.xtandemparser.interfaces.Modification> varModList
                = iXTandemFile.getModificationMap().getVariableModifications(
                        domain.getDomainKey());

        for (de.proteinms.xtandemparser.interfaces.Modification reportedMod
                : fixModList) {

            if (!reportedMod.isSubstitution()) {
                Modification mzidMod = translateToMzidModification(reportedMod,
                                                                   domain,
                                                                   fragmentIsMono);
                mzidPep.getModification().add(mzidMod);
            } else {
                SubstitutionModification mzidSubs = translateToMzidSubstitution(
                        reportedMod, domain, fragmentIsMono);
                mzidPep.getSubstitutionModification().add(mzidSubs);
            }
        }

        for (de.proteinms.xtandemparser.interfaces.Modification reportedMod
                : varModList) {

            if (!reportedMod.isSubstitution()) {
                Modification mzidMod = translateToMzidModification(reportedMod,
                                                                   domain,
                                                                   fragmentIsMono);
                mzidPep.getModification().add(mzidMod);
            } else {
                SubstitutionModification mzidSubs = translateToMzidSubstitution(
                        reportedMod, domain, fragmentIsMono);
                mzidPep.getSubstitutionModification().add(mzidSubs);
            }
        }

    }

    private SubstitutionModification translateToMzidSubstitution(
            de.proteinms.xtandemparser.interfaces.Modification reportedMod,
            Domain domain, boolean fragmentIsMono) {
        SubstitutionModification mzidSubs = new SubstitutionModification();

        int loc = Integer.parseInt(reportedMod.getLocation());

        int pepLoc = loc - domain.getDomainStart(); //location in Tandem is given as location within the whole protein
        mzidSubs.setLocation(pepLoc + 1);

        //Note: in X!Tandem, like in mzIdentML: sequence reported is the original:
        char reportedAminoAcid = domain.getDomainSequence().charAt(pepLoc);
        mzidSubs.setOriginalResidue(reportedAminoAcid + "");

        mzidSubs.setReplacementResidue(reportedMod.getSubstitutedAminoAcid());

        //TODO - the current storing of the modification mass in either
        //setMonoisotopicMassDelta or setAvgMassDelta seems wrong as the massDelta is 
        //a delta on the parent ion mass which is always monoisotopic in xtandem. On the other hand, 
        //a modification has to have one or more matching fragments as well which all are 
        //shifted by x, x being the modification mass.
        //So check with X!Tandem developer whether the modification (and substitution) masses 
        //follow the fragment mass type or the parent mass type.
        double mass = reportedMod.getMass();

        if (fragmentIsMono) {
            mzidSubs.setMonoisotopicMassDelta(mass);
        } else {
            mzidSubs.setAvgMassDelta(mass);
        }

        return mzidSubs;

    }

    private Modification translateToMzidModification(
            de.proteinms.xtandemparser.interfaces.Modification reportedMod,
            Domain domain, boolean fragmentIsMono) {
        Modification mzidMod = new Modification();
        double mass = reportedMod.getMass();
        int loc = Integer.parseInt(reportedMod.getLocation());

        CvParam modParam = new CvParam();

        if (fragmentIsMono) {
            mzidMod.setMonoisotopicMassDelta(mass);
        } else {
            mzidMod.setAvgMassDelta(mass);
        }

        int pepLoc = loc - domain.getDomainStart(); //location in Tandem is given as location within the whole protein
        char reportedAminoAcid = domain.getDomainSequence().charAt(pepLoc);
        mzidMod.setLocation(pepLoc + 1);        //mzid starts counting from 1, except for NTerm/CTerm mods which are 0 
        List<String> residueList = mzidMod.getResidues();
        residueList.add("" + reportedAminoAcid);

        //check if we can find a known modification type that matches what is reported in the given mass and aminoAcid, within
        //a given mass error tolerance: 
        ModT unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                             fragmentIsMono, reportedAminoAcid);

        //Part below is needed because it is not clear in XTandem output what are N or C-term mods, 
        //these have the same location as mods on the first aa or last aa in the peptide.
        //So if no unimod was found above, perhaps this is because it is really a N or C-term modification. 
        //Below we try to see if it fits a known N or C-term modification:
        if (unimod == null && pepLoc == 0) {
            //See if this is a possible N-terminal mod

            //System.out.println("\tModification not found in Unimod, so look to see if it is N-terminal\n");            
            unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                            fragmentIsMono, '[');
            mzidMod.setLocation(0);
        }
        if (unimod == null && loc == domain.getDomainEnd()) {
            //See if this is a possible C-terminal mod
            //System.out.println("\tNot found, so look to see if it is N-terminal\n");
            unimod = unimodDoc.getModByMass(mass, unimodMassError,
                                            fragmentIsMono, ']');
            mzidMod.setLocation(0); //also 0?
        }

        //Set the found details to modParam. If no unimod record was found, set the modification as "unknown" 
        if (unimod != null) {
            modParam.setAccession("UNIMOD:" + unimod.getRecordId());
            modParam.setCv(unimodCV);
            modParam.setName(unimod.getTitle());
        } else {
            //modification with mass not recognized:
            modParam.setName("unknown modification");
            modParam.setCv(psiCV);
            modParam.setAccession("MS:1001460");
        }

        mzidMod.getCvParam().add(modParam);
        return mzidMod;
    }

    private boolean isNewSII(String sIIKey,
                             Map<String, SpectrumIdentificationItem> sIIMap) {
        if (sIIMap.get(sIIKey) != null) {
            return false;
        } else {
            return true;
        }
    }

    private String getPeptideKey(Domain domain, XTandemFile iXTandemFile) {
        //is really the same as in getSIIKey, but the context of both maps is different (peptide map is global and siimap is local within 
        //a SepctrumIdentificationResult ):
        return getSIIKey(domain, iXTandemFile);
    }

    private String getSIIKey(Domain domain, XTandemFile iXTandemFile) {
        List<de.proteinms.xtandemparser.interfaces.Modification> fixModList
                = iXTandemFile.getModificationMap().getFixedModifications(
                        domain.getDomainKey());
        List<de.proteinms.xtandemparser.interfaces.Modification> varModList
                = iXTandemFile.getModificationMap().getVariableModifications(
                        domain.getDomainKey());

        String fixMods = "";

        for (de.proteinms.xtandemparser.interfaces.Modification fixMod
                : fixModList) {
            String name = fixMod.getName();
            if (fixMod.isSubstitution()) {
                name += "_subs_" + fixMod.getSubstitutedAminoAcid();
            }
            int loc = Integer.parseInt(fixMod.getLocation());
            int pepLoc = loc - domain.getDomainStart() + 1;
            fixMods += name + "$" + pepLoc + ";";
        }

        String varMods = "";
        for (de.proteinms.xtandemparser.interfaces.Modification varMod
                : varModList) {
            String name = varMod.getName();
            if (varMod.isSubstitution()) {
                name += "_subs_" + varMod.getSubstitutedAminoAcid();
            }

            int loc = Integer.parseInt(varMod.getLocation());
            int pepLoc = loc - domain.getDomainStart() + 1;
            varMods += name + "$" + pepLoc + ";";
        }
        String sIIKey = domain.getDomainSequence() + "_" + varMods + "_"
                + fixMods + "_";
        return sIIKey;

    }

    public void handleCVs() {

        //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
        //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
        //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>
        cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        psiCV = new Cv();

        psiCV.setUri(
                "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        psiCV.setId(psiCvID);
//        psiCV.setVersion("2.25.0");
        psiCV.setFullName("PSI-MS");

        unimodCV = new Cv();
        unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
        unimodCV.setId(unimodID);
        unimodCV.setFullName("UNIMOD");

        unitCV = new Cv();
        unitCV.setUri(
                "http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
        unitCV.setId(unitCvID);
        unitCV.setFullName("UNIT-ONTOLOGY");

        localCvList.add(psiCV);
        localCvList.add(unimodCV);
        localCvList.add(unitCV);
    }

    public void handleAnalysisSoftware(String version) {
        analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.
                getAnalysisSoftware();
        analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("xtandem");
        // analysisSoftware.setSoftwareName(makeCvParam("MS:1001476","xtandem",psiCV));

        Param tempParam = new Param();
        tempParam.setParam(MzidLibUtils.makeCvParam("MS:1001476", "X\\!Tandem",
                                                    psiCV));
        analysisSoftware.setSoftwareName(tempParam);

        analysisSoftware.setId(analysisSoftID);
        analysisSoftware.setVersion(version);

        /*
         * TO DO - need to work out how to use Param CvParam cvParam = new
         * CvParam(); cvParam.setName("xtandem"); cvParam.setCvRef(psiCvID);
         * cvParam.setAccession("MS:1001476"); ParamAlternative paramAlt = new
         * ParamAlternative(); paramAlt.setCvParam(cvParam);
         *
         * analysisSoftware.setSoftwareName(makeCvParam("MS:1001476","xtandem",psiCV));
         * analysisSoftware.setSoftwareName(paramAlt);
         */
        analysisSoftwares.add(analysisSoftware);

    }

    public void handleProvider() {
        provider = new Provider();
        provider.setId("PROVIDER");

        ContactRole contactRole = new ContactRole();
        contactRole.setContact(docOwner);

        Role role = new Role();
        role.setCvParam(MzidLibUtils.makeCvParam("MS:1001271", "researcher",
                                                 psiCV));
        contactRole.setRole(role);

        provider.setContactRole(contactRole);

    }

    public void handleAuditCollection(String firstName, String lastName,
                                      String email, String address,
                                      String affiliationName) {
        auditCollection = new AuditCollection();
        //List<Contact> contactList = auditCollection.getContactGroup();
        List<AbstractContact> contactList = auditCollection.
                getPersonOrOrganization();
        docOwner = new Person();
        docOwner.setId("PERSON_DOC_OWNER");
        docOwner.setFirstName(firstName);
        docOwner.setLastName(lastName);
        docOwner.getCvParam().add(MzidLibUtils.makeCvParam("MS:1000587",
                                                           "contact address",
                                                           psiCV, address));

        //docOwner.setEmail(email);
        Organization org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName(affiliationName);
        org.getCvParam().add(MzidLibUtils.makeCvParam("MS:1000586",
                                                      "contact name", psiCV,
                                                      address));

        //org.setAddress(address);
        List<Affiliation> affList = docOwner.getAffiliation();
        Affiliation aff = new Affiliation();
        aff.setOrganization(org);
        affList.add(aff);
        contactList.add(docOwner);
        contactList.add(org);

    }

    public void handleAnalysisSampleCollection() {
        analysisSampleCollection = new AnalysisSampleCollection();

    }

    public void handleAnalysisCollection(String activityDate) {
        analysisCollection = new AnalysisCollection();
        List<SpectrumIdentification> specIdentList = analysisCollection.
                getSpectrumIdentification();
        SpectrumIdentification specIdent = new SpectrumIdentification();
        specIdent.setId(specIdentID);
        specIdent.setSpectrumIdentificationList(siList);
        specIdent.setSpectrumIdentificationProtocol(siProtocol);
        List<SearchDatabaseRef> searchDBRefList = specIdent.
                getSearchDatabaseRef();
        SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
        searchDBRef.setSearchDatabase(searchDB);
        searchDBRefList.add(searchDBRef);

        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
        InputSpectra inputSpec = new InputSpectra();
        inputSpec.setSpectraData(spectraData);
        inputSpecList.add(inputSpec);

        specIdentList.add(specIdent);

    }

    public boolean handleAnalysisProtocolCollection(InputParams inputParams) {

        //boolean (parentIsMono, boolean fragmentIsMono, SearchModification[] searchMods, String enzymeName, double parTolPlus, double parTolMinus, double fragTolPlus, double fragTolMinus);
        analysisProtocolCollection = new AnalysisProtocolCollection();
        List<SpectrumIdentificationProtocol> sipList
                = analysisProtocolCollection.getSpectrumIdentificationProtocol();

        siProtocol = new SpectrumIdentificationProtocol();
        siProtocol.setId(siProtocolID);
        siProtocol.setAnalysisSoftware(analysisSoftware);

        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
        //siProtocol.setSearchType(makeCvParam("MS:1001083","ms-ms search",psiCV));
        Param tempParam = new Param();
        tempParam.setParam(MzidLibUtils.
                makeCvParam("MS:1001083", "ms-ms search", psiCV));
        siProtocol.setSearchType(tempParam);

        ParamList paramList = siProtocol.getAdditionalSearchParams();
        if (paramList == null) {
            paramList = new ParamList();
            siProtocol.setAdditionalSearchParams(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();

        boolean parentIsMono = true;  //does not appear to be a way in Tandem of specifying parent mass is average
        if (parentIsMono) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001211",
                                                     "parent mass type mono",
                                                     psiCV));
        } else {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001212",
                                                     "parent mass type average",
                                                     psiCV));
        }

        boolean fragmentIsMono = true;
        if ("average".equals(inputParams.getSpectrumFragMassType())) {
            fragmentIsMono = false;
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001255",
                                                     "fragment mass type average",
                                                     psiCV));
        } else {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001256",
                                                     "fragment mass type mono",
                                                     psiCV));
        }

        // Add "no special processing" cv term if this is mzid 1.2 version
        if (mzidVer.equals(MzIdentMLVersion.Version_1_2)) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1002495",
                                                     "no special processing",
                                                     psiCV));
        }

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams.
                getSearchModification();

        //residue, potential modification mass
        String varMods = inputParams.getResiduePotModMass();
        String fixedMods = inputParams.getResidueModMass();

        if (varMods != null) {
            String[] allVarMods = varMods.split(",");

            //	<note type="input" label="residue, modification mass">57.021469@C</note>
            //      <note type="input" label="residue, potential modification mass">15.99492@M</note>
            for (String varMod : allVarMods) {

                //If there are no modification parameters, this is identified with the value 'None'. 
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
                //If there are no modification parameters, this is identified with the value 'None'. 
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

 /*
         * <Enzymes independent="0"> <Enzyme id="ENZ_1" CTermGain="OH"
         * NTermGain="H" missedCleavages="1" semiSpecific="0"> <EnzymeName>
         * <cvParam accession="MS:1001251" name="Trypsin" cvRef="PSI-MS" />
         * </EnzymeName> </Enzyme> </Enzymes>
         */
        //Only add this group if there are any modifications in the list:
        if (searchModList.size() > 0) {
            siProtocol.setModificationParams(modParams);
        }
        Enzymes enzymes = siProtocol.getEnzymes();

        if (enzymes == null) {
            enzymes = new Enzymes();
            siProtocol.setEnzymes(enzymes);
        }

        enzymes.setIndependent(false);

        List<Enzyme> enzymeList = enzymes.getEnzyme();

        Enzyme enzyme = getEnzyme(inputParams.getProteinCleavageSite(),
                                  inputParams.getScoringMissCleavageSites());

        enzymeList.add(enzyme);

        Tolerance fragTol = new Tolerance();
        Tolerance parTol = new Tolerance();

        boolean isDaltons;

        if (inputParams.getSpectrumMonoIsoMassErrorUnits().equalsIgnoreCase(
                "ppm")) {
            isDaltons = false;
        } else {
            isDaltons = true;
            unimodMassError = inputParams.getSpectrumMonoIsoMassError();        //Dynamically set Unimod lookup mass error
        }

        List<CvParam> fragCvList = fragTol.getCvParam();
        CvParam fragCvPlus = getCvParamWithMassUnits(isDaltons);
        CvParam fragCvMinus = getCvParamWithMassUnits(isDaltons);


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
        fragCvPlus.setValue("" + inputParams.getSpectrumMonoIsoMassError());
        fragCvMinus.setValue("" + inputParams.getSpectrumMonoIsoMassError());
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        List<CvParam> parCvList = parTol.getCvParam();

        if (inputParams.getSpectrumParentMonoIsoMassErrorUnits().
                equalsIgnoreCase("ppm")) {
            isDaltons = false;
        } else {
            isDaltons = true;
        }

        CvParam parCvPlus = getCvParamWithMassUnits(isDaltons);
        CvParam parCvMinus = getCvParamWithMassUnits(isDaltons);

        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue("" + inputParams.
                getSpectrumParentMonoIsoMassErrorPlus());
        parCvMinus.setValue("" + inputParams.
                getSpectrumParentMonoIsoMassErrorMinus());
        parCvList.add(parCvPlus);
        parCvList.add(parCvMinus);

        siProtocol.setFragmentTolerance(fragTol);
        siProtocol.setParentTolerance(parTol);

        ParamList sip_paramList = siProtocol.getThreshold();
        if (sip_paramList == null) {
            sip_paramList = new ParamList();
            siProtocol.setThreshold(sip_paramList);
        }
        cvParamList = sip_paramList.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001494", "no threshold",
                                                 psiCV));
        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />
        sipList.add(siProtocol);

        return fragmentIsMono;
    }

    private SearchModification translateToSearchModification(String reportedMod,
                                                             boolean fragmentIsMono,
                                                             boolean isFixedMod) {
        List<String> residues = new ArrayList<String>();
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
                                                                psiCV));
        } else {

            searchMod.getCvParam().add(MzidLibUtils.makeCvParam("UNIMOD:"
                    + unimod.getRecordId(), unimod.getTitle(), unimodCV));
        }
        searchMod.setMassDelta(new Float(monoMass));

        for (String residue : residues) {
            if (residue.equals("[")) {
                searchMod.getCvParam().add(MzidLibUtils.
                        makeCvParam("MS:1001189",
                                    "modification specificity N-term", psiCV));
                residue = ".";
            } else if (residue.equals("]")) {
                searchMod.getCvParam().add(MzidLibUtils.
                        makeCvParam("MS:1001190",
                                    "modification specificity C-term", psiCV));
                residue = ".";     //The any char must be inserted into mzid
            }
            searchMod.getResidues().add(residue);
        }
        return searchMod;
    }

    public void handleInputs(String tandemFileLocation, String databaseName,
                             String databaselocation, long numProts,
                             String inputSpectraLocation) {

        inputs = new Inputs();
        List<SearchDatabase> searchDBList = inputs.getSearchDatabase();

        searchDB = new SearchDatabase();
        searchDB.setId(searchDBID);
        searchDB.setNumDatabaseSequences(numProts);
        //<cvParam accession="MS:1001401" name="xtandem xml file" cvRef="PSI-MS"/>

        UserParam param = new UserParam();
        param.setName(databaseName);
        Param tempParam = new Param();
        tempParam.setParam(param);
        searchDB.setDatabaseName(tempParam);

        //searchDB.setDatabaseName(param);
        searchDB.setLocation(databaselocation);

        FileFormat ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam(this.databaseFileFormatID,
                                               this.databaseFileFormatName,
                                               psiCV));
        searchDB.setFileFormat(ff);
        searchDBList.add(searchDB);

        List<SourceFile> sourceFileList = inputs.getSourceFile();
        SourceFile sourceFile = new SourceFile();
        sourceFile.setLocation(tandemFileLocation);
        sourceFile.setId(sourceFileID);
        ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam("MS:1001401",
                                               "X\\!Tandem xml file", psiCV));
        sourceFile.setFileFormat(ff);
        sourceFileList.add(sourceFile);

        List<SpectraData> spectraDataList = inputs.getSpectraData();
        spectraData = new SpectraData();

        SpectrumIDFormat sif = new SpectrumIDFormat();
        sif.setCvParam(MzidLibUtils.makeCvParam("MS:1000774",
                                                "multiple peak list nativeID format",
                                                psiCV));
        spectraData.setSpectrumIDFormat(sif);
        //spectraData.setSpectrumIDFormat(makeCvParam("MS:1000774","multiple peak list nativeID format",psiCV));

        ff = new FileFormat();
        ff.setCvParam(MzidLibUtils.makeCvParam(this.massSpecFileFormatID,
                                               this.massSpecFileFormatName,
                                               psiCV));
        spectraData.setFileFormat(ff);
        spectraData.setId(spectraDataID);
        spectraData.setLocation(inputSpectraLocation);
        spectraDataList.add(spectraData);

        /*
         * <fileFormat> <cvParam accession="MS:1001062" name="Mascot MGF file"
         * cvRef="PSI-MS" /> </fileFormat> <spectrumIDFormat> <cvParam
         * accession="MS:1000774" name="multiple peak list nativeID format"
         * cvRef="PSI-MS" /> </spectrumIDFormat>
         */
        //From tandem: <note label="modelling, total proteins used">22348</note>
        /*
         * <SearchDatabase
         * location="file:///C:/inetpub/mascot/sequence/SwissProt/current/SwissProt_51.6.fasta"
         * id="SDB_SwissProt" name="SwissProt" numDatabaseSequences="257964"
         * numResidues="93947433" releaseDate="SwissProt_51.6.fasta"
         * version="SwissProt_51.6.fasta"> <fileFormat> <cvParam
         * accession="MS:1001348" name="FASTA format" cvRef="PSI-MS" />
         * </fileFormat> <DatabaseName> <userParam name="SwissProt_51.6.fasta"
         * /> </DatabaseName> <cvParam accession="MS:1001073" name="database
         * type amino acid" cvRef="PSI-MS" /> </SearchDatabase>
         */
    }

    public void writeMzidFile(String outputfile) {

        try {
            Writer writer = null;
            if (outputfile.equals("")) {
                writer = new FileWriter("55merge_tandem.mzid");
            } else {
                writer = new FileWriter(outputfile);
            }

            MzIdentMLMarshaller m = new MzIdentMLMarshaller();

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
            writer.write(m.createXmlHeader() + "\n");

            // mzIdentML start tag
            writer.write(m.createMzIdentMLStartTag("12345") + "\n");

            m.marshall(cvList, writer);
            writer.write("\n");
            AnalysisSoftware analysisSoftware = new AnalysisSoftware();
            Date date = new Date();
            SimpleDateFormat dateFormat = new SimpleDateFormat(
                    "yyyy-MM-dd HH-mm-ss");
            //TODO: this seems strange...a date in the id and name fields? Check if there
            //is perhaps a date field where the date information could be stored.
            analysisSoftware.setName(this.getClass().getSimpleName() + "_"
                    + dateFormat.format(date));
            analysisSoftware.setId(this.getClass().getSimpleName() + "_"
                    + dateFormat.format(date));
            Param param = new Param();
            param.setParam(MzidLibUtils.makeCvParam("MS:1002237", "mzidLib",
                                                    psiCV));
            analysisSoftware.setSoftwareName(param);
            analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);
            m.marshal(analysisSoftwareList, writer);
            writer.write("\n");

            m.marshall(provider, writer);
            writer.write("\n");

            m.marshall(auditCollection, writer);
            writer.write("\n");

            //m.marshall(analysisSampleCollection, writer);     //TODO - complete this part
            //writer.write("\n");
            m.marshall(sequenceCollection, writer);
            writer.write("\n");

            m.marshall(analysisCollection, writer);
            writer.write("\n");

            m.marshall(analysisProtocolCollection, writer);
            writer.write("\n");

            writer.write(m.createDataCollectionStartTag() + "\n");
            m.marshall(inputs, writer);
            writer.write("\n");

            //Inputs inputs = unmarshaller.unmarshal(MzIdentMLElement.Inputs.getXpath());
            //m.marshall(inputs, writer);
            //writer.write("\n");
            writer.write(m.createAnalysisDataStartTag() + "\n");

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
            m.marshall(siList, writer);
            writer.write("\n");

            // writer.write(m.createSpectrumIdentificationListClosingTag() + "\n");
            writer.write(m.createProteinDetectionListStartTag("PDL_1", null)
                    + "\n");

            /*
             * Iterator<ProteinAmbiguityGroup> protAmbGroupIter =
             * unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
             * while (protAmbGroupIter.hasNext()) { ProteinAmbiguityGroup
             * protAmbGroup = protAmbGroupIter.next(); m.marshall(protAmbGroup,
             * writer); writer.write("\n"); }
             *
             */
            writer.write(m.createProteinDetectionListClosingTag() + "\n");

            writer.write(m.createAnalysisDataClosingTag() + "\n");

            writer.write(m.createDataCollectionClosingTag() + "\n");

            //BibliographicReference ref = unmarshaller.unmarshal(MzIdentMLElement.BibliographicReference.getXpath());
            // m.marshall(ref, writer);
            // writer.write("\n");
            writer.write(m.createMzIdentMLClosingTag());

            writer.close();
        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + ex.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    public CvParam getCvParamWithMassUnits(boolean isDaltonUnit) {
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

    public CvParam getFragmentCVParam(int iontype) {

        CvParam cvParam = new CvParam();
        cvParam.setCv(psiCV);

        //nb: don't forget the "break" when adding new ones here!
        switch (iontype) {
            case FragmentIon.A_ION:
                cvParam.setAccession("MS:1001229");
                cvParam.setName("frag: a ion");
                break;
            case FragmentIon.AH2O_ION:
                cvParam.setAccession("MS:1001234");
                cvParam.setName("frag: a ion - H2O");
                break;
            case FragmentIon.ANH3_ION:
                cvParam.setAccession("MS:1001235");
                cvParam.setName("frag: a ion - NH3");
                break;
            case FragmentIon.B_ION:
                cvParam.setAccession("MS:1001224");
                cvParam.setName("frag: b ion");
                break;
            case FragmentIon.BH2O_ION:
                cvParam.setAccession("MS:1001222");
                cvParam.setName("frag: b ion - H2O");
                break;
            case FragmentIon.BNH3_ION:
                cvParam.setAccession("MS:1001232");
                cvParam.setName("frag: b ion - NH3");
                break;
            case FragmentIon.C_ION:
                cvParam.setAccession("MS:1001231");
                cvParam.setName("frag: c ion");
                break;
            case FragmentIon.X_ION:
                cvParam.setAccession("MS:1001228");
                cvParam.setName("frag: x ion");
                break;
            case FragmentIon.Y_ION:
                cvParam.setAccession("MS:1001220");
                cvParam.setName("frag: y ion");
                break;
            case FragmentIon.YH2O_ION:
                cvParam.setAccession("MS:1001223");
                cvParam.setName("frag: y ion - H2O");
                break;
            case FragmentIon.YNH3_ION:
                cvParam.setAccession("MS:1001233");
                cvParam.setName("frag: y ion - NH3");
                break;
            case FragmentIon.Z_ION:
                cvParam.setAccession("MS:1001230");
                cvParam.setName("frag: z ion");
                break;
            case FragmentIon.MH_ION:
                cvParam.setAccession("MS:1001523");
                cvParam.setName("frag: precursor ion");
                break;
            case FragmentIon.MHNH3_ION:
                cvParam.setAccession("MS:1001522");
                cvParam.setName("frag: precursor ion - NH3");
                break;
            case FragmentIon.MHH2O_ION:
                cvParam.setAccession("MS:1001521");
                cvParam.setName("frag: precursor ion - H2O");
                break;
        }

        return cvParam;

    }

    public Enzyme getEnzyme(String tandemRegex, int missedCleavage) {

        Enzyme enzyme = new Enzyme();
        //[KR]|{P}

        //TODO only trypsin implemented - this is difficult to convert from Regex used in X!Tandem
        enzyme.setId("Enz1");
        enzyme.setCTermGain("OH");
        enzyme.setNTermGain("H");
        enzyme.setMissedCleavages(missedCleavage);
        enzyme.setSemiSpecific(false);

        if (tandemRegex.equalsIgnoreCase("[KR]|{P}")) {

            ParamList paramList = enzyme.getEnzymeName();
            if (paramList == null) {
                paramList = new ParamList();
                enzyme.setEnzymeName(paramList);
            }
            List<CvParam> cvParamList = paramList.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001251", "Trypsin",
                                                     psiCV));
        } else {
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

}
