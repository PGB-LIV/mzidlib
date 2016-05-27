package uk.ac.liv.mzidlib.converters;

import uk.ac.liv.unimod.ModT;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.lang.Math;

import de.proteinms.omxparser.*;
import de.proteinms.omxparser.util.*;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 *
 * @author jonesar
 */
public class Omssa2mzid {

    /*
     *
     * OmssaOmxFile omxFile = new OmssaOmxFile(â€œC:\\OMSSA_Files\\BSA.omxâ€�);
     *
     * HashMap<MSSpectrum, MSHitSet> results = omxFile.getSpectrumToHitSetMap();
     * Iterator<MSSpectrum> iterator = results.keySet().iterator();
     *
     * ArrayList<List<Integer>> allMzValues = new ArrayList();
     *
     * while (iterator.hasNext()) { MSSpectrum tempSpectrum = iterator.next();
     * allMzValues.add(tempSpectrum.MSSpectrum_mz.MSSpectrum_mz_E); }
     *
     */
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
    String inputOmssaFile = "example_files/55merge_omssa.omx";
    String modsFile = "build/classes/resources/mods.xml";
    //private String modsFile = getClass().getClassLoader().getResource("lib/mods.xml").getPath();
    //URL modsFileURL = getClass().getClassLoader().getResource("resources/mods.xml");
    //URL userModsFileURL = getClass().getClassLoader().getResource("resources/mods.xml");
    static String userModsFile = "build/classes/resources/usermods.xml";
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
    static String decoyRegex = null;    //Added by ARJ for setting is decoy
    //Some objects we will need globally
    Cv unimodCV;
    Cv psiCV;
    Cv unitCV;
    SpectrumIdentificationProtocol siProtocol;
    SearchDatabase searchDB;
    SpectraData spectraData;
    Person docOwner;
    //PeptideEvidenceList pepEvidList;
    AnalysisSoftware analysisSoftware;
    Map<String, DBSequence> foundProts;
    Map<String, String> pepProtMap;
    Map<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide> peptideLookup;   //lookup to get a peptide by peptideseq_varmods_fixedmods
    Map<String, PeptideEvidence> pepEvidLookup;                       //lookup to get a peptide evidence object by peptideID_proteinacc_start_end
    Map<String, Boolean> uniquePeps;
    Map<Integer, OmssaModification> intToModMap;                         //Mapping from Omssa Mod integers to objects
    int pepCounter = 0;
    int pepEvidCounter = 0;
    MSSearchSettings settings = null;
    //List<SpectrumIdentificationResult> specIdentResults;
    int sirCounter = 1; //Counter used to create unique ID for SpectrumIdentificationResult
    ReadUnimod unimodDoc;
    double unimodMassError = 0.1;
    Boolean outputFragmentation = false;
    int response_scale = 100;           // This is the scale to get correct MZ values out - get reset from the omx
    private String defline_regex = " "; //TODO  - current grabs protein accessions from the defline, by space - need to implement other options

    public Omssa2mzid(String inputfile) {

        //System.out.println("modsFileURL:" + modsFileURL);
        //OmssaOmxFile omxFile = new OmssaOmxFile(fileName,modsFileURL.getPath().replaceAll("file:/",""),userModsFileURL.getPath().replaceAll("file:/",""));
        OmssaOmxFile omxFile = new OmssaOmxFile(inputfile, modsFile, userModsFile);
        unimodDoc = new ReadUnimod();
        parseFile(omxFile);
        writeMzidFile("");

    }
    /*
     * By F. Ghali - ARJ: I have removed this constructor. Use other one with a
     * null value for decoyRegularExpression
     *
     * public Omssa2mzid(String inputfile, String outputfile, Boolean
     * outputFrags){
     *
     * try{
     *
     * //for this constructor - we will extract files back out of the jar
     * modsFile = "mods.xml"; InputStream inMods =
     * this.getClass().getResourceAsStream("/resources/mods.xml"); userModsFile
     * = "usermods.xml"; InputStream inUserMods =
     * this.getClass().getResourceAsStream("/resources/usermods.xml");
     *
     * extractFileFromJar(inMods,modsFile);
     * extractFileFromJar(inUserMods,userModsFile);
     *
     * //OmssaOmxFile omxFile = new
     * OmssaOmxFile(inputfile,modsFileURL.getPath().replaceAll("file:/",""),userModsFileURL.getPath().replaceAll("file:/",""));
     * OmssaOmxFile omxFile = new OmssaOmxFile(inputfile,modsFile,userModsFile);
     * unimodDoc = new ReadUnimod(); this.outputFragmentation = outputFrags;
     * parseFile(omxFile); writeMzidFile(outputfile);
     *
     * File modFile = new File(modsFile); File userModFile = new
     * File(userModsFile); modFile.delete(); userModFile.delete();
     *
     * }
     * catch(Exception e){ e.printStackTrace(); } }
     */

    // By F. Ghali
    public Omssa2mzid(String inputfile, String outputfile, Boolean outputFrags, String decoyRegularExpression, String omssaModsFile, String omssaUserModsFile) {

        //for this constructor - we will extract files back out of the jar
        modsFile = omssaModsFile;
        userModsFile = omssaUserModsFile;

        boolean cleanupOmssaMods = false;
        boolean cleanupUserMods = false;

        if (omssaModsFile == null || omssaModsFile.equals("")) {
            System.out.println("Using the default mods.xml file\n");
            modsFile = "mods.xml";
            InputStream inMods = ClassLoader.getSystemClassLoader().getResourceAsStream("mods.xml");
            extractFileFromJar(inMods, modsFile);
            cleanupOmssaMods = true;
        }

        if (userModsFile == null || userModsFile.equals("")) {
            System.out.println("Using the default usermods.xml file\n");
            userModsFile = "usermods.xml";
            InputStream inUserMods = ClassLoader.getSystemClassLoader().getResourceAsStream("usermods.xml");
            extractFileFromJar(inUserMods, userModsFile);
            cleanupUserMods = true;
        }

        //OmssaOmxFile omxFile = new OmssaOmxFile(inputfile,modsFileURL.getPath().replaceAll("file:/",""),userModsFileURL.getPath().replaceAll("file:/",""));
        OmssaOmxFile omxFile = new OmssaOmxFile(inputfile, modsFile, userModsFile);

        if (decoyRegularExpression != null) {
            decoyRegex = decoyRegularExpression;
        }
        outputFragmentation = outputFrags;
        unimodDoc = new ReadUnimod();
        parseFile(omxFile);
        writeMzidFile(outputfile);

        if (cleanupOmssaMods) {
            File modFile = new File(modsFile);
            modFile.delete();
        }
        if (cleanupUserMods) {
            File userModFile = new File(userModsFile);
            userModFile.delete();
        }

    }

    /*
     * Helper method to get the mod file out of the jar, since local file is
     * needed for OMXParser
     */
    private void extractFileFromJar(InputStream in, String filename) {
        try {

            StringBuilder builder = new StringBuilder();
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                builder.append(line);
            }
            br.close();
            String text = builder.toString();
            FileWriter writer = new FileWriter(filename);
            writer.write(text);
            writer.close();

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    public void parseFile(OmssaOmxFile omxFile) {

        double maxEvalue = 0.0;

        Map<MSSpectrum, MSHitSet> results = omxFile.getSpectrumToHitSetMap();
        Iterator<MSSpectrum> iterator = results.keySet().iterator();

        List<List<Integer>> allMzValues = new ArrayList<>();

        while (iterator.hasNext()) {
            MSSpectrum tempSpectrum = iterator.next();
            allMzValues.add(tempSpectrum.MSSpectrum_mz.MSSpectrum_mz_E);
        }

        // Iterate over all the spectra
        Iterator<MSSpectrum> iter = results.keySet().iterator();
        //Iterator<Spectrum> iter = iXTandemFile.getSpectraIterator();

        Map<MSSpectrum, HashSet<String>> specToPepMap = omxFile.getSpectrumToPeptideMap();
        MSSearch mssearch = omxFile.getParserResult();
        MSSearch_response responses = mssearch.MSSearch_response;
        MSSearch_request requests = mssearch.MSSearch_request;

        int reqCounter = 0;

        for (MSRequest request : requests.MSRequest) {

            MSRequest_settings req_settings = request.MSRequest_settings;
            settings = req_settings.MSSearchSettings;
            reqCounter++;

            if (reqCounter > 1) {
                System.out.println("Error: multiple requests in the OMX file - this is not currently supported");
            }
        }

        MSResponse response = null;
        int respCounter = 0;
        for (MSResponse resp : responses.MSResponse) {
            response = resp;
            response_scale = response.MSResponse_scale;

            respCounter++;

            if (respCounter > 1) {
                System.out.println("Error: multiple responses in the OMX file - this is not currently supported");
            }

            // String version = response.MSResponse_version;
            //TO DO get other running params etc.
            //String bioseq = response.MSResponse_bioseqs;
            //System.out.println("Version: " + version + " bioseqs: " + bioseq);
        }

        Map<MSSpectrum, MSHitSet> specToHitMap = omxFile.getSpectrumToHitSetMap();
        intToModMap = omxFile.getModifications();
        //HashMap<String, LinkedList<MSPepHit>> pepToProtMap = omxFile.getPeptideToProteinMap();

        // Setup the mzid objects
        handleCVs();

        foundProts = new HashMap<String, DBSequence>();
        pepProtMap = new HashMap<String, String>();
        peptideLookup = new HashMap<String, uk.ac.ebi.jmzidml.model.mzidml.Peptide>();   //lookup to get a peptide by peptideseq_varmods_fixedmods
        pepEvidLookup = new HashMap<String, PeptideEvidence>();
        uniquePeps = new HashMap<String, Boolean>();
        sequenceCollection = new SequenceCollection();

        List<PeptideEvidence> peptideEvidenceList = sequenceCollection.getPeptideEvidence();

        siList = new SpectrumIdentificationList();
        siList.setId(siiListID);

        handleAnalysisSoftware("");      //TODO - not clear that we can get the software version from the results file

        handleAuditCollection("firstname", "secondName", "email@place.com", "address", "myworkplace");
        handleProvider();                //Performed after auditcollection, since contact is needed for provider

        handleAnalysisProtocolCollection();  //This method is to fix TODO

        handleInputs(response, settings);
        handleAnalysisCollection("TIME:TODO");

        // Create fragmentation table
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
        FragmentationTable fragTable = new FragmentationTable();
        List<Measure> measureList = fragTable.getMeasure();
        Measure mzMeasure = new Measure();
        mzMeasure.setId(measureMzID);
        List<CvParam> cvParamList = mzMeasure.getCvParam();
        cvParamList.add(makeCvParam("MS:1001225", "product ion m/z", psiCV, "MS:1000040", "m/z", psiCV));
        Measure intMeasure = new Measure();
        intMeasure.setId(measureIntID);
        cvParamList = intMeasure.getCvParam();
        cvParamList.add(makeCvParam("MS:1001226", "product ion intensity", psiCV, "MS:1000131", "number of counts", psiCV));
        Measure errorMeasure = new Measure();
        errorMeasure.setId(measureErrorID);
        cvParamList = errorMeasure.getCvParam();
        cvParamList.add(makeCvParam("MS:1001227", "product ion m/z error", psiCV, "MS:1000040", "m/z", psiCV));

        measureList.add(mzMeasure);
        measureList.add(intMeasure);
        measureList.add(errorMeasure);

        siList.setFragmentationTable(fragTable);

        List<SpectrumIdentificationResult> specIdentResults = siList.getSpectrumIdentificationResult();
        List<DBSequence> dbSeqList = sequenceCollection.getDBSequence();
        List<Peptide> peptideList = sequenceCollection.getPeptide();
        spectraData = new SpectraData();
        spectraData.setId(spectraDataID);

        while (iter.hasNext()) {

            // Get the next spectrum.
            MSSpectrum spectrum = iter.next();
            //int spectrumNumber = spectrum.getSpectrumNumber();
            int spectrumNumber = spectrum.MSSpectrum_number;

            Set<String> peptides = specToPepMap.get(spectrum);

            //MSSpectrum_ids msspec_ids = spectrum.MSSpectrum_ids;
            Iterator<String> pepIter = peptides.iterator();

            //System.out.println("Spectrum " + spectrumNumber);

            /*
             * while (pepIter.hasNext()) {
             *
             * String pepSeq = (String)pepIter.next(); System.out.println("\t" +
             * pepSeq); }
             */
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
                MSHitSet_error ms_error = msHitSet.MSHitSet_error;
                //int error = ms_error.MSHitError;
                //System.out.println("Error:" + ms_error.MSHitError);                
                specIdentRes.setSpectraData(spectraData);
                specIdentRes.setSpectrumID("index=" + spectrum.MSSpectrum_number);

                MSSpectrum_ids specIDs = spectrum.MSSpectrum_ids;

                if (!specIDs.MSSpectrum_ids_E.isEmpty()) {
                    List<CvParam> sir_cvParamList = specIdentRes.getCvParam();
                    for (String id : specIDs.MSSpectrum_ids_E) {
                        CvParam cvp = makeCvParam("MS:1000796", "spectrum title", psiCV, id);
                        sir_cvParamList.add(cvp);
                    }
                }

                int rank = 1;

                int siiCounter = 1; //Counter used to create unique ID for SpectrumIdentificationItem

                for (MSHits hits : msHits.MSHits) {

                    SpectrumIdentificationItem sii = new SpectrumIdentificationItem();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = sii.getPeptideEvidenceRef();
                    siiList.add(sii);

                    int charge = hits.MSHits_charge;
                    double evalue = hits.MSHits_evalue;
                    double pvalue = hits.MSHits_pvalue;
                    String pre = hits.MSHits_pepstart;
                    String post = hits.MSHits_pepstop;
                    String pepSeq = hits.MSHits_pepstring;

                    sii.setId("SII_" + sirCounter + "_" + siiCounter);
                    siiCounter++;

                    sii.setChargeState(charge);
                    sii.setRank(rank);
                    sii.setPassThreshold(true);              //by default these are supposed to be set to true
                    rank++;
                    //sii.setExperimentalMassToCharge(hits.MSHits_mass);

                    long expMZ = java.lang.Math.round(hits.MSHits_mass / hits.MSHits_charge);
                    long theoMass = java.lang.Math.round(hits.MSHits_theomass / hits.MSHits_charge);

                    //System.out.println("expMZ:" + expMZ + "hits.MSHits_mass " + hits.MSHits_mass);
                    if (evalue > maxEvalue) {
                        maxEvalue = evalue;
                    }

                    sii.setExperimentalMassToCharge((double) expMZ / 1000);
                    sii.setCalculatedMassToCharge((double) theoMass / 1000);

                    cvParamList = sii.getCvParam();
                    cvParamList.add(makeCvParam("MS:1001328", "OMSSA:evalue", psiCV, "" + evalue));
                    cvParamList.add(makeCvParam("MS:1001329", "OMSSA:pvalue", psiCV, "" + pvalue));

                    MSHits_mods mods = hits.MSHits_mods;

                    Peptide mzidPep = new Peptide();
                    mzidPep.setPeptideSequence(pepSeq);

                    String modString = convertVarMods(mzidPep, mods);
                    convertFixedMods(mzidPep);                       //Uses the fixed mods in the search settings to guess the fixed modifications

                    List<IonType> ionTypeList = null;
                    if (outputFragmentation) {
                        MSHits_mzhits mzhits = hits.MSHits_mzhits;
                        Map<String, List<MSMZHit>> mapIonsToHits = new HashMap<>();

                        //ionTypeList = sii.getFragmentation();
                        Fragmentation frag = new Fragmentation();
                        sii.setFragmentation(frag);
                        ionTypeList = frag.getIonType();

                        for (MSMZHit mzhit : mzhits.MSMZHit) {

                            int index = mzhit.MSMZHit_index;
                            int ionCharge = mzhit.MSMZHit_charge;
                            double ionMz = (double) mzhit.MSMZHit_mz;
                            int number = mzhit.MSMZHit_number;

                            MSMZHit_annotation ionAnnot = mzhit.MSMZHit_annotation;

                            MSIonAnnot ionA = ionAnnot.MSIonAnnot;
                            double massError = ionA.MSIonAnnot_massdiff;

                            MSMZHit_ion ion = mzhit.MSMZHit_ion;
                            int ionType = ion.MSIonType;
                            String ionString = OmssaEnumerators.getMSIonTypeAsText(ionType);
                            //System.out.println("\t" + index + "\t" + number + "\t" + ionCharge + "\t" + ionString + "\t" + ionMz + "\t" + massError);

                            //Logic is slightly different from Tandem parser, since all ion types appear to be mixed up in one structure
                            //1. Work out all the ion types we have, put in a temporary HashMap?
                            String testString = ionType + "_" + charge;

                            List<MSMZHit> mzHitList = null;
                            if (mapIonsToHits.containsKey(testString)) {
                                mzHitList = mapIonsToHits.get(testString);
                            } else {
                                mzHitList = new ArrayList<MSMZHit>();
                                mapIonsToHits.put(testString, mzHitList);
                            }
                            mzHitList.add(mzhit);

                        }

                        for (String ionKey : mapIonsToHits.keySet()) {

                            IonType mzidIon = new IonType();
                            List<Integer> ionIndexList = mzidIon.getIndex();
                            List<FragmentArray> fragmentList = mzidIon.getFragmentArray();

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
                                String ionString = OmssaEnumerators.getMSIonTypeAsText(ionType);

                                if (j == 0) {
                                    int ionCharge = mzhit.MSMZHit_charge;
                                    mzidIon.setCharge(ionCharge);
                                    //System.out.println("Lookup:" + ionType);

                                    CvParam cvParam = getFragmentCVParam(ionType);
                                    mzidIon.setCvParam(cvParam);
                                }
                                j++;

                                double theoMZ = (double) (mzhit.MSMZHit_mz) / response_scale;
                                double[] matchedPeak = getMatchedIon(spectrum, theoMZ);
                                mzValues.add((float) matchedPeak[0]);
                                intValues.add((float) matchedPeak[1]);
                                errorValues.add((float) (matchedPeak[0] - theoMZ));
                                ionIndexList.add(mzhit.MSMZHit_number);  //index position

                            }
                            fragmentList.add(mzArray);
                            fragmentList.add(intArray);
                            fragmentList.add(errorArray);

                            ionTypeList.add(mzidIon);
                        }
                    }

                    String uniquePep = pepSeq + modString;
                    if (!uniquePeps.containsKey(uniquePep)) {
                        peptideList.add(mzidPep);
                    }
                    uniquePeps.put(uniquePep, true);

                    //System.out.println("\t" + charge + " " + evalue + " " + pre  + " " + post + " " +pepSeq);

                    /*
                     *********************************************************
                     ******** These are peptide evidences ****************
                     * ********************************************************
                     */
                    MSHits_pephits pephits = hits.MSHits_pephits;
                    for (MSPepHit pephit : pephits.MSPepHit) {
                        int start = pephit.MSPepHit_start + 1; //start and end positions are mapped from zero in Omssa but from one in the rest of the known universe
                        int end = pephit.MSPepHit_stop + 1;   //Note this is not an issue in the CSV version

                        String defline = pephit.MSPepHit_defline;
                        String protObjID = pephit.MSPepHit_accession;
                        // Added By FG
                        String protAcc;
                        if (defline.indexOf(defline_regex) >= 0) {
                            protAcc = defline.substring(0, defline.indexOf(defline_regex));
                        } else {
                            protAcc = defline;
                        }

                        String protSeq = "";
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
                            dbSeq.setSearchDatabase(searchDB);
                            List<CvParam> db_cvParamList = dbSeq.getCvParam();
                            String desc = defline;

                            db_cvParamList.add(makeCvParam("MS:1001088", "protein description", psiCV, desc));
                            dbSeqList.add(dbSeq);
                        } else {
                            dbSeq = foundProts.get(protAcc);
                        }

                        String testPepMods = pepSeq + "_" + modString + "_" + start + "_" + end;
                        String testProt = protAcc + "_" + start + "_" + end;

                        //Check this is a unique peptide
                        //TODO - Wasteful of memory, since mzidPep doesn't need to be created if this is not a unique peptide
                        boolean newPepEvid = false;
                        String pepID = null;

                        if (!pepProtMap.containsKey(testPepMods)) {

                            //Add peptide sequence to mzid Peptide object
                            pepID = "Peptide" + pepCounter;
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
                            if (mappedProts.indexOf(testProt) == -1) {  //test if testProt is within mapped prots - if yes, do nothing, if not, create a new peptideevidence
                                pepProtMap.put(testPepMods, mappedProts + ";" + testProt);
                                newPepEvid = true;
                            }
                            mzidPep = peptideLookup.get(testPepMods);
                        }

                        PeptideEvidence pepEvid = null;
                        if (newPepEvid) {
                            pepEvid = new PeptideEvidence();
                            pepEvidLookup.put(mzidPep.getId() + "_" + testProt, pepEvid);
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
                            pepEvid.setId("PE" + sirCounter + "_" + siiCounter + "_" + pepEvidCounter);
                            pepEvidCounter++;
                            pepEvid.setIsDecoy(Boolean.FALSE);
                            if (decoyRegex != null) {
                                if (protAcc.contains(decoyRegex)) {
                                    pepEvid.setIsDecoy(Boolean.TRUE);
                                }
                            }
                            peptideEvidenceList.add(pepEvid);
                        } else {
                            pepEvid = pepEvidLookup.get(mzidPep.getId() + "_" + testProt);

                        }

                        PeptideEvidenceRef peptideEvidenceRef = new PeptideEvidenceRef();
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
                        List<CvParam> paramList = mzidmod.getCvParam();
                        CvParam modParam = new CvParam();
                        double monoMass = intToModMap.get(fixedModification).getModMonoMass();
                        mzidmod.setMonoisotopicMassDelta(monoMass);
                        int mzidModLocation = index + 1;              //If second res is modified, index would return 1, but mzid position should be 2
                        mzidmod.setLocation(mzidModLocation);
                        boolean isMono = true;
                        ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, pepSeq.charAt(index));
                        if (unimod != null) {
                            modParam.setAccession("UNIMOD:" + unimod.getRecordId());
                            modParam.setCv(unimodCV);
                            modParam.setName(unimod.getTitle());
                        } else {
                            System.out.println("Error: modification with mass not recognized");
                            modParam.setName("unknown modification");
                            modParam.setCv(psiCV);
                            modParam.setAccession("MS:1001460");
                        }

                        paramList.add(modParam);
                        allMods.add(mzidmod);
                        index = pepSeq.indexOf(modifiedResidue, index + 1);
                    }
                }

                //Candidate N or C terminal mods
                if (modifiedResidues == null || modifiedResidues.isEmpty() || modifiedResidues.size() == 0) {

                    boolean isPepNTerminalMod = false;
                    boolean isPepCTerminalMod = false;
                    boolean isProtNTerminalMod = false;
                    boolean isProtCTerminalMod = false;

                    int modType = omod.getModType();
                    if (modType == 1 || modType == 2) {
                        isProtNTerminalMod = true;
                    } else if (modType == 5 || modType == 6) {
                        isPepNTerminalMod = true;
                    } else if (modType == 3 || modType == 4) {
                        isProtCTerminalMod = true;
                    } else if (modType == 7 || modType == 8) {
                        isPepCTerminalMod = true;
                    }
                    boolean isMono = true;
                    ModT unimod = null;

                    if (isPepNTerminalMod || isPepCTerminalMod) {

                        Modification mzidmod = new Modification();
                        List<CvParam> paramList = mzidmod.getCvParam();
                        CvParam modParam = new CvParam();
                        double monoMass = intToModMap.get(fixedModification).getModMonoMass();
                        mzidmod.setMonoisotopicMassDelta(monoMass);

                        if (isPepNTerminalMod) {
                            unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, '[');
                            mzidmod.setLocation(0);
                        } else if (isPepCTerminalMod) {
                            unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, ']');
                            mzidmod.setLocation(pepSeq.length() + 1);
                        }

                        if (unimod != null) {
                            modParam.setAccession("UNIMOD:" + unimod.getRecordId());
                            modParam.setCv(unimodCV);
                            modParam.setName(unimod.getTitle());
                        } else {
                            System.out.println("Error: modification with mass not recognized");
                            modParam.setName("unknown modification");
                            modParam.setCv(psiCV);
                            modParam.setAccession("MS:1001460");
                        }

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

            List<CvParam> paramList = mzidmod.getCvParam();
            CvParam modParam = new CvParam();
            OmssaModification omod = intToModMap.get(modNum);
            double monoMass = omod.getModMonoMass();

            mzidmod.setMonoisotopicMassDelta(monoMass);
            mzidmod.setLocation(modSite + 1);                  //+1 since Omssa counts from zero, mzid counts first position in peptide as 1

            boolean isMono = true;
            ModT unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, pepSeq.charAt(modSite));

            if (unimod != null) {
                mzidmod.getResidues().add("" + pepSeq.charAt(modSite));
            }

            if (unimod == null && modSite == 0) {
                //See if this is a possible N-terminal mod
                System.out.println("\tNot found, so look to see if it is N-terminal\n");
                unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, '[');
                mzidmod.setLocation(0);
            }
            if (unimod == null && modSite == pepSeq.length()) {
                //See if this is a possible C-terminal mod
                System.out.println("\tNot found, so look to see if it is C-terminal\n");
                unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, ']');
                mzidmod.setLocation(pepSeq.length() + 1);
            }

            if (unimod != null) {
                modParam.setAccession("UNIMOD:" + unimod.getRecordId());
                modParam.setCv(unimodCV);
                modParam.setName(unimod.getTitle());

            } else {
                System.out.println("Error: modification with mass not recognized");
                modParam.setName("unknown modification");
                modParam.setCv(psiCV);
                modParam.setAccession("MS:1001460");
            }

            //paramList.add(getModCV(mass));
            paramList.add(modParam);
            allMods.add(mzidmod);
        }
        return modString;
    }

    private void handleCVs() {

        //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
        //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
        //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>
        cvList = new CvList();
        List<Cv> localCvList = cvList.getCv();
        psiCV = new Cv();
        psiCV.setUri("https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        psiCV.setId(psiCvID);
//        psiCV.setVersion("2.25.0");
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

    public CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
    }

    public CvParam makeCvParam(String accession, String name, Cv cv, String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        cvParam.setValue(value);
        return cvParam;
    }

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

    public void handleAnalysisSoftware(String version) {
        analysisSoftwareList = new AnalysisSoftwareList();
        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("OMSSA");
        Param tempParam = new Param();
        //tempParam.setParamGroup(makeCvParam("MS:1001475","OMSSA",psiCV));
        tempParam.setParam(makeCvParam("MS:1001475", "OMSSA", psiCV));
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
        role.setCvParam(makeCvParam("MS:1001271", "researcher", psiCV));
        contactRole.setRole(role);

        provider.setContactRole(contactRole);

    }

    public void handleAuditCollection(String firstName, String lastName, String email, String address, String affiliationName) {
        auditCollection = new AuditCollection();
        //List<Contact> contactList = auditCollection.getContactGroup();
        List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();
        docOwner = new Person();
        docOwner.setId("PERSON_DOC_OWNER");
        docOwner.setFirstName(firstName);
        docOwner.setLastName(lastName);

        docOwner.getCvParam().add(makeCvParam("MS:1000587", "contact address", psiCV, address));

        //docOwner.setEmail(email);
        Organization org = new Organization();
        org.setId("ORG_DOC_OWNER");
        org.setName(affiliationName);
        org.getCvParam().add(makeCvParam("MS:1000586", "contact name", psiCV, address));

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

    public void handleAnalysisCollection(String activityDate) {
        analysisCollection = new AnalysisCollection();
        List<SpectrumIdentification> specIdentList = analysisCollection.getSpectrumIdentification();
        SpectrumIdentification specIdent = new SpectrumIdentification();
        specIdent.setId(specIdentID);
        specIdent.setSpectrumIdentificationList(siList);
        specIdent.setSpectrumIdentificationProtocol(siProtocol);
        List<SearchDatabaseRef> searchDBRefList = specIdent.getSearchDatabaseRef();
        SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
        searchDBRef.setSearchDatabase(searchDB);
        searchDBRefList.add(searchDBRef);

        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
        InputSpectra inputSpec = new InputSpectra();
        inputSpec.setSpectraData(spectraData);
        inputSpecList.add(inputSpec);

        specIdentList.add(specIdent);

    }

    public void handleAnalysisProtocolCollection() {

        //Boolean (parentIsMono, Boolean fragmentIsMono, SearchModification[] searchMods, String enzymeName, Double parTolPlus, Double parTolMinus, Double fragTolPlus, Double fragTolMinus);
       /*
         * <MSSearchSettings> <MSSearchSettings_precursorsearchtype>
         * <MSSearchType>0</MSSearchType>
         * </MSSearchSettings_precursorsearchtype>
         * <MSSearchSettings_productsearchtype> <MSSearchType>0</MSSearchType>
         * </MSSearchSettings_productsearchtype> <MSSearchSettings_ionstosearch>
         * <MSIonType>1</MSIonType> <MSIonType>4</MSIonType>
         * </MSSearchSettings_ionstosearch>
         * <MSSearchSettings_peptol>1.5</MSSearchSettings_peptol>
         * <MSSearchSettings_msmstol>0.8</MSSearchSettings_msmstol>
         * <MSSearchSettings_zdep> <MSZdependence>0</MSZdependence>
         * </MSSearchSettings_zdep>
         * <MSSearchSettings_cutoff>10</MSSearchSettings_cutoff>
         * <MSSearchSettings_cutlo>0</MSSearchSettings_cutlo>
         * <MSSearchSettings_cuthi>0.2</MSSearchSettings_cuthi>
         * <MSSearchSettings_cutinc>0.0005</MSSearchSettings_cutinc>
         * <MSSearchSettings_singlewin>20</MSSearchSettings_singlewin>
         * <MSSearchSettings_doublewin>14</MSSearchSettings_doublewin>
         * <MSSearchSettings_singlenum>2</MSSearchSettings_singlenum>
         * <MSSearchSettings_doublenum>2</MSSearchSettings_doublenum>
         * <MSSearchSettings_fixed> <MSMod>3</MSMod> </MSSearchSettings_fixed>
         * <MSSearchSettings_variable> <MSMod>1</MSMod>
         * </MSSearchSettings_variable> <MSSearchSettings_enzyme>
         * <MSEnzymes>0</MSEnzymes> </MSSearchSettings_enzyme>
         * <MSSearchSettings_missedcleave>1</MSSearchSettings_missedcleave>
         * <MSSearchSettings_hitlistlen>30</MSSearchSettings_hitlistlen>
         * <MSSearchSettings_db>D:/Software/Databases/Neospora_3rndTryp/Neo_rndTryp_3times.fasta</MSSearchSettings_db>
         * <MSSearchSettings_tophitnum>6</MSSearchSettings_tophitnum>
         * <MSSearchSettings_minhit>2</MSSearchSettings_minhit>
         * <MSSearchSettings_minspectra>4</MSSearchSettings_minspectra>
         * <MSSearchSettings_scale>1000</MSSearchSettings_scale>
         * <MSSearchSettings_maxmods>128</MSSearchSettings_maxmods>
         * <MSSearchSettings_chargehandling> <MSChargeHandle>
         * <MSChargeHandle_calcplusone> <MSCalcPlusOne>1</MSCalcPlusOne>
         * </MSChargeHandle_calcplusone> <MSChargeHandle_calccharge>
         * <MSCalcCharge>2</MSCalcCharge> </MSChargeHandle_calccharge>
         * <MSChargeHandle_mincharge>1</MSChargeHandle_mincharge>
         * <MSChargeHandle_maxcharge>3</MSChargeHandle_maxcharge>
         * <MSChargeHandle_considermult>3</MSChargeHandle_considermult>
         * <MSChargeHandle_plusone>0.95</MSChargeHandle_plusone>
         * <MSChargeHandle_maxproductcharge>2</MSChargeHandle_maxproductcharge>
         * </MSChargeHandle> </MSSearchSettings_chargehandling>
         * <MSSearchSettings_pseudocount>1</MSSearchSettings_pseudocount>
         * <MSSearchSettings_searchb1>1</MSSearchSettings_searchb1>
         * <MSSearchSettings_searchctermproduct>0</MSSearchSettings_searchctermproduct>
         * <MSSearchSettings_maxproductions>100</MSSearchSettings_maxproductions>
         * <MSSearchSettings_minnoenzyme>4</MSSearchSettings_minnoenzyme>
         * <MSSearchSettings_maxnoenzyme>40</MSSearchSettings_maxnoenzyme>
         * <MSSearchSettings_exactmass>1446.94</MSSearchSettings_exactmass>
         * <MSSearchSettings_settingid>0</MSSearchSettings_settingid>
         * <MSSearchSettings_iterativesettings> <MSIterativeSettings>
         * <MSIterativeSettings_researchthresh>0.01</MSIterativeSettings_researchthresh>
         * <MSIterativeSettings_subsetthresh>0</MSIterativeSettings_subsetthresh>
         * <MSIterativeSettings_replacethresh>0</MSIterativeSettings_replacethresh>
         * </MSIterativeSettings> </MSSearchSettings_iterativesettings>
         * <MSSearchSettings_precursorcull>0</MSSearchSettings_precursorcull>
         * <MSSearchSettings_infiles> <MSInFile>
         * <MSInFile_infile>D:/TestSpace/NeoTestMarch2011/55merge.mgf</MSInFile_infile>
         * <MSInFile_infiletype> <MSSpectrumFileType>7</MSSpectrumFileType>
         * </MSInFile_infiletype> </MSInFile> </MSSearchSettings_infiles>
         * <MSSearchSettings_outfiles> <MSOutFile>
         * <MSOutFile_outfile>D:/TestSpace/NeoTestMarch2011/55merge_omssa.omx</MSOutFile_outfile>
         * <MSOutFile_outfiletype> <MSSerialDataFormat>3</MSSerialDataFormat>
         * </MSOutFile_outfiletype> <MSOutFile_includerequest value="true"/>
         * </MSOutFile> </MSSearchSettings_outfiles>
         * <MSSearchSettings_nocorrelationscore>0</MSSearchSettings_nocorrelationscore>
         * <MSSearchSettings_probfollowingion>0.5</MSSearchSettings_probfollowingion>
         * <MSSearchSettings_nmethionine value="true"/>
         * <MSSearchSettings_automassadjust>1</MSSearchSettings_automassadjust>
         * <MSSearchSettings_noprolineions></MSSearchSettings_noprolineions>
         * </MSSearchSettings>
         */
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

        int msSearchType = settings.MSSearchSettings_precursorsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact

        if (msSearchType == 0) {
            cvParamList.add(makeCvParam("MS:1001211", "parent mass type mono", psiCV));
        } else if (msSearchType == 1) {
            cvParamList.add(makeCvParam("MS:1001212", "parent mass type average", psiCV));
        } else if (msSearchType == 2) {
            cvParamList.add(makeCvParam("MS:No_acc", "monoisotopic N15", psiCV));
            System.out.println("Warning: No CV term for monoisotopic N15 search type");
        } else if (msSearchType == 3) {
            cvParamList.add(makeCvParam("MS:No_acc", "exact", psiCV));
            System.out.println("Warning: No CV term for exact mass search type");
        } else {
            System.out.println("Error search type not recognised");

        }

        int prodSearchType = settings.MSSearchSettings_productsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact

        if (prodSearchType == 0) {
            cvParamList.add(makeCvParam("MS:1001256", "fragment mass type mono", psiCV));
        } else if (prodSearchType == 1) {
            cvParamList.add(makeCvParam("MS:1001255", "fragment mass type average", psiCV));
        } else if (prodSearchType == 2) {
            cvParamList.add(makeCvParam("MS:No_acc", "monoisotopic N15", psiCV));
            System.out.println("Warning: No CV term for monoisotopic N15 search type");
        } else if (prodSearchType == 3) {
            cvParamList.add(makeCvParam("MS:No_acc", "exact", psiCV));
            System.out.println("Warning: No CV term for exact mass search type");
        } else {
            System.out.println("Error search type not recognised");
        }

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams.getSearchModification();

        for (int fixedMod : settings.MSSearchSettings_fixed.MSMod) {
            OmssaModification omod = intToModMap.get(fixedMod);
            SearchModification searchMod = handleOmssaMod(omod, true);
            searchModList.add(searchMod);
        }

        for (int varMod : settings.MSSearchSettings_variable.MSMod) {
            OmssaModification omod = intToModMap.get(varMod);
            SearchModification searchMod = handleOmssaMod(omod, false);
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
        List<Integer> msEnzymeList = settings.MSSearchSettings_enzyme.MSEnzymes;
        for (Integer msEnzyme : msEnzymeList) {
            //OmssaEnumerators omssaEnums = new OmssaEnumerators();
            Enzyme enzyme = getEnzyme(OmssaEnumerators.getEnzymeAsText(msEnzyme), settings.MSSearchSettings_missedcleave);
            enzymeList.add(enzyme);
        }

        Tolerance fragTol = new Tolerance();
        Tolerance parTol = new Tolerance();

        List<CvParam> fragCvList = fragTol.getCvParam();
        CvParam fragCvPlus = getCvParamWithMassUnits(true);
        CvParam fragCvMinus = getCvParamWithMassUnits(true);


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
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        List<CvParam> parCvList = parTol.getCvParam();
        CvParam parCvPlus = getCvParamWithMassUnits(true);
        CvParam parCvMinus = getCvParamWithMassUnits(true);

        parCvPlus.setAccession("MS:1001412");
        parCvPlus.setName("search tolerance plus value");
        parCvMinus.setAccession("MS:1001413");
        parCvMinus.setName("search tolerance minus value");
        parCvPlus.setValue("" + settings.MSSearchSettings_peptol);
        parCvMinus.setValue("" + settings.MSSearchSettings_peptol);
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

        cvParamList.add(makeCvParam("MS:1001494", "no threshold", psiCV));
        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />

        sipList.add(siProtocol);

    }

    private SearchModification handleOmssaMod(OmssaModification omod, boolean isFixed) {
        SearchModification searchMod = new SearchModification();
        double monoMass = omod.getModMonoMass();
        List<String> residues = omod.getModResidues();

        boolean isPepNTerminalMod = false;
        boolean isPepCTerminalMod = false;
        boolean isProtNTerminalMod = false;
        boolean isProtCTerminalMod = false;

        int modType = omod.getModType();

        if (modType == 1 || modType == 2) {
            isProtNTerminalMod = true;
        } else if (modType == 5 || modType == 6) {
            isPepNTerminalMod = true;
        } else if (modType == 3 || modType == 4) {
            isProtCTerminalMod = true;
        } else if (modType == 7 || modType == 8) {
            isPepCTerminalMod = true;
        }
        boolean isMono = true;
        ModT unimod;
        if (isPepNTerminalMod || isProtNTerminalMod) {
            residues.add("[");
            unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, residues);
            residues.remove(residues.size() - 1); //Remove the ] or [ char
        } else if (isPepCTerminalMod || isProtCTerminalMod) {
            residues.add("]");
            unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, residues);
            residues.remove(residues.size() - 1); //Remove the ] or [ char
        } else {
            unimod = unimodDoc.getModByMass(monoMass, unimodMassError, isMono, residues);
        }

        searchMod.setFixedMod(isFixed);

        List<CvParam> modCvParamList = searchMod.getCvParam();
        CvParam modParam = new CvParam();

        if (unimod != null) {
            modParam = makeCvParam("UNIMOD:" + unimod.getRecordId(), unimod.getTitle(), unimodCV);

        } else {
            modParam.setName("unknown modification");
            modParam.setCv(psiCV);
            modParam.setAccession("MS:1001460");
        }

        modCvParamList.add(modParam);
        searchMod.setMassDelta(new Float(monoMass));
        List<String> residueList = searchMod.getResidues();

        if (isPepNTerminalMod) {
            modCvParamList.add(makeCvParam("MS:1001189", "modification specificity peptide N-term", psiCV));
        } else if (isPepCTerminalMod) {
            modCvParamList.add(makeCvParam("MS:1001190", "modification specificity peptide C-term", psiCV));
        } else if (isProtNTerminalMod) {
            modCvParamList.add(makeCvParam("MS:1002057", "modification specificity protein N-term", psiCV));
        } else if (isProtCTerminalMod) {
            modCvParamList.add(makeCvParam("MS:1002058", "modification specificity protein C-term", psiCV));
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

    public void handleInputs(MSResponse response, MSSearchSettings settings) {

        inputs = new Inputs();
        List<SearchDatabase> searchDBList = inputs.getSearchDatabase();

        searchDB = new SearchDatabase();
        searchDB.setId(searchDBID);
        searchDB.setNumDatabaseSequences((long) response.MSResponse_dbversion);
        //<cvParam accession="MS:1001401" name="xtandem xml file" cvRef="PSI-MS"/>

        UserParam param = new UserParam();
        param.setName(settings.MSSearchSettings_db);
        Param tempParam = new Param();
        tempParam.setParam(param);
        searchDB.setDatabaseName(tempParam);
        searchDB.setLocation(settings.MSSearchSettings_db);

        FileFormat ff = new FileFormat();
        ff.setCvParam(makeCvParam("MS:1001348", "FASTA format", psiCV));
        searchDB.setFileFormat(ff);   //TODO - this should not be hard coded <cvParam accession="MS:1001348" name="FASTA format" cvRef="PSI-MS"/>
        searchDBList.add(searchDB);

        List<SourceFile> sourceFileList = inputs.getSourceFile();
        SourceFile sourceFile = new SourceFile();
        sourceFile.setLocation(inputOmssaFile);
        sourceFile.setId(sourceFileID);
        ff = new FileFormat();
        ff.setCvParam(makeCvParam("MS:1001400", "OMSSA xml file", psiCV));

        sourceFile.setFileFormat(ff);
        sourceFileList.add(sourceFile);

        List<SpectraData> spectraDataList = inputs.getSpectraData();
        spectraData = new SpectraData();
        SpectrumIDFormat sif = new SpectrumIDFormat();
        sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", psiCV));
        spectraData.setSpectrumIDFormat(sif);

        int fileType = settings.MSSearchSettings_infiles.MSInFile.MSInFile_infiletype.MSSpectrumFileType;
        ff = new FileFormat();
        ff.setCvParam(getFileFormatCVParam(fileType));
        spectraData.setFileFormat(ff);

        spectraData.setId(spectraDataID);
        spectraData.setLocation(settings.MSSearchSettings_infiles.MSInFile.MSInFile_infile);
        spectraDataList.add(spectraData);

        /*
         * <fileFormat> <cvParam accession="MS:1001062" name="Mascot MGF file"
         * cvRef="PSI-MS" /> </fileFormat> <spectrumIDFormat> <cvParam
         * accession="MS:1000774" name="multiple peak list nativeID format"
         * cvRef="PSI-MS" /> </spectrumIDFormat>
         */
        //<note label="modelling, total proteins used">22348</note>
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

    public void writeMzidFile(String output) {
        try {
            FileWriter writer = null;
            if (output.equals("")) {
                writer = new FileWriter("55merge_omssa.mzid");
            } else {
                writer = new FileWriter(output);
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
            SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
            analysisSoftware.setName(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
            analysisSoftware.setId(this.getClass().getSimpleName() + "_" + dateFormat.format(date));

            Param param = new Param();
            param.setParam(makeCvParam("MS:1002237", "mzidLib", psiCV));
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
            writer.write(m.createProteinDetectionListStartTag("PDL_1", null) + "\n");

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
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

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

    public CvParam getFragmentCVParam(int iontype) {

        CvParam cvParam = new CvParam();
        cvParam.setCv(psiCV);

        /*
         * id: MS:1001108 name: param: a ion id: MS:1001118 name: param: b ion
         * id: MS:1001119 name: param: c ion id: MS:1001146 name: param: a
         * ion-NH3 id: MS:1001148 name: param: a ion-H2O id: MS:1001149 name:
         * param: b ion-NH3 id: MS:1001150 name: param: b ion-H2O id: MS:1001151
         * name: param: y ion-NH3 id: MS:1001152 name: param: y ion-H2O id:
         * MS:1001257 name: param: v ion id: MS:1001258 name: param: d ion id:
         * MS:1001259 name: param: immonium ion id: MS:1001260 name: param: w
         * ion id: MS:1001261 name: param: x ion id: MS:1001262 name: param: y
         * ion id: MS:1001263 name: param: z ion id: MS:1001406 name: param:
         * internal yb ion id: MS:1001407 name: param: internal ya ion id:
         * MS:1001408 name: param: z+1 ion id: MS:1001409 name: param: z+2 ion
         *
         * index = 0; MSIonType.put(new Integer(index++), "a");
         * MSIonType.put(new Integer(index++), "b"); MSIonType.put(new
         * Integer(index++), "c"); MSIonType.put(new Integer(index++), "x");
         * MSIonType.put(new Integer(index++), "y"); MSIonType.put(new
         * Integer(index++), "z"); MSIonType.put(new Integer(index++),
         * "parent"); MSIonType.put(new Integer(index++), "internal");
         * MSIonType.put(new Integer(index++), "immonium"); MSIonType.put(new
         * Integer(index++), "unknown"); MSIonType.put(new Integer(index++),
         * "adot"); MSIonType.put(new Integer(index++), "x-CO2");
         * MSIonType.put(new Integer(index++), "adot-CO2"); MSIonType.put(new
         * Integer(index++), "max");
         *
         */
        switch (iontype) {
            case 0:
                cvParam.setAccession("MS:1001108");
                cvParam.setName("param: a ion");
                break;
            case 1:
                cvParam.setAccession("MS:1001118");
                cvParam.setName("param: b ion");
                break;
            case 2:
                cvParam.setAccession("MS:1001108");
                cvParam.setName("param: a ion");
                break;
            case 3:
                cvParam.setAccession("MS:1001261");
                cvParam.setName("param: x ion");
                break;
            case 4:
                cvParam.setAccession("MS:1001262");
                cvParam.setName("param: y ion");
                break;
            case 5:
                cvParam.setAccession("MS:1001263");
                cvParam.setName("param: z ion");
                break;
            case 6:
                cvParam.setAccession("MS:No CV term");
                cvParam.setName("param: parent");
                break;
            case 7:
                cvParam.setAccession("MS:No CV term");
                cvParam.setName("param: internal");
                break;
            case 8:
                cvParam.setAccession("MS:MS:1001259");
                cvParam.setName("param: immonium");
                break;
            default:
                cvParam.setAccession("MS:Unknown");
                cvParam.setName("param: unknown ion type");

        }

        return cvParam;

    }

    public CvParam getFileFormatCVParam(int omssaFileType) {

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
        cvParam.setCv(psiCV);

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

    public double[] getMatchedIon(MSSpectrum spectrum, double ion_mz) {

        double[] matchedPeak = new double[2];
        MSSpectrum_mz spec_mz = spectrum.MSSpectrum_mz;
        List<Integer> mzVals = spec_mz.MSSpectrum_mz_E;
        MSSpectrum_abundance spec_int = spectrum.MSSpectrum_abundance;
        List<Integer> abundValues = spec_int.MSSpectrum_abundance_E;

        double error = 100;

        int i = 0;
        int foundPos = -1;

        double foundMZ = -1;
        for (int mzVal : mzVals) {
            double mz = (double) mzVal / response_scale;

            double tempError = mz - ion_mz;

            if (Math.abs(tempError) < Math.abs(error)) {
                error = tempError;
                foundPos = i;
                foundMZ = mz;
            }

            i++;
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
            double intensity_scale = spectrum.MSSpectrum_iscale;
            foundInt = (double) (abundValues.get(foundPos));
            foundInt = foundInt / intensity_scale;

        }

        matchedPeak[0] = foundMZ;
        matchedPeak[1] = foundInt;
        return matchedPeak;
    }

    public Enzyme getEnzyme(String omssaEnzyme, int missedCleavage) {

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
            cvParamList.add(makeCvParam("MS:1001251", "Trypsin", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Arg-C")) {
            cvParamList.add(makeCvParam("MS:1001303", "Arg-C", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("CNBr")) {
            cvParamList.add(makeCvParam("MS:1001307", "CNBr", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Chymotrypsin")) {
            cvParamList.add(makeCvParam("MS:1001306", "Chymotrypsin", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Formic Acid")) {
            cvParamList.add(makeCvParam("TODO", "TODO", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Lys-C")) {
            cvParamList.add(makeCvParam("MS:1001309", "Lys-C", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Lys-C, no P rule")) {
            cvParamList.add(makeCvParam("MS:1001310", "Lys-C/P", psiCV));
        } else if (omssaEnzyme.equalsIgnoreCase("Pepsin A")) {
            cvParamList.add(makeCvParam("TODO", "TODO", psiCV));
        } else {
            cvParamList.add(makeCvParam("TODO", "TODO", psiCV));
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
}
