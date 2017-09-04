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
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.proteinms.omxparser.OmssaOmxFile;
import de.proteinms.omxparser.util.MSHitSet;
import de.proteinms.omxparser.util.MSHitSet_error;
import de.proteinms.omxparser.util.MSHitSet_hits;
import de.proteinms.omxparser.util.MSHits;
import de.proteinms.omxparser.util.MSHits_mods;
import de.proteinms.omxparser.util.MSHits_mzhits;
import de.proteinms.omxparser.util.MSHits_pephits;
import de.proteinms.omxparser.util.MSIonAnnot;
import de.proteinms.omxparser.util.MSMZHit;
import de.proteinms.omxparser.util.MSMZHit_annotation;
import de.proteinms.omxparser.util.MSMZHit_ion;
import de.proteinms.omxparser.util.MSPepHit;
import de.proteinms.omxparser.util.MSRequest;
import de.proteinms.omxparser.util.MSRequest_settings;
import de.proteinms.omxparser.util.MSResponse;
import de.proteinms.omxparser.util.MSSearch;
import de.proteinms.omxparser.util.MSSearchSettings;
import de.proteinms.omxparser.util.MSSearch_request;
import de.proteinms.omxparser.util.MSSearch_response;
import de.proteinms.omxparser.util.MSSpectrum;
import de.proteinms.omxparser.util.MSSpectrum_ids;
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
import uk.ac.ebi.jmzidml.model.mzidml.FragmentArray;
import uk.ac.ebi.jmzidml.model.mzidml.Fragmentation;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.Measure;
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
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.Tolerance;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.ReadUnimod;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 02-Sep-2017 00:12:41
 */
public class Omssa2mzidMzidContainer implements MzidContainer {

    private String modsFile;
    private String userModsFile;
    private MzIdentMLVersion version;
    private OmssaOmxFile omxFile;
    private Boolean outputFragmentation;
    private String decoyRegularExpression;
    private AnalysisSoftware analysisSoftwareOmssa;
    private SequenceCollection sequenceCollection;
    private SpectrumIdentificationList siList;

    private Person docOwner;

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

    private ReadUnimod unimodDoc;
    private MSSearchSettings settings;

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

    public Omssa2mzidMzidContainer(String input,
                                   Boolean outputFrags,
                                   String decoyRegex,
                                   String omssaModsFile,
                                   String omssaUserModsFile,
                                   MzIdentMLVersion ver) {
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

        init();

        parseFile(this.omxFile);

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
            FileWriter writer = new FileWriter(filename);
            writer.write(text);
            writer.close();
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

        Map<MSSpectrum, HashSet<String>> specToPepMap = omxFile
                .getSpectrumToPeptideMap();
        MSSearch mssearch = omxFile.getParserResult();
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

        MSResponse response = null;
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

        Map<MSSpectrum, MSHitSet> specToHitMap = omxFile
                .getSpectrumToHitSetMap();
        intToModMap = omxFile.getModifications();

        foundProts = new HashMap<>();
        pepProtMap = new HashMap<>();
        //lookup to get a peptide by peptideseq_varmods_fixedmods
        peptideLookup = new HashMap<>();
        pepEvidLookup = new HashMap<>();
        uniquePeps = new HashMap<>();
        sequenceCollection = new SequenceCollection();

        List<PeptideEvidence> peptideEvidenceList = sequenceCollection
                .getPeptideEvidence();

        siList = new SpectrumIdentificationList();
        siList.setId(SI_LIST_ID);

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
        mzMeasure.setId(MEASURE_MZ_ID);
        List<CvParam> cvParamList = mzMeasure.getCvParam();
        cvParamList.add(MzidLibUtils
                .makeCvParam("MS:1001225", "product ion m/z", CvConstants.PSI_CV,
                             "MS:1000040",
                             "m/z", CvConstants.UNIT_CV));
        Measure intMeasure = new Measure();
        intMeasure.setId(MEASURE_INT_ID);
        cvParamList = intMeasure.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001226",
                                                 "product ion intensity",
                                                 CvConstants.PSI_CV,
                                                 "MS:1000131",
                                                 "number of counts",
                                                 CvConstants.UNIT_CV));
        Measure errorMeasure = new Measure();
        errorMeasure.setId(MEASURE_ERROR_ID);
        cvParamList = errorMeasure.getCvParam();
        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001227",
                                                 "product ion m/z error",
                                                 CvConstants.PSI_CV,
                                                 "MS:1000040", "m/z",
                                                 CvConstants.UNIT_CV));

        measureList.add(mzMeasure);
        measureList.add(intMeasure);
        measureList.add(errorMeasure);

        siList.setFragmentationTable(fragTable);

        List<SpectrumIdentificationResult> specIdentResults = siList
                .getSpectrumIdentificationResult();
        List<DBSequence> dbSeqList = sequenceCollection.getDBSequence();
        List<Peptide> peptideList = sequenceCollection.getPeptide();
        SpectraData spectraData = new SpectraData();
        spectraData.setId(SPECTRA_DATA_ID);

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
                MSHitSet_error msError = msHitSet.MSHitSet_error;

                specIdentRes.setSpectraData(spectraData);
                specIdentRes
                        .setSpectrumID("index=" + spectrum.MSSpectrum_number);

                MSSpectrum_ids specIDs = spectrum.MSSpectrum_ids;

                if (!specIDs.MSSpectrum_ids_E.isEmpty()) {
                    List<CvParam> sir_cvParamList = specIdentRes.getCvParam();
                    for (String id : specIDs.MSSpectrum_ids_E) {
                        CvParam cvp = MzidLibUtils.makeCvParam("MS:1000796",
                                                               "spectrum title",
                                                               CvConstants.PSI_CV,
                                                               id);
                        sir_cvParamList.add(cvp);
                    }
                }

                int rank = 1;

                //Counter used to create unique ID for SpectrumIdentificationItem
                int siiCounter = 1;
                for (MSHits hits : msHits.MSHits) {

                    SpectrumIdentificationItem sii
                            = new SpectrumIdentificationItem();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = sii
                            .getPeptideEvidenceRef();
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

                    long expMZ = java.lang.Math.round(hits.MSHits_mass
                            / hits.MSHits_charge);
                    long theoMass = java.lang.Math.round(hits.MSHits_theomass
                            / hits.MSHits_charge);

                    if (evalue > maxEvalue) {
                        maxEvalue = evalue;
                    }

                    sii.setExperimentalMassToCharge((double) expMZ / 1000);
                    sii.setCalculatedMassToCharge((double) theoMass / 1000);

                    cvParamList = sii.getCvParam();
                    cvParamList.add(CvConstants.OMSSA_EVALUE);
                    cvParamList.add(CvConstants.OMSSA_PVALUE);

                    MSHits_mods mods = hits.MSHits_mods;

                    Peptide mzidPep = new Peptide();
                    mzidPep.setPeptideSequence(pepSeq);

                    String modString = convertVarMods(mzidPep, mods);
                    convertFixedMods(mzidPep);                       //Uses the fixed mods in the search settings to guess the fixed modifications

                    List<IonType> ionTypeList = null;
                    if (outputFragmentation) {
                        MSHits_mzhits mzhits = hits.MSHits_mzhits;
                        Map<String, List<MSMZHit>> mapIonsToHits
                                = new HashMap<>();

                        //ionTypeList = sii.getFragmentation();
                        Fragmentation frag = new Fragmentation();
                        sii.setFragmentation(frag);
                        ionTypeList = frag.getIonType();

                        for (MSMZHit mzhit : mzhits.MSMZHit) {

                            int index = mzhit.MSMZHit_index;
                            int ionCharge = mzhit.MSMZHit_charge;
                            double ionMz = (double) mzhit.MSMZHit_mz;
                            int number = mzhit.MSMZHit_number;

                            MSMZHit_annotation ionAnnot
                                    = mzhit.MSMZHit_annotation;

                            MSIonAnnot ionA = ionAnnot.MSIonAnnot;
                            double massError = ionA.MSIonAnnot_massdiff;

                            MSMZHit_ion ion = mzhit.MSMZHit_ion;
                            int ionType = ion.MSIonType;
                            String ionString = OmssaEnumerators.
                                    getMSIonTypeAsText(ionType);
                            //System.out.println("\t" + index + "\t" + number + "\t" + ionCharge + "\t" + ionString + "\t" + ionMz + "\t" + massError);

                            //Logic is slightly different from Tandem parser, since all ion types appear to be mixed up in one structure
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
                            List<FragmentArray> fragmentList = mzidIon
                                    .getFragmentArray();

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
                                String ionString = OmssaEnumerators
                                        .getMSIonTypeAsText(ionType);

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

                                double theoMZ = (double) (mzhit.MSMZHit_mz)
                                        / responseScale;
                                double[] matchedPeak = getMatchedIon(spectrum,
                                                                     theoMZ);
                                mzValues.add((float) matchedPeak[0]);
                                intValues.add((float) matchedPeak[1]);
                                errorValues.add(
                                        (float) (matchedPeak[0] - theoMZ));
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
                            protAcc = defline.substring(0, defline.indexOf(
                                                        defline_regex));
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

                            db_cvParamList.add(MzidLibUtils.makeCvParam(
                                    "MS:1001088", "protein description", CvConstants.PSI_CV,
                                    desc));
                            dbSeqList.add(dbSeq);
                        } else {
                            dbSeq = foundProts.get(protAcc);
                        }

                        String testPepMods = pepSeq + "_" + modString + "_"
                                + start + "_" + end;
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
                                if (protAcc.contains(this.decoyRegularExpression)) {
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

    @Override
    public AnalysisCollection getAnalysisCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        AnalysisProtocolCollection analysisProtocolCollection
                = new AnalysisProtocolCollection();
        List<SpectrumIdentificationProtocol> sipList
                = analysisProtocolCollection.getSpectrumIdentificationProtocol();

        SpectrumIdentificationProtocol siProtocol
                = new SpectrumIdentificationProtocol();
        siProtocol.setId(SI_PROTOCOL_ID);
        siProtocol.setAnalysisSoftware(analysisSoftwareOmssa);

        //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
        Param tempParam = new Param();
        tempParam.setParam(CvConstants.MS_MS_SEARCH);
        siProtocol.setSearchType(tempParam);

        //List<CvParam> cvParamList = siProtocol.getAdditionalSearchCvParams();
        ParamList paramList = siProtocol.getAdditionalSearchParams();
        if (paramList == null) {
            paramList = new ParamList();
            siProtocol.setAdditionalSearchParams(paramList);
        }
        List<CvParam> cvParamList = paramList.getCvParam();

        int msSearchType
                = settings.MSSearchSettings_precursorsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact

        if (msSearchType == 0) {
            cvParamList.add(CvConstants.PARENT_MASS_TYPE_MONO);
        } else if (msSearchType == 1) {
            cvParamList.add(CvConstants.PARENT_MASS_TYPE_AVERAGE);
        } else if (msSearchType == 2) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:No_acc",
                                                     "monoisotopic N15",
                                                     CvConstants.PSI_CV));
            System.out.println(
                    "Warning: No CV term for monoisotopic N15 search type");
        } else if (msSearchType == 3) {
            cvParamList.add(MzidLibUtils
                    .makeCvParam("MS:No_acc", "exact", CvConstants.PSI_CV));
            System.out.println("Warning: No CV term for exact mass search type");
        } else {
            System.out.println("Error search type not recognised");

        }

        int prodSearchType
                = settings.MSSearchSettings_productsearchtype.MSSearchType;    //with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact

        if (prodSearchType == 0) {
            cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_MONO);
        } else if (prodSearchType == 1) {
            cvParamList.add(CvConstants.FRAGMENT_MASS_TYPE_AVERAGE);
        } else if (prodSearchType == 2) {
            cvParamList.add(MzidLibUtils.makeCvParam("MS:No_acc",
                                                     "monoisotopic N15",
                                                     CvConstants.PSI_CV));
            System.out.println(
                    "Warning: No CV term for monoisotopic N15 search type");
        } else if (prodSearchType == 3) {
            cvParamList.add(MzidLibUtils
                    .makeCvParam("MS:No_acc", "exact", CvConstants.PSI_CV));
            System.out.println("Warning: No CV term for exact mass search type");
        } else {
            System.out.println("Error search type not recognised");
        }

        // Add "no special processing" cv term if this is mzid 1.2 version
        if (version.equals(MzIdentMLVersion.Version_1_2)) {
            cvParamList.add(CvConstants.NO_SPECIAL_PROCESSING);
        }

        ModificationParams modParams = new ModificationParams();
        List<SearchModification> searchModList = modParams
                .getSearchModification();

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
            Enzyme enzyme
                    = getEnzyme(OmssaEnumerators.getEnzymeAsText(msEnzyme),
                                settings.MSSearchSettings_missedcleave);
            enzymeList.add(enzyme);
        }

        Tolerance fragTol = new Tolerance();
        Tolerance parTol = new Tolerance();

        List<CvParam> fragCvList = fragTol.getCvParam();
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
        fragCvList.add(fragCvPlus);
        fragCvList.add(fragCvMinus);

        List<CvParam> parCvList = parTol.getCvParam();
        CvParam parCvPlus = MzidLibUtils.getCvParamWithMassUnits(true);
        CvParam parCvMinus = MzidLibUtils.getCvParamWithMassUnits(true);

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

        cvParamList.add(MzidLibUtils.makeCvParam("MS:1001494", "no threshold",
                                                 CvConstants.PSI_CV));
        //<cvParam accession="MS:1001494" name="no threshold" cvRef="PSI-MS" />

        sipList.add(siProtocol);
        
        return analysisProtocolCollection;
    }

    @Override
    public AnalysisSampleCollection getAnalysisSampleCollection() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
    public Inputs getInputs() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public ProteinDetectionList getProteinDetectionList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
        SequenceCollection sequenceCollection = new SequenceCollection();

        List<PeptideEvidence> peptideEvidenceList = sequenceCollection
                .getPeptideEvidence();
        List<DBSequence> dbSeqList = sequenceCollection.getDBSequence();
        List<Peptide> peptideList = sequenceCollection.getPeptide();

        return sequenceCollection;
    }

    @Override
    public SpectrumIdentificationList getSpectrumIdentificationList() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private void init() {

        //analysisSoftwareXtandem
        this.analysisSoftwareOmssa = MzidLibUtils.createAnalysisSoftware(
                "xtandem",
                ANALYSIS_SOFT_ID,
                CvConstants.OMSSA,
                "");

        //docOwner
        docOwner = createDocOwner("firstName", "secondName", "address");

        unimodDoc = new ReadUnimod();
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

}
