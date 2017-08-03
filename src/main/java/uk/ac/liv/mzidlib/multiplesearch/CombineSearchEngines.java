
package uk.ac.liv.mzidlib.multiplesearch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Vector;

import uk.ac.ebi.jmzidml.MzIdentMLElement;

import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.CombinedSearchEngines2Mzid;
import uk.ac.liv.mzidlib.fdr.FalseDiscoveryRate;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

public class CombineSearchEngines {

    private static boolean verbose = true;

    private static MzIdentMLVersion mzidVer;
    int noOfSearchEngines = 0;

    Map<String, List<List<Object>>> combinedResultContainer;
    FdrAndMzIdentInformationContainer[] singleFDRInformation;
    String[] searchEngineIdentifiers; // The ordered entries correspond to the array entries in singleFDRInformation[] 
    String[] shortRepresentationOfSearchEngineIdentiifers; // stores 'fdr' for 'mascot', 'o' for 'omssa', 't' for tandem
    int lineNumber = 1;
    // FG
    private CombinedSearchEngines2Mzid combinedSearchEngines2Mzid
            = new CombinedSearchEngines2Mzid();
    private String decoyRegex = null;
    //metadata
    private FalseDiscoveryRate fdr;
    private AnalysisSoftwareList analysisSoftwareList;
    private AuditCollection auditCollection;
    private Provider provider;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private CvList cvList;
    // private SequenceCollection sequenceCollection;
    private AnalysisCollection analysisCollection;
    private Inputs inputs;
    private int sirCounter = 0;
    static String csvFileName;

    private Map<String, SpectraData> spectraDataHashMap = new HashMap<>();
    private Map<String, String> spectraIDLocation = new HashMap<>();

    private AnalysisProtocolCollection analysisProtocolCollectionTandem;

    public CombineSearchEngines() {
    }

    public CombineSearchEngines(String[] searchEngineNames) {

        int totalSearchEngines = searchEngineNames.length;

        noOfSearchEngines = totalSearchEngines;
        searchEngineIdentifiers = searchEngineNames.clone();
        shortRepresentationOfSearchEngineIdentiifers
                = new String[totalSearchEngines];
        singleFDRInformation
                = new FdrAndMzIdentInformationContainer[totalSearchEngines];
        for (int i = 0; i < totalSearchEngines; i++) {
            singleFDRInformation[i] = new FdrAndMzIdentInformationContainer();
        }

        int totalContainersNeeded = 0;

        switch (totalSearchEngines) {
            case 2:
                totalContainersNeeded = 3;
                break;
            case 3:
                totalContainersNeeded = 7;
                break;
            default:
                System.out.println(
                        "Min 2 and maximum 3 search engines allowed in the current version");
        }

        String[] namesOfContainers = new String[totalContainersNeeded];
        //	combinedResultContainer = new HashMap<String,ArrayList<ArrayList<String>>>(totalContainersNeeded);
        combinedResultContainer = new HashMap<>(totalContainersNeeded);

        switch (totalSearchEngines) {
            case 2:
                if ((searchEngineNames[0].equalsIgnoreCase("s1")
                        && searchEngineNames[1].equalsIgnoreCase("s2"))
                        || (searchEngineNames[1].equalsIgnoreCase("s1")
                        && searchEngineNames[0].equalsIgnoreCase("s2"))) {
                    namesOfContainers[0] = "1";
                    namesOfContainers[1] = "2";
                    namesOfContainers[2] = "12";
                }
                if ((searchEngineNames[0].equalsIgnoreCase("s1")
                        && searchEngineNames[1].equalsIgnoreCase("s3"))
                        || (searchEngineNames[1].equalsIgnoreCase("s1")
                        && searchEngineNames[0].equalsIgnoreCase("s3"))) {
                    namesOfContainers[0] = "1";
                    namesOfContainers[1] = "3";
                    namesOfContainers[2] = "13";
                }
                if ((searchEngineNames[0].equalsIgnoreCase("s2")
                        && searchEngineNames[1].equalsIgnoreCase("s3"))
                        || (searchEngineNames[1].equalsIgnoreCase("s2")
                        && searchEngineNames[0].equalsIgnoreCase("s3"))) {
                    namesOfContainers[0] = "2";
                    namesOfContainers[1] = "3";
                    namesOfContainers[2] = "23";
                }
                break;
            case 3:
                namesOfContainers[0] = "1";
                namesOfContainers[1] = "2";
                namesOfContainers[2] = "3";
                namesOfContainers[3] = "12";
                namesOfContainers[4] = "13";
                namesOfContainers[5] = "23";
                namesOfContainers[6] = "123";
            default:
                System.out.println(
                        "Min 2 and maximum 3 search engines allowed in the current version");
        }

        for (int i = 0; i < totalContainersNeeded; i++) {
            combinedResultContainer.put(namesOfContainers[i],
                                        new ArrayList<List<Object>>());
        }

        for (int i = 0; i < searchEngineIdentifiers.length; i++) {
            if (searchEngineIdentifiers[i].equalsIgnoreCase("s1")) {
                shortRepresentationOfSearchEngineIdentiifers[i] = "1";
            }
            if (searchEngineIdentifiers[i].equalsIgnoreCase("s2")) {
                shortRepresentationOfSearchEngineIdentiifers[i] = "2";
            }
            if (searchEngineIdentifiers[i].equalsIgnoreCase("s3")) {
                shortRepresentationOfSearchEngineIdentiifers[i] = "3";
            }
        }

    }

    public void computeFDRForSingleSearchEngine(int i, String xmlToRead,
                                                String searchEngine,
                                                FdrAndMzIdentInformationContainer fdrObj,
                                                int decoyRatio,
                                                String decoyRegex, String cvterm,
                                                String betterScore)
            throws Exception {

        fdr
                = new FalseDiscoveryRate(new File(xmlToRead), decoyRatio,
                                         decoyRegex, cvterm, Boolean.valueOf(
                                                 betterScore));

        //FG
        this.decoyRegex = decoyRegex;
        if (searchEngine.equals("s1")) {
            combinedSearchEngines2Mzid.setMascotDbSequenceHashMap(fdr.
                    getdBSequenceHashMap());
            combinedSearchEngines2Mzid.setMascotPeptideEvidenceHashMap(fdr.
                    getPeptideEvidenceMap());
            combinedSearchEngines2Mzid.setMascotPeptideHashMap(fdr.
                    getPeptideHashMap());
            combinedSearchEngines2Mzid.
                    setMascotSpectrumIdentificationItemHashMap(fdr.
                            getSpectrumIdentificationItemHashMap());
//            combinedSearchEngines2Mzid.setMascotSpectrumIdentificationResultHashMap(fdr.getSpectrumIdentificationResultHashMap());
            combinedSearchEngines2Mzid.getUnimodHashmap().putAll(fdr.
                    getUnimodHashmap());
        } else if (searchEngine.equals("s2")) {
            combinedSearchEngines2Mzid.setOmssaDbSequenceHashMap(fdr.
                    getdBSequenceHashMap());
            combinedSearchEngines2Mzid.setOmssaPeptideEvidenceHashMap(fdr.
                    getPeptideEvidenceMap());
            combinedSearchEngines2Mzid.setOmssaPeptideHashMap(fdr.
                    getPeptideHashMap());
            combinedSearchEngines2Mzid.
                    setOmssaSpectrumIdentificationItemHashMap(fdr.
                            getSpectrumIdentificationItemHashMap());
//            combinedSearchEngines2Mzid.setOmssaSpectrumIdentificationResultHashMap(fdr.getSpectrumIdentificationResultHashMap());
            combinedSearchEngines2Mzid.getUnimodHashmap().putAll(fdr.
                    getUnimodHashmap());

        } else if (searchEngine.equals("s3")) {
            combinedSearchEngines2Mzid.setTandemDbSequenceHashMap(fdr.
                    getdBSequenceHashMap());
            combinedSearchEngines2Mzid.setTandemPeptideEvidenceHashMap(fdr.
                    getPeptideEvidenceMap());
            combinedSearchEngines2Mzid.setTandemPeptideHashMap(fdr.
                    getPeptideHashMap());
            combinedSearchEngines2Mzid.
                    setTandemSpectrumIdentificationItemHashMap(fdr.
                            getSpectrumIdentificationItemHashMap());
//            combinedSearchEngines2Mzid.setTandemSpectrumIdentificationResultHashMap(fdr.getSpectrumIdentificationResultHashMap());
            combinedSearchEngines2Mzid.getUnimodHashmap().putAll(fdr.
                    getUnimodHashmap());
        }

        fdr.computeFDRusingJonesMethod();
        //m.writeTheSortedDataToFile("tmp/debug_singleSE_"+searchEngine+".txt");

        //m.readMzIdentMLData(xmlToRead, searchEngine);
        //m.computeFDRusingJonesMethod();
        Map<String, List<List<String>>> pepMod = new HashMap<>(fdr.
                getFromXMLPeptideModificationHash());
        Map<String, String> pepSeq = new HashMap<>(fdr.
                getFromXMLPeptideSequenceHash());
        Map<String, List<List<Object>>> specInfo = new HashMap<>(fdr.
                getFromXMLSpectrumInfoHash());
        List<String> sorted_spec = new ArrayList<>(fdr.
                getSorted_spectrumResult());
        List<String> sorted_pepID
                = new ArrayList<>(fdr.getSorted_peptideNames());
        List<Double> sorted_evalues = new ArrayList<>(fdr.getSorted_evalues());
        List<Double> sorted_scores = new ArrayList<>(fdr.getSorted_scores());
        List<String> sorted_decoy = new ArrayList<>(fdr.getSorted_decoyOrNot());
        List<Double> sortedFDR = new ArrayList<>(fdr.getSorted_simpleFDR());
        List<Double> sorted_qValues = new ArrayList<>(fdr.getSorted_qValues());
        List<Double> sorted_estFDR = new ArrayList<>(fdr.
                getSorted_estimatedFDR());

        // clear the data, we have all we need
        fdr.clearAllData();

        // Populate the FdrAndMzIdentInformationContainer object
        fdrObj.populateData(xmlToRead, searchEngine, pepMod, pepSeq,
                            specInfo, sorted_spec, sorted_pepID, sorted_evalues,
                            sorted_scores, sorted_decoy, sortedFDR,
                            sorted_qValues, sorted_estFDR);
        // Here the codes only takes AnalysisProtocolCollection result from the first search engine result. 
        if (i == 0) {
            analysisProtocolCollectionTandem = fdr.
                    getAnalysisProtocolCollection();
        }
    }


    /*
     * // ** Overloaded method to deal with Mascot file
     *
     * public void computeFDRForSingleSearchEngine(String xmlToRead, String
     * searchEngine, FdrAndMzIdentInformationContainer fdrObj,String
     * mascotDecoyIdentifier,int decoyRatio) throws Exception{
     *
     * FalseDiscoveryRate fdr = new FalseDiscoveryRate(decoyRatio);
     * fdr.readMzIdentMLData(xmlToRead, searchEngine,mascotDecoyIdentifier);
     * fdr.computeFDRusingJonesMethod();
     *
     * HashMap <String, ArrayList<ArrayList<String>>> pepMod = (HashMap<String,
     * ArrayList<ArrayList<String>>>)
     * fdr.getFromXMLPeptideModificationHash().clone(); HashMap <String, String>
     * pepSeq = (HashMap<String, String>)
     * fdr.getFromXMLPeptideSequenceHash().clone(); HashMap <String,
     * ArrayList<ArrayList<Object>>> specInfo = (HashMap<String,
     * ArrayList<ArrayList<Object>>>) fdr.getFromXMLSpectrumInfoHash().clone();
     * ArrayList <String> sorted_spec = (ArrayList<String>)
     * fdr.getSorted_spectrumResult().clone(); ArrayList <String> sorted_pepID =
     * (ArrayList<String>) fdr.getSorted_peptideNames().clone(); ArrayList
     * <Double> sorted_evalues = (ArrayList<Double>)
     * fdr.getSorted_evalues().clone(); ArrayList <Double> sorted_scores =
     * (ArrayList<Double>) fdr.getSorted_scores().clone(); ArrayList <String>
     * sorted_decoy = (ArrayList<String>) fdr.getSorted_decoyOrNot().clone();
     * ArrayList <Double> sortedFDR = (ArrayList<Double>)
     * fdr.getSorted_simpleFDR().clone(); ArrayList <Double> sorted_qValues =
     * (ArrayList<Double>) fdr.getSorted_qValues().clone(); ArrayList <Double>
     * sorted_estFDR = (ArrayList<Double>) fdr.getSorted_estimatedFDR().clone();
     *
     * // clear the data, we have all we need fdr.clearAllData();
     *
     * // Populate the FdrAndMzIdentInformationContainer object
     * fdrObj.populateData(xmlToRead, searchEngine, pepMod, pepSeq, specInfo,
     * sorted_spec, sorted_pepID, sorted_evalues, sorted_scores, sorted_decoy,
     * sortedFDR, sorted_qValues, sorted_estFDR); }
     */
    /**
     * Do the container assignment with AFS values for each spectrum + sequence
     * pair
     *
     * @param rank Minimum peptide rank for combination.
     *
     * @throws java.lang.Exception Exception
     */
    public void combinePeptidesAcrossSearchEngines(int rank)
            throws Exception {

        // RK 19-02-13
        // Extract ALL spectrum IDs found across all the search engines
        String[] allspectrumIDs = findAllSpectrumIdsFromSearchEngines();
        for (int i = 0; i < allspectrumIDs.length; i++) {
            // debug
            //if(commonSpectrumIDs[i].toString().equalsIgnoreCase("43.3409.3412.3.dta")){
            //	System.out.println("Processing -- " + commonSpectrumIDs[i].toString());
            //}

            // collect the peptides above a given rank and the corresponding sequences for a given spectrum ID	
            String[][] peptideIds = collectPeptideIdentifiersForGivenSpectrumID(
                    allspectrumIDs[i], rank);
            // for each SE, get the peptide sequences
            String[][] peptideSeqs = collectPeptideSequences(peptideIds);
            // For each sequence, identify the SEs which report this sequence
            Map<String, List<Integer>> sequenceSearchEngineMapping
                    = compareSequences(peptideSeqs);
            // Create a 'mo', 'ot' type string representation for the above information
            Map<String, String> seqAndMultipleSeIdentifier
                    = createMultipleSeIdentifier(sequenceSearchEngineMapping);
            // extract the fdr score etc information from peptideID
            Map<String, List<List<Object>>> fdrRelatedInfo
                    = extractFdrRelatedInformationForSeq(
                            sequenceSearchEngineMapping, peptideSeqs, peptideIds,
                            allspectrumIDs[i]);
            // Compute AFS score
            Map<String, Double> pepSeqAndAFS_score = computeAFS_score(
                    fdrRelatedInfo);

            addInformationToCombinedResultContainer(allspectrumIDs[i],
                                                    seqAndMultipleSeIdentifier,
                                                    fdrRelatedInfo,
                                                    pepSeqAndAFS_score);
        }

    }

    /**
     * Update AnalysisProtocolCollection if write out mzid version 1.2.
     * Adding CV term "peptide-level scoring" to thte AdditionalSearchParams.
     * Removing CV term "no special processing" if exists.
     */
    private void handleAnalysisProtocolCollection() {
        List<SpectrumIdentificationProtocol> sipList
                = this.analysisProtocolCollection.
                getSpectrumIdentificationProtocol();
        for (SpectrumIdentificationProtocol sip : sipList) {
            List<CvParam> sipCvParamList = sip.getAdditionalSearchParams().
                    getCvParam();
            boolean flag = false; //flag of existance of the CV term "peptide-level scoring"
            for (CvParam cp : sipCvParamList) {
                if (cp.getAccession().equals(CvConstants.PEPTIDE_LEVEL_SCORING.
                        getAccession())) {
                    flag = true;
                    break;
                }
            }
            if (!flag) {
                sipCvParamList.remove(CvConstants.NO_SPECIAL_PROCESSING);
                sipCvParamList.add(CvConstants.PEPTIDE_LEVEL_SCORING);
            }
        }
    }

    String[] findAllSpectrumIdsFromSearchEngines()
            throws Exception {
        Map<String, String> tempMap = new HashMap<>();

        for (int i = 0; i < singleFDRInformation.length; i++) {
            String[] specNames = singleFDRInformation[i].spectrumInfo.keySet().
                    toArray(new String[0]);

            for (String specID : specNames) {
                if (!tempMap.containsKey(specID)) {
                    tempMap.put(specID, "");
                }
            }
        }

        return tempMap.keySet().toArray(new String[0]).clone();
    }

    void addInformationToCombinedResultContainer(String spectrumId,
                                                 Map<String, String> seqAndMultipleSe,
                                                 Map<String, List<List<Object>>> fdrInfo,
                                                 Map<String, Double> pepSeqAndAFS) {

        Iterator<String> seqSet = seqAndMultipleSe.keySet().iterator();
        while (seqSet.hasNext()) {
            String seq = seqSet.next();
            String containerKey = seqAndMultipleSe.get(seq);
            List<List<Object>> fdrAndOtherInfo = fdrInfo.get(seq);

            //ArrayList<String> infoForThisSeq = new ArrayList<String>();
            List<Object> infoForThisSeq = new ArrayList<>();

            infoForThisSeq.add(spectrumId);
            for (int i = 0; i < fdrAndOtherInfo.size(); i++) //infoForThisSeq.add(fdrAndOtherInfo.elementAt(i).toString());
            {
                infoForThisSeq.add(fdrAndOtherInfo.get(i));
            }

            infoForThisSeq.add(pepSeqAndAFS.get(seq));

            // Finally, put everything into the global container
            //ArrayList<ArrayList<String>> storage = combinedResultContainer.get(containerKey);
            List<List<Object>> storage = combinedResultContainer.get(
                    containerKey);

            storage.add(infoForThisSeq);
            combinedResultContainer.put(containerKey, storage);
        }
    }

    Map<String, Double> computeAFS_score(
            Map<String, List<List<Object>>> fdrRelatedInfo) {

        Map<String, Double> afs_score_hash = new HashMap<>();

        Iterator<String> seqSet = fdrRelatedInfo.keySet().iterator();
        while (seqSet.hasNext()) {
            String seq = seqSet.next();
            List<List<Object>> fdrInfo = fdrRelatedInfo.get(seq);
            int totalVectorsInsideThisVector = fdrInfo.size();
            double AFS_score = 1;

            for (int i = 0; i < totalVectorsInsideThisVector; i++) {
                List<Object> insideVec = fdrInfo.get(i);
                AFS_score = AFS_score * (double) insideVec.get(2);
            }
            // Take the square/cube root according to the no. of SEs for this sequence
            AFS_score = Math.pow(AFS_score, 1.0 / totalVectorsInsideThisVector);

            afs_score_hash.put(seq, AFS_score);
        }

        return afs_score_hash;
    }

    Map<String, List<List<Object>>> extractFdrRelatedInformationForSeq(
            Map<String, List<Integer>> sequenceSearchEngineMapping,
            String[][] peptideSeqs, String[][] peptideIds, String spectrumId) {

        Map<String, List<List<Object>>> pepSeqAndFdrDecoyInfo = new HashMap<>();
        Iterator<String> it = sequenceSearchEngineMapping.keySet().iterator();

        while (it.hasNext()) {
            String seq = it.next();
            List<Integer> seMap = sequenceSearchEngineMapping.get(seq);

            List<List<Object>> infoForThisSeq = new ArrayList<>(1);

            for (int i = 0; i < seMap.size(); i++) {
                int seIndex = seMap.get(i);
                // get all the seq from that SE to check for the corresponding index

                List<String> seqValuesForThisSe = new ArrayList<>(Arrays.asList(
                        peptideSeqs[seIndex]));
                // the position of the seq within that array
                int j = seqValuesForThisSe.indexOf(seq);
                // the peptide ID corresponding to the sequence
                //String correspondingPeptideId = peptideIds[i][j];
                String correspondingPeptideId = peptideIds[seIndex][j];

                // We know the SE and the peptide ID. We can extract rest of the information form
                //singleFDRInformation structure
                //Vector allThePeptideIdFromThisSe = new Vector(Arrays.asList(singleFDRInformation[seIndex].sorted_peptideID));
                Vector<String> allThePeptideIdFromThisSe = new Vector<>(
                        singleFDRInformation[seIndex].sorted_peptideID);

                //RK 19-02-13
                // There might be repeated peptides for different Spec IDs.
                int indexOfThisPeptideId = -1;
                int start = 0;
                boolean found = false;
                while (start < allThePeptideIdFromThisSe.size() && !found) {
                    int idx = 0;
                    try {
                        //      System.out.println("seindex "+seIndex +"  correspondingPeptideId "+ correspondingPeptideId);
                        idx = allThePeptideIdFromThisSe.indexOf(
                                correspondingPeptideId, start);

                        String spec
                                = singleFDRInformation[seIndex].sorted_spectrumID.
                                get(idx);
                        if (spec.equals(spectrumId)) {
                            found = true;
                            indexOfThisPeptideId = idx;
                        }
                        start = idx + 1;
                    } catch (Exception e) {

                        System.out.println(seIndex + "\t"
                                + correspondingPeptideId);

                        e.printStackTrace();
                    }

                }
                if (!found) {
                    System.out.println("No peptide found for " + spectrumId);
                }

                //int indexOfThisPeptideId = allThePeptideIdFromThisSe.indexOf(correspondingPeptideId);
                double fdrEstimated
                        = singleFDRInformation[seIndex].sorted_estimatedFDR.get(
                                indexOfThisPeptideId);
                double qValue = singleFDRInformation[seIndex].sorted_qValues.
                        get(indexOfThisPeptideId);
                String decoyOrNot
                        = singleFDRInformation[seIndex].sorted_decoyornot.get(
                                indexOfThisPeptideId);

                List<Object> fdrRelatedInfo = new ArrayList<>(1);
                fdrRelatedInfo.add(seIndex);				// The identifier SE
                fdrRelatedInfo.add(correspondingPeptideId); // the peptide ID

                fdrRelatedInfo.add(fdrEstimated);

                // Estimated FDR
                fdrRelatedInfo.add(qValue);                 // Estimated q value
                fdrRelatedInfo.add(decoyOrNot);             // decoy or not

                infoForThisSeq.add(fdrRelatedInfo);
            }

            pepSeqAndFdrDecoyInfo.put(seq, infoForThisSeq);
        }
        return pepSeqAndFdrDecoyInfo;
    }

    Map<String, String> createMultipleSeIdentifier(
            Map<String, List<Integer>> sequenceSearchEngineMap) {

        Map<String, String> seqAndMappedSe = new HashMap<>();

        Iterator<String> sequencesIterator = sequenceSearchEngineMap.keySet().
                iterator();
        while (sequencesIterator.hasNext()) {
            String seCombination = "";
            String sequence = sequencesIterator.next();
            List<Integer> searchEnginesFound = sequenceSearchEngineMap.get(
                    sequence);
            int totalSe = searchEnginesFound.size();

            // build identifier string
            for (int k = 0; k < totalSe; k++) {
                int index = searchEnginesFound.get(k);
                seCombination = seCombination.concat(
                        shortRepresentationOfSearchEngineIdentiifers[index]);
            }
            // certain combinations need to be corrected to be recognised by the container
            if (seCombination.equalsIgnoreCase("21")) {
                seCombination = "12";
            }
            if (seCombination.equalsIgnoreCase("32")) {
                seCombination = "23";
            }
            if (seCombination.equalsIgnoreCase("31")) {
                seCombination = "13";
            }
            if (seCombination.equalsIgnoreCase("132") || seCombination.
                    equalsIgnoreCase("213")
                    || seCombination.equalsIgnoreCase("231") || seCombination.
                    equalsIgnoreCase("312")
                    || seCombination.equalsIgnoreCase("321")) {
                seCombination = "123";
            }

            seqAndMappedSe.put(sequence, seCombination);
        }
        return seqAndMappedSe;
    }

    String[][] collectPeptideIdentifiersForGivenSpectrumID(String specID,
                                                           int rank) {

        String[][] peptideIds = new String[noOfSearchEngines][];

        try {
            for (int i = 0; i < noOfSearchEngines; i++) {
                //System.out.println("::Debug - " + searchEngineIdentifiers[i] + "\t" + specID);
                List<List<Object>> dataForThisSpectrumID
                        = singleFDRInformation[i].spectrumInfo.get(specID);

                // RK - 19-02-13, if this spec is not identified by the SE, put a null
                if (dataForThisSpectrumID == null) {
                    peptideIds[i] = new String[1];
                    peptideIds[i] = null;
                    continue;
                }

                List<String> pepId = new ArrayList<>(1);
                for (int j = 0; j < dataForThisSpectrumID.size(); j++) {
                    //System.out.println("::Debug 2 - " + dataForThisSpectrumID.get(j).get(1).toString() + "\t" +dataForThisSpectrumID.get(j).get(2).toString() + "\t" + dataForThisSpectrumID.get(j).get(3).toString());

                    if (Integer.parseInt(dataForThisSpectrumID.get(j).get(2).
                            toString()) <= rank) {
                        pepId.add(
                                dataForThisSpectrumID.get(j).get(1).toString().
                                trim());
                    }
                }
                peptideIds[i] = new String[pepId.size()];
                peptideIds[i] = pepId.toArray(new String[0]);
            }
        } catch (NullPointerException n) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + n.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 06 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

        return peptideIds.clone();
    }

    String[][] collectPeptideSequences(String[][] peptideId) {
        String[][] peptideSequences = new String[noOfSearchEngines][];
        try {

            for (int i = 0; i < peptideId.length; i++) {

                // RK - 19-02-13, if this peptide is not present, because there was no corresponding specID, put a null
                if (peptideId[i] == null) {
                    peptideSequences[i] = new String[1];
                    peptideSequences[i] = null;
                    continue;
                }

                List<String> seqs = new ArrayList<>(1);
                for (int j = 0; j < peptideId[i].length; j++) {
                    seqs.add(singleFDRInformation[i].peptideSequence.get(
                            peptideId[i][j]).trim());
                }
                peptideSequences[i] = new String[seqs.size()];
                peptideSequences[i] = seqs.toArray(new String[0]);
            }

        } catch (NullPointerException n) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + n.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 06 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        return peptideSequences.clone();
    }

    Map<String, List<Integer>> compareSequences(String[][] peptideSeqs) {

        Map<String, List<Integer>> sequenceMap = new HashMap<>();
        // extract the unique sequences
        for (int i = 0; i < peptideSeqs.length; i++) {

            // RK 19-02-13, skip the null string
            if (peptideSeqs[i] == null) {
                continue;
            }

            for (int j = 0; j < peptideSeqs[i].length; j++) {
                sequenceMap.put(peptideSeqs[i][j], new ArrayList<Integer>(1));
            }
        }
        // check which sequence lies in which search engine
        Iterator<String> seqIterator = sequenceMap.keySet().iterator();
        while (seqIterator.hasNext()) {
            String pepSeq = seqIterator.next();
            for (int i = 0; i < peptideSeqs.length; i++) {

                // RK 19-02-13, skip the null string
                if (peptideSeqs[i] == null) {
                    continue;
                }

                List<String> v1 = new ArrayList<>(Arrays.asList(peptideSeqs[i]));
                if (v1.contains(pepSeq)) {
                    List<Integer> vtemp = sequenceMap.get(pepSeq);
                    vtemp.add(i);
                    sequenceMap.put(pepSeq, vtemp);
                }
            }
        }
        // Added by Fawaz Ghali 30/7/2015 Debug info
//        Set keys = sequenceMap.keySet();
//        for (Iterator i = keys.iterator(); i.hasNext();) {
//            String key = (String) i.next();
//            List value = (List) sequenceMap.get(key);
//            String s="";
//            for (int j = 0; j < value.size(); j++) {
//                Object object = value.get(j);
//                s= s+"; "+object;
//                
//            }
//            System.out.println( key + ", " + s);
//        }

        return sequenceMap;
    }

    public void writeToFileForDiagnostics(String fileName)
            throws Exception {

        Writer out = new BufferedWriter(new FileWriter(fileName));

        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            String key = it.next();
            out.write(" [ " + key + " ] \n");
            //	ArrayList<ArrayList<String>> info = combinedResultContainer.get(key);
            List<List<Object>> info = combinedResultContainer.get(key);

            for (int i = 0; i < info.size(); i++) {
                List<Object> insideOne = info.get(i);
                out.write(insideOne.toString() + "\n");
            }
        }
        out.close();
    }

    public void writeToFileForDiagnostics_singleFDRObj(String fileName)
            throws Exception {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        // reassign the line number for a new file
        this.lineNumber = 1;

        for (int i = 0; i < singleFDRInformation.length; i++) {
            out.write(searchEngineIdentifiers[i] + "\n-----------------\n");

            FdrAndMzIdentInformationContainer fdrInfo = singleFDRInformation[i];
            List<String> peptides = fdrInfo.sorted_peptideID;

            for (int k = 0; k < peptides.size(); k++) {
                //String seIndex = String.valueOf(i);
                String peptideId = peptides.get(k);
                String specId = fdrInfo.sorted_spectrumID.get(k);
                String simpleFDR = fdrInfo.sortedSimpleFDR.get(k).toString();
                String estFDR = fdrInfo.sorted_estimatedFDR.get(k).toString();

                String csvDescription = retrieveInformationToWriteInCsvFile(i,
                                                                            peptideId,
                                                                            specId,
                                                                            simpleFDR,
                                                                            estFDR,
                                                                            "");
                out.write(csvDescription);
            }

        }
        out.close();
    }

    void computeSortedIndicesForSingleContainer(String key) {

        //ArrayList<ArrayList<String>> content = (ArrayList<ArrayList<String>>) combinedResultContainer.get(key);
        List<List<Object>> content
                = (List<List<Object>>) combinedResultContainer.get(key);

        List<Double> afs_values = new ArrayList<>(1);
        for (int i = 0; i < content.size(); i++) //afs_values.add( Double.parseDouble(content.get(i).get(content.get(i).size() - 1)));
        {
            afs_values.add(Double.parseDouble(content.get(i).get(content.get(i).
                    size() - 1).toString()));
        }

        // Call the sorting routine to find the indices of sorted evalues
        TreeSortForIndices sortClass = new TreeSortForIndices();
        Integer[] sortOrderForEvalues = sortClass.sortTheValueColumn(afs_values.
                toArray(new Double[0]), true);

        List<List<Object>> sort_content = new ArrayList<>();

        for (int i = 0; i < sortOrderForEvalues.length; i++) {
            //ArrayList<String> tempRow = content.get(sortOrderForEvalues[i]);
            List<Object> tempRow = content.get(sortOrderForEvalues[i]);
            sort_content.add(new ArrayList<>(tempRow));
        }

        combinedResultContainer.put(key, new ArrayList<>(sort_content));
    }

    public void sortWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            computeSortedIndicesForSingleContainer(it.next());
        }
    }

    public void insertFakeDecoyInTheEndInSingleContainer(String key) {

        Double epsilon_afs = 0.1d;

        List<List<Object>> content = combinedResultContainer.get(key);

        List<Object> lastRec = content.get(content.size() - 1);
        Double lastAfs = Double.parseDouble(lastRec.get(lastRec.size() - 1).
                toString());
        Double fakeAfs = new Double(lastAfs + epsilon_afs);

        //create a fake data-structure
        List<Object> fake_ds = new ArrayList<>();
        fake_ds.add("Fake_spectrum");

        for (int i = 0; i < key.length(); i++) {
            List<Object> fake_vec = new ArrayList<>();
            fake_vec.add(i);
            fake_vec.add("Fake Peptide");
            fake_vec.add(epsilon_afs);
            fake_vec.add(epsilon_afs);
            fake_vec.add(true);

            fake_ds.add(fake_vec);
        }

        fake_ds.add(fakeAfs);

        content.add(fake_ds);
        combinedResultContainer.put(key, new ArrayList<>(content));
    }

    public void insertFakeDecoyInWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            String key = it.next();
            if (!combinedResultContainer.get(key).isEmpty()) {
                insertFakeDecoyInTheEndInSingleContainer(key);
            }
        }
    }

    public void removeFakeDecoyFromTheEndInSingleContainer(String key) {
        List<List<Object>> content = combinedResultContainer.get(key);
        combinedResultContainer.get(key).remove(content.size() - 1);
    }

    /**
     * Remove fake decoy from all the containers
     */
    public void removeFakeDecoyFromWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            String key = it.next();
            if (!combinedResultContainer.get(key).isEmpty()) {
                removeFakeDecoyFromTheEndInSingleContainer(key);
            }
        }
    }

    void computeSimpleFDRForSingleContainer(String key) {
        int falsePositiveCount = 0;

        for (int i = 0; i < combinedResultContainer.get(key).size(); i++) {
            List<Object> vec = (List<Object>) combinedResultContainer.get(key).
                    get(i).get(1);
            if (vec.get(4).toString().equals("true")) {
                falsePositiveCount++;
            }

            combinedResultContainer.get(key).get(i).add(
                    (double) falsePositiveCount / (double) (i + 1));
            combinedResultContainer.get(key).get(i).add((double) 0);
        }
    }

    public void simpleFdrForWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            computeSimpleFDRForSingleContainer(it.next());
        }
    }

    private void computeQvalueForSingleContainer(String key) {

        if (combinedResultContainer.get(key).isEmpty()) {
            return;
        }

        int indexForFdrValue = combinedResultContainer.get(key).get(0).size()
                - 2;
        int indexForQValue = combinedResultContainer.get(key).get(0).size() - 1;

        Double immediateMinFdr = 0.0;
        //try{
        immediateMinFdr = (Double) combinedResultContainer.get(key).get(
                combinedResultContainer.get(key).size() - 1).get(
                indexForFdrValue);
        combinedResultContainer.get(key).get(combinedResultContainer.get(key).
                size() - 1).set(indexForQValue, immediateMinFdr);
        //}catch(Exception e){
        //	e.printStackTrace();
        //}
        for (int i = combinedResultContainer.get(key).size() - 1; i > 0; i--) {
            Double currentFDR = (Double) combinedResultContainer.get(key).get(i
                    - 1).get(indexForFdrValue);

            if (currentFDR < immediateMinFdr) {
                immediateMinFdr = currentFDR;
            }

            combinedResultContainer.get(key).get(i - 1).set(indexForQValue,
                                                            immediateMinFdr);
        }
    }

    /**
     * Compute q value for all the containers
     */
    public void qValueForWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            computeQvalueForSingleContainer(it.next());
        }
    }

    /**
     * Compute FDR Score
     *
     */
    private void computeFDRScoreForSingleContainer(String key) {

        if (combinedResultContainer.get(key).isEmpty()) {
            return;
        }

        int indexForAFSValue = combinedResultContainer.get(key).get(0).size()
                - 3;
        int indexForFdrValue = combinedResultContainer.get(key).get(0).size()
                - 2;
        int indexForQValue = combinedResultContainer.get(key).get(0).size() - 1;

//		for (int i = 0 ; i < combinedResultContainer.get(key).size(); i++){
//			combinedResultContainer.get(key).get(i).add((Double)0);
//		}
        /*
         * RK -22-10-10 double prev_afsvalue =
         * (Double)combinedResultContainer.get(key).get(0).get(indexForAFSValue);
         * double prev_qvalue =
         * (Double)combinedResultContainer.get(key).get(0).get(indexForQValue);
         * double prev_prev_afsvalue = prev_afsvalue; // previous to previous
         * evalue in case of straight vertical rise in q value without any
         * change in evalue
         */
        double prev_afsvalue = 0d; 			// RK -22-10-10
        double prev_qvalue = 0d;   			// RK -22-10-10
        double prev_prev_afsvalue = 0d; 	// RK -22-10-10

        int counter_backwardStep = 0;
        //int i = 1 ;
        int i = 0; // RK 22-10-10
        for (; i < combinedResultContainer.get(key).size(); i++) {
            double current_afsvalue = (Double) combinedResultContainer.get(key).
                    get(i).get(indexForAFSValue);
            double current_qvalue = (Double) combinedResultContainer.get(key).
                    get(i).get(indexForQValue);

            if (current_qvalue > prev_qvalue) {

                double slope;
                double intercept;

                int id = 0;
                if (current_afsvalue != prev_afsvalue) {
                    slope = (current_qvalue - prev_qvalue) / (current_afsvalue
                            - prev_afsvalue);
                    id = 1; //RK 22-10-10
                } else {
                    slope = (current_qvalue - prev_qvalue) / (current_afsvalue
                            - prev_prev_afsvalue);
                    id = 2; //RK 22-10-10
                }
                //RK 22-10-10 - Tells us which co-ordinates to use for calculating intercepts.
                if (id == 1) {
                    intercept = prev_qvalue - slope * prev_afsvalue;
                } else {
                    intercept = prev_qvalue - slope * prev_prev_afsvalue;
                }

                if (counter_backwardStep > 0) { // compute the FDR score for flat q-value region
                    for (int k = 0; k <= counter_backwardStep; k++) {
                        int index = i - counter_backwardStep + k;
                        double fdrScore = slope
                                * ((Double) combinedResultContainer.get(key).
                                get(index).get(indexForAFSValue)) + intercept;
                        if (Double.isNaN(fdrScore)) {
                            System.out.print(fdrScore);
                        }
                        combinedResultContainer.get(key).get(index).
                                add(fdrScore);
                        //System.out.println("i = " + i + "Count = " + counter_backwardStep +  " k = " + k + " index::1 = " + index + " slope = " + slope + " intercept = " + intercept + " e-val = " + ((Double)combinedResultContainer.get(key).get(index).get(indexForAFSValue)) + " FDR = " + fdrScore);
                    }
                } else { 							// In case an immediate increment in q value is found
                    double fdrScore = slope * current_afsvalue + intercept;
                    combinedResultContainer.get(key).get(i).add(fdrScore);
                    //System.out.println("index::2 i = " + i + " FDR = " + fdrScore);
                }

                // Re-initialise the variables
                counter_backwardStep = 0;
                if (current_afsvalue > prev_afsvalue) { // the previous e-value will change only if the current e-value is different
                    prev_prev_afsvalue = prev_afsvalue;
                    prev_afsvalue = current_afsvalue;
                }
                prev_qvalue = current_qvalue;

            } else {
                counter_backwardStep++;
            }
        }

    }

    public void estFDRForWholeCombinedResultContainer() {
        Iterator<String> it = combinedResultContainer.keySet().iterator();
        while (it.hasNext()) {
            computeFDRScoreForSingleContainer(it.next());
        }
    }

    public void prepareTheCSVFileForMzIdentMLParser(String fileName)
            throws Exception {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        // reassign the line number for a new file
        this.lineNumber = 1;

        Iterator<String> seComboKeys = combinedResultContainer.keySet().
                iterator();

        while (seComboKeys.hasNext()) {
            String key = seComboKeys.next();

            out.write("\n\n Sequences in the container : [" + key + "]\n\n");

            for (int i = 0; i < combinedResultContainer.get(key).size(); i++) {

                String specID = combinedResultContainer.get(key).get(i).get(0).
                        toString();
                Double estFDR_value = (Double) combinedResultContainer.get(key).
                        get(i).get(combinedResultContainer.get(key).get(i).
                        size() - 1);
                //String estFDR = combinedResultContainer.get(key).get(i).get(combinedResultContainer.get(key).get(i).combinedResultContainer.get(key).get(i).get(combinedResultContainer.get(key).get(i).size() - 1)size() - 1).toString();
                String estFDR = Double.toString(estFDR_value);
                if (Double.valueOf(estFDR).isNaN()) {
                    System.out.println(estFDR);
                }

                int noOfSearchEngineInThisContainer = key.length();
                for (int k = 1; k < noOfSearchEngineInThisContainer + 1; k++) { // k < noOfSearchEngines is not right !
                    //this.lineNumber++;
                    //System.out.println(key + "  " + specID + "  " + estFDR);
                    //System.out.println(i + " <-> " + k  );
                    //System.out.println(	" -- " + ((Vector)combinedResultContainer.get(key).get(i).get(k)).toString());

                    int seIndex
                            = (Integer) ((List<Object>) combinedResultContainer.
                            get(key).get(i).get(k)).get(0);

                    String pepId = ((List<Object>) combinedResultContainer.get(
                            key).get(i).get(k)).get(1).toString();
                    String simpleFDR = ((List<Object>) combinedResultContainer.
                            get(key).get(i).get(k)).get(2).toString();
                    String toWrite
                            = retrieveInformationToWriteInCsvFile(seIndex, pepId,
                                                                  specID,
                                                                  simpleFDR,
                                                                  estFDR, key);

                    out.write(toWrite);
                }
            }
        } //end of while
        out.close();// Updated by Fawaz Ghali 04/03/2014 closing BufferedWriter before writing MzidFile
        createMetadata();
        writeMzidFile(); // Commented by RK

    }

    /**
     * combinedResultContainer and singleFDRInformation should be queried for
     * the information needed to prepare data for CSV file. Each line in CSV
     * file should have the following in the mentioned order -
     *
     * lineCounter, dtaFileName, seqString, Estimated fdr
     * score,calculatedMass,groupID, protAcc,
     * startIndex,$endIndex,peptideSpecificFDR,
     * modString,chargeState,experimentalMass "\n"
     *
     * @return the comma separated string with the above information
     */
    String retrieveInformationToWriteInCsvFile(int searchEngineIndex,
                                               String pepId, String spectrumID,
                                               String simpleFDR, String estFDR,
                                               String key)
            throws Exception {

        String outString = new String();

        //System.out.println(pepId + "\t " + singleFDRInformation[searchEngineIndex].peptideModifications.containsKey(pepId));
        List<List<String>> modArray
                = singleFDRInformation[searchEngineIndex].peptideModifications.
                get(pepId);
        String modStr = createModString(modArray);

        String pepSequence
                = singleFDRInformation[searchEngineIndex].peptideSequence.get(
                        pepId);

        List<List<Object>> specInfo
                = singleFDRInformation[searchEngineIndex].spectrumInfo.get(
                        spectrumID);

        List<Object> relevantPepInfo = new ArrayList<>();
        for (int i = 0; i < specInfo.size(); i++) {
            if (specInfo.get(i).get(1).equals(pepId)) {
                relevantPepInfo = specInfo.get(i);
                break;
            }
        }

        if (relevantPepInfo.size() < 5) {
            System.out.println("Wrong peptide info : " + searchEngineIndex + " "
                    + pepId + " " + spectrumID);
            //System.in.read();
            return "\n ";
        }

        String charge = relevantPepInfo.get(3).toString();
        String calcMass = relevantPepInfo.get(4).toString();
        String expMAss = relevantPepInfo.get(5).toString();

        //FG
        List<List<Object>> peptideEvd = new ArrayList<>();
        List<String> improvedPSM = new ArrayList<>();

        if (relevantPepInfo.size() > 8) {
            peptideEvd = (List<List<Object>>) relevantPepInfo.get(8);
            for (int i = 0; i < peptideEvd.size(); i++) {
                outString = outString
                        + this.lineNumber + "," + spectrumID + "," + pepSequence
                        + ","
                        + estFDR + "," + calcMass + "," + "GROUP_ID" + ","
                        + peptideEvd.get(i).get(3).toString()
                        + "," + peptideEvd.get(i).get(0).toString() + ","
                        + peptideEvd.get(i).get(1).toString() + "," + simpleFDR
                        + "," + modStr
                        + "," + charge + "," + expMAss + "\n";
                String proteinAccession = peptideEvd.get(i).get(3).toString();
                int start = Integer.
                        parseInt(peptideEvd.get(i).get(0).toString());
                int end = Integer.parseInt(peptideEvd.get(i).get(1).toString());
                // Added by FG 09/5/2014 Missing pre post
                String pre = peptideEvd.get(i).get(4).toString();
                String post = peptideEvd.get(i).get(5).toString();
                //Added by FG 13/8/2014 fixing multiple spectra ref
                String location = peptideEvd.get(i).get(6).toString();
                String improvedPSMString = spectrumID + "," + pepSequence + ","
                        + peptideEvd.get(i).get(0).toString() + ","
                        + peptideEvd.get(i).get(1).toString() + ","
                        + proteinAccession;

                if (!improvedPSM.contains(improvedPSMString) && !key.equals("")) {

                    createPSM(spectrumID, proteinAccession, pepSequence, modStr,
                              modArray, start, end, charge, expMAss, calcMass,
                              estFDR, key, pre, post, location);
                    improvedPSM.add(improvedPSMString);
                }

                this.lineNumber++;
            }
        }

        return outString;
    }

    public void createMetadata() {

        analysisSoftwareList = fdr.getAnalysisSoftwareList();
        auditCollection = fdr.getAuditCollection();
        provider = fdr.getProvider();
        analysisProtocolCollection = analysisProtocolCollectionTandem;
        cvList = fdr.getCvList();
        analysisCollection = fdr.getAnalysisCollection();
        //analysisCollection.getSpectrumIdentification().get(0).getInputSpectra().clear();
        List<SpectrumIdentification> spectrumIdentificationList
                = analysisCollection.getSpectrumIdentification();
        for (int i = 0; i < spectrumIdentificationList.size(); i++) {
            SpectrumIdentification spectrumIdentification
                    = spectrumIdentificationList.get(i);
            List<InputSpectra> inputSpectraList = spectrumIdentification.
                    getInputSpectra();
            for (int j = 0; j < inputSpectraList.size(); j++) {
                InputSpectra inputSpectra = inputSpectraList.get(j);
                String ref = inputSpectra.getSpectraDataRef();
                String location = spectraIDLocation.get(ref);
                // 
                // 2015-7-20 by Bo Wen: according to the file name to map the spectraData
                File locFile = new File(location);
                location = locFile.getName();

                SpectraData sd = spectraDataHashMap.get(location);
                inputSpectra.setSpectraData(sd);

            }

        }

        inputs = fdr.getInputs();
        inputs.getSpectraData().clear();
        inputs.getSpectraData().addAll(spectraDataHashMap.values());
    }

    public void createPSM(String spectrumID, String proteinAccession,
                          String pepSequence, String modStr,
                          List<List<String>> modArray, int start, int end,
                          String charge, String expMAss, String calcMass,
                          String estFDR, String key, String pre, String post,
                          String location)
            throws IOException {
        //FG

        SpectrumIdentificationResult sir;
        if (combinedSearchEngines2Mzid.
                getCombinedSpectrumIdentificationResultHashMap().containsKey(
                        spectrumID)) {
            sir = combinedSearchEngines2Mzid.
                    getCombinedSpectrumIdentificationResultHashMap().get(
                            spectrumID);
        } else {
            sir = new SpectrumIdentificationResult();
            sir.setId("SIR_" + sirCounter);
            sirCounter = sirCounter + 1;
            sir.setSpectrumID(spectrumID);
            //Added by FG 13/8/2014 fixing multiple spectra ref

            // 2015-7-20 by Bo Wen: according to the file name to map the spectraData
            File locFile = new File(location);
            location = locFile.getName();
            SpectraData sd = spectraDataHashMap.get(location);
            if (sd != null) {
                sir.setSpectraData(sd);
            } else {
                throw new RuntimeException(
                        "Spectra data is missing for location: " + sd);
            }

            combinedSearchEngines2Mzid.
                    getCombinedSpectrumIdentificationResultHashMap().put(
                            spectrumID, sir);
        }

        DBSequence dbSequence = null;
        String newDBsequenceID = proteinAccession;//proteinAccession

        if (combinedSearchEngines2Mzid.getCombinedDbSequenceHashMap().
                containsKey(newDBsequenceID)) {
            dbSequence = combinedSearchEngines2Mzid.
                    getCombinedDbSequenceHashMap().get(newDBsequenceID);
        } else {
            dbSequence = new DBSequence();
            dbSequence.setId("dbseq_" + newDBsequenceID);
            dbSequence.setAccession(newDBsequenceID);
            if (fdr.getInputs().getSearchDatabase().get(0) != null) {
                dbSequence.setSearchDatabase(
                        fdr.getInputs().getSearchDatabase().get(0));
            }
            combinedSearchEngines2Mzid.getCombinedDbSequenceHashMap().put(
                    newDBsequenceID, dbSequence);
        }

        Peptide pep;
        String pepID = pepSequence + "_" + modStr;

        if (combinedSearchEngines2Mzid.getCombinedPeptideHashMap().containsKey(
                pepID)) {
            pep = combinedSearchEngines2Mzid.getCombinedPeptideHashMap().get(
                    pepID);
        } else {
            pep = new Peptide();
            pep.setPeptideSequence(pepSequence);
            pep.setId(pepID);
            combinedSearchEngines2Mzid.getCombinedPeptideHashMap().put(pepID,
                                                                       pep);
            List<Modification> modList = pep.getModification();
            //System.out.println("New pep:" + pepSeq + modString);

            for (int z = 0; z < modArray.size(); z++) {
                Modification mzidMod = new Modification();
                Boolean foundOkay = true;
                mzidMod.setLocation(Integer.valueOf(modArray.get(z).get(0)));
                if (modArray.get(z).get(2) != null) {
                    mzidMod.setMonoisotopicMassDelta(Double.valueOf(modArray.
                            get(z).get(2)));
                } else {
                    foundOkay = false;
                }
                if (modArray.get(z).get(1) != null) {
                    mzidMod.getResidues().add(modArray.get(z).get(1));
                } else {
                    foundOkay = false;
                }
                if (modArray.get(z).get(3) != null && foundOkay) {
                    CvParam modParam = new CvParam();
                    String accession = combinedSearchEngines2Mzid.
                            getUnimodHashmap().get(modArray.get(z).get(3));
                    modParam.setAccession(accession);

                    modParam.setName(modArray.get(z).get(3));
                    Cv unimodCV;
                    unimodCV = new Cv();
                    unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
                    unimodCV.setId("UNIMOD");
                    unimodCV.setFullName("UNIMOD");

                    modParam.setCv(unimodCV);
                    mzidMod.getCvParam().add(modParam);
                } else {
                    foundOkay = false;
                }

                if (!foundOkay) {
                    CvParam modParam = new CvParam();
                    modParam.setAccession("MS:1001460");

                    modParam.setName("unknown modification");
                    Cv psiCV = new Cv();
                    psiCV.setUri(
                            "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
                    psiCV.setId("PSI-MS");
//                    psiCV.setVersion("2.25.0");
                    psiCV.setFullName("PSI-MS");
                    modParam.setCv(psiCV);
                    mzidMod.getCvParam().add(modParam);
                }
                modList.add(mzidMod);
            }

        }

        String pepEvidID = pepSequence + "_" + newDBsequenceID + "_" + start
                + "_" + end;

        PeptideEvidence peptideEvidence = null;
        if (combinedSearchEngines2Mzid.getCombinedPeptideEvidenceHashMap().
                containsKey(pepEvidID)) {
            peptideEvidence = combinedSearchEngines2Mzid.
                    getCombinedPeptideEvidenceHashMap().get(pepEvidID);

        } else {
            peptideEvidence = new PeptideEvidence();
            peptideEvidence.setId(pepEvidID);
            peptideEvidence.setDBSequence(dbSequence);
            peptideEvidence.setPeptide(pep);

            peptideEvidence.setStart(start);
            peptideEvidence.setEnd(end);

            // Added by FG 09/5/2014 Missing pre post
            peptideEvidence.setPre(pre);
            peptideEvidence.setPost(post);

            peptideEvidence.setIsDecoy(Boolean.FALSE);
            if (decoyRegex != null) {
                if (newDBsequenceID.contains(decoyRegex)) {
                    peptideEvidence.setIsDecoy(Boolean.TRUE);
                }
            }

            combinedSearchEngines2Mzid.getCombinedPeptideEvidenceHashMap().put(
                    pepEvidID, peptideEvidence);

        }

        List<SpectrumIdentificationItem> siiList = sir.
                getSpectrumIdentificationItem();

        /*
         * Cases: If PepID is the same, then this is another PeptideEvidence,
         * use same sii If PepID is different
         */
        SpectrumIdentificationItem sii = null;

        for (SpectrumIdentificationItem currentSii : siiList) {
            String currentPepSeq = currentSii.getPeptide().getPeptideSequence();
            if (currentPepSeq.equals(pepSequence)) {
                sii = currentSii;
                //TODO - how to exit foreach loop in Java?
            }
        }

        if (sii == null) {
            sii = new SpectrumIdentificationItem();

            sii.setChargeState(Integer.parseInt(charge));
            sii.setExperimentalMassToCharge(Double.parseDouble(expMAss));
            sii.setCalculatedMassToCharge(Double.parseDouble(calcMass));

            sii.setPeptide(pep);
            List<CvParam> cvParamList = sii.getCvParam();

            cvParamList.add(combinedSearchEngines2Mzid.makeCvParam("MS:1002356",
                                                                   "PSM-level combined FDRScore",
                                                                   estFDR));

            Boolean orderLowToHigh = true;

            combinedSearchEngines2Mzid.addSIIToListAndSetRank(siiList, sii,
                                                              "MS:1002356",
                                                              orderLowToHigh,
                                                              sir.getId());
            UserParam e = new UserParam();
            e.setName("search engines identifying PSM");
            e.setValue(key);
            sii.getUserParam().add(e);

        }
        PeptideEvidenceRef pepEvidRef = new PeptideEvidenceRef();
        pepEvidRef.setPeptideEvidence(peptideEvidence);
        List<PeptideEvidenceRef> peptideEvidenceRefList = sii.
                getPeptideEvidenceRef();
        boolean peptideEvidenceRefBool = false;
        for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.
                    get(i);

            if (peptideEvidenceRef.getPeptideEvidenceRef().equals(pepEvidRef.
                    getPeptideEvidenceRef())) {
                peptideEvidenceRefBool = true;
                break;
            }
        }
        if (!peptideEvidenceRefBool) {
            sii.getPeptideEvidenceRef().add(pepEvidRef);
        }

        // Remove FDRScore cv param
        List<CvParam> cvParamList = sii.getCvParam();
        for (int i = 0; i < cvParamList.size(); i++) {
            CvParam cvParam = cvParamList.get(i);
            if (cvParam.getAccession().equals("MS:1001874")) {
                cvParamList.remove(cvParam);
            }

        }

    }

    String createModString(List<List<String>> modArray) {
        String modString = new String();

        for (int i = 0; i < modArray.size(); i++) {
            String forThisMod = new String();
            if (modArray.get(i).get(3).equals("unknown modification")) {
                forThisMod = modArray.get(i).get(2) + "_" + modArray.get(i).get(
                        1) + ":" + modArray.get(i).get(0);
            } else {
                forThisMod = modArray.get(i).get(3) + "(" + modArray.get(i).get(
                        1) + "):" + modArray.get(i).get(0);
            }

//            forThisMod = forThisMod.replaceAll("&gt;", "_");
//            forThisMod = forThisMod.replaceAll("&lt;", "_");
//            forThisMod = forThisMod.replaceAll("&quot;", "_");
//            forThisMod = forThisMod.replaceAll("&apos;", "_");
//            forThisMod = forThisMod.replaceAll("&amp;", "_");
            forThisMod = forThisMod.replaceAll("->", "_");

            modString = modString + "##" + forThisMod;
        }

        return modString;
    }

    public static void runTwoSearchEngines(String[] args)
            throws Exception {

        String file_2 = args[0];
        String searchEngine_2 = args[1];
        String secondcvTerm = args[2];
        String secondbetterScoresAreLower = args[3];

        String file_3 = args[4];
        String searchEngine_3 = args[5];
        String thirdcvTerm = args[6];
        String thirdbetterScoresAreLower = args[7];

        int rank = Integer.parseInt(args[8]);
        int decoyRatio = Integer.parseInt(args[9]);

        csvFileName = args[10];
        //csvFileName = csvFileName + ".txt";

        String debugFileName = null;
        if (args[11] != null) {
            debugFileName = args[11];
        }

        String[] searchEngine = {searchEngine_3, searchEngine_2};
        String[] inputFiles = {file_3, file_2};
        String[] cvTerms = {thirdcvTerm, secondcvTerm};
        String[] betterScoresAreLower = {thirdbetterScoresAreLower,
            secondbetterScoresAreLower};

        CombineSearchEngines C = new CombineSearchEngines(searchEngine);
        // Added by FG 8/10/2014
        C.checkInputFiles(inputFiles);
        String decoyRegex = args[12];
        mzidVer = MzIdentMLVersion.getVersion(args[13]);
        long startTime;
        long stopTime;
        long elapsedTime;
        // Perform algorithm 1 and fill all the values in FdrAndMzIdentInformationContainer of C
        for (int i = 0; i < searchEngine.length; i++) {
            startTime = System.currentTimeMillis();
            C.computeFDRForSingleSearchEngine(i, inputFiles[i], searchEngine[i],
                                              C.singleFDRInformation[i],
                                              decoyRatio, decoyRegex, cvTerms[i],
                                              betterScoresAreLower[i]);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("computeFDRForSingleSearchEngine time for "
                        + searchEngine[i] + " " + elapsedTime / 1000
                        + " Seconds");
            }

            System.out.println("Done - " + searchEngine[i]);
        }

        if (args[11] != null) {
            startTime = System.currentTimeMillis();
            C.writeToFileForDiagnostics_singleFDRObj(debugFileName);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println(
                        "writeToFileForDiagnostics_singleFDRObj time "
                        + +elapsedTime / 1000 + " Seconds");
            }

        }
        startTime = System.currentTimeMillis();
        // The common peptides across search engines processing starts here
        C.combinePeptidesAcrossSearchEngines(rank);
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("combinePeptidesAcrossSearchEngines time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        // do the sorting according to afs value
        startTime = System.currentTimeMillis();

        C.sortWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("sortWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();
        C.insertFakeDecoyInWholeCombinedResultContainer(); //RK 15-08-12
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println(
                    "insertFakeDecoyInWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.simpleFdrForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("simpleFdrForWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.qValueForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("qValueForWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.estFDRForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("estFDRForWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }
        startTime = System.currentTimeMillis();

        C.removeFakeDecoyFromWholeCombinedResultContainer(); //RK 15-08-12
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println(
                    "removeFakeDecoyFromWholeCombinedResultContainer time "
                    + +elapsedTime / 1000 + " Seconds");
        }

        // if(args[7] != null)
        //	 C.writeToFileForDiagnostics(debugFileName);
        startTime = System.currentTimeMillis();

        C.prepareTheCSVFileForMzIdentMLParser(csvFileName);
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("prepareTheCSVFileForMzIdentMLParser time "
                    + +elapsedTime / 1000 + " Seconds");
        }
    }

    /**
     * FG The method to call the pipeline from other java programs, in case if
     * we want to tinker around in main()
     *
     * @param args Input arguments
     *
     * @throws java.lang.Exception Exception
     */
    public static void runThreeSearchEngines(String[] args)
            throws Exception {

        /**
         * newArgs[0] = Utils.getCmdParameter(args, "firstFile", true);
         * newArgs[1] = Utils.getCmdParameter(args, "firstSearchEngine", true);
         * newArgs[2] = Utils.getCmdParameter(args, "firstcvTerm", true);
         * newArgs[3] = Utils.getCmdParameter(args, "firstbetterScoresAreLower",
         * true); newArgs[4] = Utils.getCmdParameter(args, "secondFile", true);
         * newArgs[5] = Utils.getCmdParameter(args, "secondSearchEngine", true);
         * newArgs[6] = Utils.getCmdParameter(args, "secondcvTerm", true);
         * newArgs[7] = Utils.getCmdParameter(args,
         * "secondbetterScoresAreLower", true); newArgs[8] =
         * Utils.getCmdParameter(args, "thirdFile", true); newArgs[9] =
         * Utils.getCmdParameter(args, "thirdSearchEngine", true); newArgs[10] =
         * Utils.getCmdParameter(args, "thirdcvTerm", true); newArgs[11] =
         * Utils.getCmdParameter(args, "thirdbetterScoresAreLower", true);
         * newArgs[12] = Utils.getCmdParameter(args, "rank", true); newArgs[13]
         * = Utils.getCmdParameter(args, "decoyRatio", true); newArgs[14] =
         * Utils.getCmdParameter(args, "outputFile", true); newArgs[15] =
         * Utils.getCmdParameter(args, "debugFile", true); newArgs[16] =
         * Utils.getCmdParameter(args, "decoyRegex", true);
         */
        String file_1 = args[0];
        String searchEngine_1 = args[1];
        String firstcvTerm = args[2];
        String firstbetterScoresAreLower = args[3];

        String file_2 = args[4];
        String searchEngine_2 = args[5];
        String secondcvTerm = args[6];
        String secondbetterScoresAreLower = args[7];

        String file_3 = args[8];
        String searchEngine_3 = args[9];
        String thirdcvTerm = args[10];
        String thirdbetterScoresAreLower = args[11];

        int rank = Integer.parseInt(args[12]);
        int decoyRatio = Integer.parseInt(args[13]);

        csvFileName = args[14];
        //csvFileName = csvFileName + ".txt";

        String debugFileName = "";
        if (args[15] != null) {
            debugFileName = args[15];
        }

        String[] searchEngine = {searchEngine_1, searchEngine_2, searchEngine_3};
        String[] inputFiles = {file_1, file_2, file_3};
        String[] cvTerms = {firstcvTerm, secondcvTerm, thirdcvTerm};
        String[] betterScoresAreLower = {firstbetterScoresAreLower,
            secondbetterScoresAreLower, thirdbetterScoresAreLower};

        CombineSearchEngines C = new CombineSearchEngines(searchEngine);
        long startTime;
        long stopTime;
        long elapsedTime;
        // Added by FG 8/10/2014
        C.checkInputFiles(inputFiles);
        String decoyRegex = args[16];
        mzidVer = MzIdentMLVersion.getVersion(args[17]);
        for (int i = 0; i < searchEngine.length; i++) {
            startTime = System.currentTimeMillis();

            C.computeFDRForSingleSearchEngine(i, inputFiles[i], searchEngine[i],
                                              C.singleFDRInformation[i],
                                              decoyRatio, decoyRegex, cvTerms[i],
                                              betterScoresAreLower[i]);

            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("computeFDRForSingleSearchEngine time for "
                        + searchEngine[i] + " " + elapsedTime / 1000
                        + " Seconds");
            }

            System.out.println("Done - " + searchEngine[i]);
        }

        if (args[15] != null) {
            startTime = System.currentTimeMillis();
            C.writeToFileForDiagnostics_singleFDRObj(debugFileName);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println(
                        "writeToFileForDiagnostics_singleFDRObj time "
                        + elapsedTime / 1000 + " Seconds");
            }

        }
        startTime = System.currentTimeMillis();

        // The common peptides across search engines processing starts here
        C.combinePeptidesAcrossSearchEngines(rank);
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("combinePeptidesAcrossSearchEngines time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();
        // do the sorting according to afs value
        C.sortWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("sortWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.insertFakeDecoyInWholeCombinedResultContainer(); //RK 15-08-12
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println(
                    "insertFakeDecoyInWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();
        C.simpleFdrForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("simpleFdrForWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();
        C.qValueForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("qValueForWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.estFDRForWholeCombinedResultContainer();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("estFDRForWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.removeFakeDecoyFromWholeCombinedResultContainer(); //RK 15-08-12

        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println(
                    "removeFakeDecoyFromWholeCombinedResultContainer time "
                    + elapsedTime / 1000 + " Seconds");
        }

        startTime = System.currentTimeMillis();

        C.prepareTheCSVFileForMzIdentMLParser(csvFileName);
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("prepareTheCSVFileForMzIdentMLParser time "
                    + elapsedTime / 1000 + " Seconds");
        }

    }

//    public static void main(String[] args) throws Exception {
//        if (args.length == 9) {
//            runTwoSearchEngines(args);
//        } else if (args.length == 11) {
//            runThreeSearchEngines(args);
//        } else {
//            System.out.println("Check the Class arguments");
//        }
//
//    }
    public void writeMzidFile() {

        try {
            String outFile = csvFileName;
            if (!outFile.endsWith(".mzid")) {
                outFile = outFile + ".mzid";
            }
            Writer writer = new FileWriter(outFile);

            if (mzidVer == null) {
                mzidVer = MzIdentMLVersion.Version_1_2;
            }

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller(mzidVer);

            writer.write(marshaller.createXmlHeader() + "\n");

            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");

            if (cvList != null) {
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");

            AnalysisSoftware analysisSoftware = new AnalysisSoftware();
            Date date = new Date();
            SimpleDateFormat dateFormat = new SimpleDateFormat(
                    "yyyy-MM-dd HH-mm-ss");
            analysisSoftware.setName(this.getClass().getSimpleName() + "_"
                    + dateFormat.format(date));
            analysisSoftware.setId(this.getClass().getSimpleName() + "_"
                    + dateFormat.format(date));
            Param param = new Param();
            param.setParam(MzidLibUtils.makeCvParam("MS:1002237", "mzidLib",
                                                    CvConstants.PSI_CV));
            analysisSoftware.setSoftwareName(param);
            analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);
            marshaller.marshal(analysisSoftwareList, writer);
            writer.write("\n");

            if (provider != null) {
                marshaller.marshal(provider, writer);
            }
            writer.write("\n");

            if (auditCollection != null) {
                marshaller.marshal(auditCollection, writer);
            }
            writer.write("\n");

            SequenceCollection sequenceCollection = new SequenceCollection();
            Iterator<Entry<String, DBSequence>> itDbSeq
                    = combinedSearchEngines2Mzid.getCombinedDbSequenceHashMap().
                    entrySet().iterator();
            List<DBSequence> dbSeqList = new ArrayList<>();
            while (itDbSeq.hasNext()) {
                Entry<String, DBSequence> pairs = itDbSeq.next();
                dbSeqList.add((DBSequence) pairs.getValue());
                itDbSeq.remove();
            }

            Iterator<Entry<String, Peptide>> itPeptide
                    = combinedSearchEngines2Mzid.getCombinedPeptideHashMap().
                    entrySet().iterator();
            List<Peptide> peptideList = new ArrayList<>();
            while (itPeptide.hasNext()) {
                Entry<String, Peptide> pairs = itPeptide.next();
                peptideList.add((Peptide) pairs.getValue());
                itPeptide.remove();
            }

            Iterator<Entry<String, PeptideEvidence>> itPeptideEvidence
                    = combinedSearchEngines2Mzid.
                    getCombinedPeptideEvidenceHashMap().entrySet().iterator();
            List<PeptideEvidence> peptideEvidenceList = new ArrayList<>();
            while (itPeptideEvidence.hasNext()) {
                Entry<String, PeptideEvidence> pairs = itPeptideEvidence.next();
                //pepEvidID.equals("GASLIQTVLRIASR_generic|TGME49_286932|organism=Toxoplasma_gondii_ME49_3856_3869");
                peptideEvidenceList.add((PeptideEvidence) pairs.getValue());
                itPeptideEvidence.remove();
            }

            sequenceCollection.getPeptideEvidence().addAll(peptideEvidenceList);
            sequenceCollection.getDBSequence().addAll(dbSeqList);
            sequenceCollection.getPeptide().addAll(peptideList);

            marshaller.marshal(sequenceCollection, writer);

            writer.write("\n");
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId("SII_LIST_1");
            // Added by Fawaz Ghali to handle mzid 1.2 version 20/01/2015

            CvParam cvParam = new CvParam();
            cvParam.setAccession("MS:1002439");
            cvParam.setName("final PSM list UNDER DISCUSSION");
            cvParam.setCv(CvConstants.PSI_CV);
            siList.getCvParam().add(cvParam);

            if (analysisCollection != null) {
                List<SpectrumIdentification> spectrumIdentificationList
                        = analysisCollection.getSpectrumIdentification();
                for (int i = 0; i < spectrumIdentificationList.size(); i++) {
                    SpectrumIdentification spectrumIdentification
                            = spectrumIdentificationList.get(i);
                    spectrumIdentification.setSpectrumIdentificationList(siList);

                }
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
                if (mzidVer.equals(MzIdentMLVersion.Version_1_2)) {
                    handleAnalysisProtocolCollection();
                }
                marshaller.marshal(analysisProtocolCollection, writer);
            }
            writer.write("\n");

            writer.write(marshaller.createDataCollectionStartTag() + "\n");

            writer.write("\n");

            if (inputs != null) {
                marshaller.marshal(inputs, writer);
            }
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");

            Iterator<Entry<String, SpectrumIdentificationResult>> itSpectrumIdentificationResult
                    = combinedSearchEngines2Mzid.
                    getCombinedSpectrumIdentificationResultHashMap().entrySet().
                    iterator();
            List<SpectrumIdentificationResult> spectrumIdentificationResultList
                    = new ArrayList<>();
            while (itSpectrumIdentificationResult.hasNext()) {
                Entry<String, SpectrumIdentificationResult> pairs
                        = itSpectrumIdentificationResult.next();
                SpectrumIdentificationResult spectrumIdentificationResult
                        = (SpectrumIdentificationResult) pairs.getValue();
                String spectrumID
                        = spectrumIdentificationResult.getSpectrumID().
                        split(";")[0];
                // 2015-7-20 by Bo Wen: fixed the bug of spectraData mapping
                String location = spectrumIdentificationResult.getSpectrumID().
                        split(";")[1];
                spectrumIdentificationResult.setSpectraData(spectraDataHashMap.
                        get(location));
                spectrumIdentificationResult.setSpectrumID(spectrumID);
                spectrumIdentificationResultList.add(
                        spectrumIdentificationResult);
                itSpectrumIdentificationResult.remove();
            }

            siList.getSpectrumIdentificationResult().addAll(
                    spectrumIdentificationResultList);

            marshaller.marshal(siList, writer);

            writer.write("\n");
            // Update by Fawaz Ghali 17/02/2015 for mzid 1.2 valdiation
            // writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");
            //writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + e.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    private void checkInputFiles(String[] inputFiles) {
        System.out.
                println("Check if the input files have the same spectra data");
        for (int i = 0; i < inputFiles.length; i++) {
            String string = inputFiles[i];
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller
                    = new MzIdentMLUnmarshaller(new File(string));
            Inputs in = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
            List<SpectraData> spectraDataList = in.getSpectraData();
            if (i == 0) {
                for (int j = 0; j < spectraDataList.size(); j++) {
                    SpectraData spectraData = spectraDataList.get(j);
                    String oldID = spectraData.getId();
                    // 2015-7-20 by Bo Wen: change i to j
                    String newID = "SD_COMBINED_SE_" + j;
                    spectraData.setId(newID);
                    File msFile = new File(spectraData.getLocation());
                    // 2015-7-20 by Bo Wen: change to use the spectra file name to map spectraData, 
                    // this change makes it work when the location is different from different search engines.
                    spectraIDLocation.put(oldID, msFile.getName());
                    spectraDataHashMap.put(msFile.getName(), spectraData);
                }
            } else {
                boolean identical = true;

                if (spectraDataList.size() == spectraDataHashMap.keySet().size()) {
                    for (int j = 0; j < spectraDataList.size(); j++) {
                        SpectraData spectraData = spectraDataList.get(j);
                        File msFile = new File(spectraData.getLocation());
                        // 2015-7-20 by Bo Wen: save the spectraID for other mzid files
                        spectraIDLocation.put(spectraData.getId(), msFile.
                                              getName());
                        if (!spectraDataHashMap.keySet().contains(msFile.
                                getName())) {
                            identical = false;

                            break;
                        }

                    }
                } else {
                    identical = false;
                }

                if (!identical) {
                    throw new RuntimeException(
                            "The input files do not have same spectra data.");
                }

            }

        }
    }

}
