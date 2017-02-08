/**
 * @author Fawaz Ghali, Ritesh Krishna , University of Liverpool, 2011
 */
package uk.ac.liv.mzidlib.fdr;

import java.io.BufferedWriter;
import java.io.File;
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

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.multiplesearch.TreeSortForIndices;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

public class FalseDiscoveryRate {

    private boolean verbose = true;
    private int decoyRatio = 1;
    private String decoy;
    //private int siListIndex = 0;
    boolean usingFileDecoyAttribute = false;
    public String allowedEvalues = "MS:1001330" + ";" + "MS:1001172" + ";" + "MS:1001159" + ";" + "MS:1001328" + ";" + "MS:1002045" + ";" + "MS:1002053";
    private String cvTerm = "";
    private boolean betterScoresAreLower = true;
    /*
     * For each Spectrum Result id we need to store :: - peptides associated-
     * scores - evalue	- decoy = true/false
     */
    private List<String> spectrumResult = new ArrayList<>();
    private List<String> spectrumItem = new ArrayList<>();
    private List<String> peptideNames = new ArrayList<>();
    private List<Double> evalues = new ArrayList<>();
    //private ArrayList<Double> scores = new ArrayList<Double>();
    private List<String> decoyOrNot = new ArrayList<>();
    /*
     * Once the evalues/scores list is sorted, we need to remember the original
     * indices of the entries (pre-sorted order) so that we can map back that
     * which value belongs to which spectrumResult,peptideNames etc.
     */
    private Integer[] sortOrderForEvalues;

    /*
     * the above information in sorted order
     */
    private List<String> sorted_spectrumResult = new ArrayList<>();
    private List<String> sorted_spectrumItem = new ArrayList<>();
    private List<String> sorted_peptideNames = new ArrayList<>();
    private List<Double> sorted_evalues = new ArrayList<>();
    private List<Double> sorted_scores = new ArrayList<>();
    private List<String> sorted_decoyOrNot = new ArrayList<>();

    /*
     * Store the estimated FDR, q-Value, FDR Score here
     */
    private List<Double> estimated_simpleFDR = new ArrayList<>();
    private List<Double> estimated_qvalue = new ArrayList<>();
    private List<Double> estimated_fdrscore = new ArrayList<>();
    // F. Ghali adding TP and FP 31/08/2011
    private List<Double> tp = new ArrayList<>();
    private List<Double> fp = new ArrayList<>();
    /*
     * MzIdentML elements
     */
    // private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    //private List<SpectrumIdentificationList> spectrumIdentificationList = new ArrayList();
    private List<SpectrumIdentificationItem> spectrumIdentificationItemList = new ArrayList<>();
    private List<DBSequence> dBSequenceList = new ArrayList<>();
    private List<Peptide> peptideList = new ArrayList<>();
    private Map<String, DBSequence> dBSequenceHashMap = new HashMap<>();
    private Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
    private Map<String, Peptide> peptideHashMap = new HashMap<>();
    private Map<String, SpectrumIdentificationItem> spectrumIdentificationItemHashMap = new HashMap<>();
//    private HashMap<String, SpectrumIdentificationResult> spectrumIdentificationResultHashMap = new HashMap();

    /*
     * The data read from MzIdentML file after parsing
     */
    public Map<String, List<List<String>>> peptideModificationHash;
    public Map<String, String> peptideIdAndSequenceHash;
    public Map<String, List<List<Object>>> spectrumInformationHash;
    private Map<String, String> dbReferenceHash;
    String searchEngineIdentifierScore;
    String searchEngineIdentifierExpectation; // The value used in the original algorithm
    String inputMzIdentMLFile;
    String searchEngine;
    String mascotDecoyTag;
    //metadata
    private AnalysisSoftwareList analysisSoftwareList;
    private AuditCollection auditCollection;
    private Provider provider;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private CvList cvList;
    private AnalysisCollection analysisCollection;
    private Inputs inputs;
    private String searchDatabase_Ref;
    private Map<String, String> unimodHashmap;
    private String spectraData_ref;
    private MzidLibUtils mzidLibUtils;

    public static void main(String args[]) {
    }

    public Map<String, List<List<String>>> getFromXMLPeptideModificationHash() {
        return peptideModificationHash;
    }

    public Map<String, String> getFromXMLPeptideSequenceHash() {
        return peptideIdAndSequenceHash;
    }

    public Map<String, List<List<Object>>> getFromXMLSpectrumInfoHash() {
        return spectrumInformationHash;
    }

    public FalseDiscoveryRate(String mzid, String searchEngine, String decoyRatio, String decoy) {

        this.decoyRatio = Integer.valueOf(decoyRatio);
        this.decoy = decoy;

        this.searchEngine = searchEngine;
        if (decoy == null || decoy.equals("")) {
            usingFileDecoyAttribute = true;
        }

        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid));
            //mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);

            readMzIdentML();

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

//    public FalseDiscoveryRate(String mzid, String searchEngine, String decoyRatio, String decoy, String cvTerm, boolean betterScore) {
//        this.decoyRatio = Integer.valueOf(decoyRatio);
//        this.decoy = decoy;
//        
//        this.searchEngine = searchEngine;
//        this.cvTerm = cvTerm;
//        
//        if(!cvTerm.equals("") && cvTerm!=null){
//            allowedEvalues = cvTerm;            
//        }
//this.betterScoresAreLower =betterScore;
//        
//        if (decoy == null || decoy.equals("")) {
//            usingFileDecoyAttribute = true;
//        }
//        
//        if (searchEngine.equals("mascot")) {
//            searchEngineIdentifierExpectation = "Mascot:expectation value";
//            searchEngineIdentifierScore = "Mascot:score";
//            mascotDecoyTag = null;
//        } else if (searchEngine.equals("omssa")) {
//            searchEngineIdentifierExpectation = "OMSSA:evalue";
//            searchEngineIdentifierScore = "OMSSA:pvalue";
//            mascotDecoyTag = null;
//        } else if (searchEngine.equals("X!Tandem")) {
//            searchEngineIdentifierExpectation = "X\\!Tandem:expect";
//            searchEngineIdentifierScore = "X\\!Tandem:hyperscore";
//            mascotDecoyTag = null;
//        }else{
//            searchEngineIdentifierScore = cvTerm;
//        }
//        
//        try {
//            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid));
//            //mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
//
//            
//            readMzIdentML();
//            System.out.println(cvTerm);
//            
//            
//        } catch (OutOfMemoryError error) {
//            String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
//             String className = this.getClass().getName();
//             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ error.getMessage()+"."+
//                "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
//             System.out.println (message);
//        } 
//    }

    private void readMzIdentML() {
        long startTime = System.currentTimeMillis();
        try {
            mzidLibUtils = new MzidLibUtils();
            peptideIdAndSequenceHash = new HashMap<>();
            peptideModificationHash = new HashMap<>();
            spectrumInformationHash = new HashMap<>();
            dbReferenceHash = new HashMap<>();

            unimodHashmap = new HashMap<>();
            dBSequenceList.clear();
            dBSequenceHashMap.clear();
            peptideHashMap.clear();
            peptideEvidenceHashMap.clear();
            spectrumIdentificationItemHashMap.clear();
//            spectrumIdentificationResultHashMap.clear();
            getUnimodHashmap().clear();
            //cvList = mzIdentML.getCvList();
            cvList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            auditCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            provider = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            analysisCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);

            //Stores SpectraData id and location attribute (one to one).
            HashMap<String, String> spectraDataHaspMap = new HashMap<>();
            List<SpectraData> spectraDataList = inputs.getSpectraData();
            for (int i = 0; i < spectraDataList.size(); i++) {
                SpectraData spectraData = spectraDataList.get(i);
                spectraDataHaspMap.put(spectraData.getId(), spectraData.getLocation());

            }
            searchDatabase_Ref = inputs.getSearchDatabase().get(0).getId();

            Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {

                dBSequenceList.add(iterDBSequence.next());

            }

            for (DBSequence dBSequence : dBSequenceList) {
                dBSequenceHashMap.put(dBSequence.getId(), dBSequence);
                dbReferenceHash.put(dBSequence.getId(), dBSequence.getAccession());
            }

            getPeptideEvidenceMap().clear();
            Iterator<PeptideEvidence> iterPeptideEvidence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {

                PeptideEvidence pe = iterPeptideEvidence.next();
                peptideEvidenceHashMap.put(pe.getId(), pe);
                //peptideEvidenceHashMap.put(iterPeptideEvidence.next().getId(), iterPeptideEvidence.next());
            }

            peptideList.clear();
            getPeptideHashMap().clear();
            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            List<List<String>> modArray = new ArrayList<>(); // There can be multiple modifications for a peptide
            List<String> tempSingleMod = new ArrayList<>();
            List<Modification> modificationmList = new ArrayList<>();

            List<List<Object>> multiSpectrum = new ArrayList<>();
            List<Object> temp = new ArrayList<>();

            String spectrumItemId = new String(), peptideRef = new String();
            int rank = 0, chargedState = 0;
            double calculatedMass = 0.0, experimentalMass = 0.0;
            double expectationValue = 0.0, scoreValue = 0.0;
            int pepStart = 0, pepEnd = 0;
            String pepDBSequence_Ref = new String();
            //boolean isDecoy = false;
            // NOTE - This also has to be a Arraylist, as there can be many peptides
            List<List<Object>> peptideEvd = new ArrayList<>();

            List<Object> singlePeptideEvd = new ArrayList<>();

            boolean isdecoy = false;
            String proteinAccession = "";

            while (iterPeptide.hasNext()) {

                Peptide peptide = iterPeptide.next();
                peptideList.add(peptide);
                String pepId = peptide.getId();
                String pepSeq = peptide.getPeptideSequence();
                modArray.clear();
                tempSingleMod.clear();
                modificationmList.clear();
                modificationmList = peptide.getModification();

                for (Modification modification : modificationmList) {
                    String location = modification.getLocation().toString();
                    String residue = "";
                    if (modification.getResidues() != null && !modification.getResidues().isEmpty()) {
                        residue = modification.getResidues().get(0);
                    }
                    String monoisotopicMassDelta = modification.getMonoisotopicMassDelta().toString();

                    tempSingleMod.add(0, location);
                    tempSingleMod.add(1, residue);
                    tempSingleMod.add(2, monoisotopicMassDelta);

                    List<CvParam> cvParamList = modification.getCvParam();
                    if (cvParamList != null && cvParamList.get(0).getName() != null) {
                        tempSingleMod.add(3, cvParamList.get(0).getName());
                        getUnimodHashmap().put(cvParamList.get(0).getName(), cvParamList.get(0).getAccession());
                    }
                    modArray.add(new ArrayList<>(tempSingleMod));
                    tempSingleMod.clear();
                }
                // Peptide ID and sequence pair; Peptide ID and Modification arrays
                // pepId = pepSeq;
                peptideIdAndSequenceHash.put(pepId, pepSeq);
                peptideModificationHash.put(pepId, new ArrayList<>(modArray));
            }

            for (Peptide peptide1 : peptideList) {
                peptideHashMap.put(peptide1.getId(), peptide1);
            }

//            spectrumIdentificationList.clear();
//            iterSpectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
//            while (iterSpectrumIdentificationList.hasNext()) {
//                spectrumIdentificationList.add();
//            SpectrumIdentificationList spectrumIdentificationList = iterSpectrumIdentificationList.next();
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {

//            if (spectrumIdentificationList != null ) {
//                sirListTemp = iterSpectrumIdentificationResult.next();
//                for (int i = 0; i < sirListTemp.size(); i++) {
                SpectrumIdentificationResult spectrumIdentificationResult = iterSpectrumIdentificationResult.next();
                String spectraRef1 = spectrumIdentificationResult.getSpectraDataRef();
                String location = spectraDataHaspMap.get(spectraRef1);
                // 2015-07-20 by Bo Wen : take the SpectrumID and the spectra file name as the id for each PSM
                // This is very useful for combining results from multiple spectra file searching.
                File locFile = new File(location);
                String spectrumResultId = spectrumIdentificationResult.getSpectrumID() + ";" + locFile.getName();

//                    spectrumIdentificationResultHashMap.put(spectrumResultId, spectrumIdentificationResult);
                temp.clear();
                multiSpectrum.clear();

                peptideEvd.clear();
                for (SpectrumIdentificationItem spectrumIdentItem : spectrumIdentificationResult.getSpectrumIdentificationItem()) {
                    spectrumIdentificationItemList.add(spectrumIdentItem);
                    spectrumIdentificationItemHashMap.put(spectrumIdentItem.getId(), spectrumIdentItem);
                    getSpectrumIdentificationItemHashMap().put(spectrumIdentItem.getId(), spectrumIdentItem);

                    spectrumItemId = spectrumIdentItem.getId();
                    experimentalMass = spectrumIdentItem.getExperimentalMassToCharge();
                    calculatedMass = spectrumIdentItem.getCalculatedMassToCharge();
                    chargedState = spectrumIdentItem.getChargeState();
                    peptideRef = spectrumIdentItem.getPeptideRef();

                    rank = spectrumIdentItem.getRank();

                    List<CvParam> cvParamListSpectrumIdentificationItem = spectrumIdentItem.getCvParam();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentItem.getPeptideEvidenceRef();

                    singlePeptideEvd.clear();
                    isdecoy = false;
                    //String peptideEv;
                    for (PeptideEvidenceRef peptideEvidenceRef : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence = getPeptideEvidenceMap().get(peptideEvidenceRef.getPeptideEvidenceRef());

                        if (peptideEvidence != null) {
                            pepStart = peptideEvidence.getStart();
                            pepEnd = peptideEvidence.getEnd();
                            pepDBSequence_Ref = peptideEvidence.getDBSequenceRef();

                            isdecoy = peptideEvidence.isIsDecoy();

                            //peptideEv = peptideEvidence.getId();

                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence.isIsDecoy()) {
                                    //decoyOrNot.add("false");
                                    isdecoy = true;
                                    break;
                                }
                            } else {

                                if (decoy != null && !decoy.equals("") && (decoy.length() > 1)) {
                                    if (dBSequenceHashMap.get(peptideEvidence.getDBSequenceRef()).getAccession().contains(decoy)) {
                                        //decoyOrNot.add("false");
                                        isdecoy = true;
                                        break;
                                        //                                        peptideEvidence.setIsDecoy(Boolean.TRUE);

                                    }
                                } else {
                                    System.out.println("Error - no decoy value set, need to use alternative constructor");
                                }
                            }
                        }
                    } // end of for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                    if (isdecoy == true) {
                        decoyOrNot.add("true");
                    } else {
                        decoyOrNot.add("false");

                    }

                    // For filliing up the PeptideEvidence - RK
                    for (PeptideEvidenceRef peptideEvidenceRef : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence = getPeptideEvidenceMap().get(peptideEvidenceRef.getPeptideEvidenceRef());

                        if (peptideEvidence != null) {
                            singlePeptideEvd.clear();
                            pepStart = peptideEvidence.getStart();
                            pepEnd = peptideEvidence.getEnd();
                            pepDBSequence_Ref = peptideEvidence.getDBSequenceRef();
                            // Added by FG 09/5/2014 Missing pre post
                            String pre = peptideEvidence.getPre();
                            String post = peptideEvidence.getPost();

                            if (getdBSequenceHashMap().containsKey(pepDBSequence_Ref)) {
                                DBSequence dbSeq = getdBSequenceHashMap().get(pepDBSequence_Ref);
                                proteinAccession = dbSeq.getAccession();
                            }
                            //Removed by FG
                            //isdecoy = peptideEvidence2.isIsDecoy();

                            //peptideEv = peptideEvidence.getId();

                            singlePeptideEvd.add(0, pepStart);
                            singlePeptideEvd.add(1, pepEnd);
                            singlePeptideEvd.add(2, isdecoy);
                            singlePeptideEvd.add(3, proteinAccession);
                            singlePeptideEvd.add(4, pre);
                            singlePeptideEvd.add(5, post);
                            //Added by FG 13/8/2014 fixing multiple spectra ref
                            singlePeptideEvd.add(6, location);
                            peptideEvd.add(new ArrayList<>(singlePeptideEvd));
                        }
                    }

                    spectrumResult.add(spectrumResultId);

                    // RK
                    //if(peptideHashMap.containsKey(peptideRef))
                    //	peptideRef = peptideHashMap.get(peptideRef).getPeptideSequence();
                    peptideNames.add(peptideRef);
                    spectrumItem.add(spectrumIdentItem.getId());
                    for (CvParam cvParam : cvParamListSpectrumIdentificationItem) {
                        String accession = cvParam.getAccession();
                        String scoreExpectationIdentifier = cvParam.getAccession();
                        if (allowedEvalues.contains(accession)) {
                            //System.out.println("Adding e-value");
                            evalues.add(Double.valueOf(cvParam.getValue()));
                            expectationValue = Double.valueOf(cvParam.getValue());
                        }

                        if (scoreExpectationIdentifier.equals(searchEngineIdentifierScore)) {
                            scoreValue = Double.parseDouble(cvParam.getValue());
                        }
                    }
//                        if(spectrumResultId.equals("index=2945")){
//                            System.out.println("found");
//                            System.in.read();
//                        }
                    temp.add(0, spectrumItemId);
                    temp.add(1, peptideRef);
                    temp.add(2, rank);
                    temp.add(3, chargedState);
                    temp.add(4, calculatedMass);
                    temp.add(5, experimentalMass);
                    temp.add(6, expectationValue);
                    temp.add(7, scoreValue);

                    if (!peptideEvd.isEmpty()) {
                        temp.add(8, new ArrayList<>(peptideEvd));
                        peptideEvd.clear();
                    }

                    if (!temp.isEmpty()) {
                        multiSpectrum.add(new ArrayList<>(temp));
                        temp.clear();
                    }

                }

                spectrumInformationHash.put(spectrumResultId, new ArrayList<>(multiSpectrum));
//                }

            }

//        }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("readMzIdentML time " + elapsedTime / 1000 + " Seconds");
        }
    }

    /**
     * Make simple arrayLists from complicated data structures returned after
     * reading mzIdentML file. These simple arrayLists will be used for further
     * calculations of various scores.
     */
    private void getEvalueSortedPeptideList() {

        // RK - The missing copy-paste part
//        Set<String> spectraID = spectrumInformationHash.keySet();
//        Iterator spectraIdIterator = spectraID.iterator();
//        String specRes;
//        ArrayList<ArrayList<Object>> spectrumRes;
//        
//        String pepName;
//        String checkForDecoy;
//        ArrayList<ArrayList<Object>> DBRefs;
//        String decoyCheckForFirstDBRef = "";
//        ArrayList<Object> singleDBRefs;
//        // FG: 04/02/2013 Is it needed? All arrays are already set
//        while (spectraIdIterator.hasNext()) {
//
//            specRes = (String) spectraIdIterator.next();
//            spectrumRes = (ArrayList<ArrayList<Object>>) spectrumInformationHash.get(specRes);
//
//            for (int i = 0; i < spectrumRes.size(); i++) {
//
//                ArrayList<Object> spectrumItem1 = spectrumRes.get(i);
//
//                spectrumResult.add(specRes);
//                spectrumItem.add(spectrumItem1.get(0).toString());
//                pepName = spectrumItem1.get(1).toString();
//                peptideNames.add(pepName);
//                evalues.add(Double.parseDouble(spectrumItem1.get(6).toString()));
//                //scores.add(Double.parseDouble(spectrumItem.get(7).toString()));
//
//                checkForDecoy = "false";
//
//                if (spectrumItem1.size() > 8) {
//                    DBRefs = (ArrayList<ArrayList<Object>>) spectrumItem1.get(8);
//
//
//
//                    int flag = 0;
//                    for (int p = 0; p < DBRefs.size(); p++) {
//                        singleDBRefs = DBRefs.get(p);
//                        checkForDecoy = singleDBRefs.get(2).toString();
//
//                        if (p == 0) {
//                            decoyCheckForFirstDBRef = checkForDecoy;
//                        }
//
//                        if (!checkForDecoy.equalsIgnoreCase(decoyCheckForFirstDBRef)) {
//                          //  System.out.println("Different DB REferences for same Spectrum Item have peptide " + specRes + " --" + pepName + "labelled as Decoy and Non-decoy !");
//                           // System.out.println("Assigning the peptide to be a decoy in this case. Press enter to continue....");
//                            //System.in.read();
//                            flag = 1;
//                        }
//
//                        if (flag == 1) {
//                            checkForDecoy = "false";
//                            break;
//                        }
//                    }
//                } // end of - if( spectrumItem.size() > 8){
//
//                decoyOrNot.add(checkForDecoy);
//            }
//        } // End of while loop
        // RK
        // Call the sorting routine to find the indices of sorted evalues
        TreeSortForIndices sortClass = new TreeSortForIndices();
        sortOrderForEvalues = sortClass.sortTheValueColumn(evalues.toArray(), betterScoresAreLower);

        // Arrange the values in the order determined by the sorting operation, such that, each index an arraylist can be map to
        // the entries in other arraylists for the "same" index
        for (int index : sortOrderForEvalues) {
//            if(spectrumResult.get(index).equals("index=2945"))
//                System.out.println("index=2945");
            sorted_spectrumResult.add(spectrumResult.get(index));
            sorted_spectrumItem.add(spectrumItem.get(index));
            sorted_peptideNames.add(peptideNames.get(index));
            sorted_evalues.add(evalues.get(index));
            //sorted_scores.add(scores.get(index));
            sorted_decoyOrNot.add(decoyOrNot.get(index));
            //System.out.println(spectrumResult.get(index) + " " + evalues.get(index) + " " +decoyOrNot.get(index));

            //System.out.println(sorted_evalues.get(i));
        }

        // Clear the memory for the items no more needed
        spectrumResult.clear();
        peptideNames.clear();
        spectrumItem.clear();
        evalues.clear();
        //scores.clear();
        decoyOrNot.clear();

    }

    // RK - Write the sorted data into a file
    public void writeTheSortedDataToFile(String fileName) throws Exception {

        // Basic check to see that each peptide has a evalue
        if (sorted_peptideNames.size() != sorted_evalues.size()) {
            throw (new Exception("Number of entries = " + sorted_peptideNames.size()
                    + "in sorted_peptideNames don't match with the number of entries = " + sorted_evalues.size() + "in sorted_evalues"));
        }

//        	System.out.println("Peptide = " + sorted_peptideNames.size() + " Evalue = " + sorted_evalues.size()
//        			+ " DecoyorNot = " + sorted_decoyOrNot.size() + " Score = " + sorted_scores.size()
//        			+ " Simple FDR = " + estimated_simpleFDR.size() + " Est QVal = " + estimated_qvalue.size()
//        			+ " Est FDR = " + estimated_fdrscore.size());
        //System.in.read();
        Writer out = new BufferedWriter(new FileWriter(fileName));

        for (int i = 0; i < sorted_peptideNames.size(); i++) {

            String outStr = sorted_spectrumResult.get(i) + "\t"
                    + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i) + "\t"
                    + sorted_evalues.get(i).toString() + "\t"
                    + estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.get(i) + "\t"
                    + // "\n";
                    estimated_fdrscore.get(i) + "\n";

            out.write(outStr);
        }

        out.write("\n\n The Peptide Map");
        String seq;
        for (String key : peptideIdAndSequenceHash.keySet()) {
            seq = peptideIdAndSequenceHash.get(key);
            out.write(key + "\t" + seq + "\n");
        }

        out.close();
    }

    public void writeMzidFile(String csvFileName) {
        try {
            String outFile = csvFileName;
            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller = new MzIdentMLMarshaller();

            writer.write(marshaller.createXmlHeader() + "\n");

            String mzID = mzIdentMLUnmarshaller.getMzIdentMLId();
            if (mzID != null) {
                writer.write(marshaller.createMzIdentMLStartTag(mzID) + "\n");
            } else {
                writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            }

            if (cvList != null) {
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");

            AnalysisSoftware analysisSoftware = new AnalysisSoftware();
            Date date = new Date();
            SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
            analysisSoftware.setName(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
            analysisSoftware.setId(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
            Param param = new Param();
            Cv psiCV = mzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);
            param.setParam(mzidLibUtils.makeCvParam("MS:1002237", "mzidLib", psiCV));
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
            Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {

                sequenceCollection.getDBSequence().add(iterDBSequence.next());

            }

            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide pe = iterPeptide.next();
                sequenceCollection.getPeptide().add(pe);
            }

            Iterator<PeptideEvidence> iterPeptideEvidence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence pe = iterPeptideEvidence.next();
                sequenceCollection.getPeptideEvidence().add(pe);
            }

            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
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

            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Map<String, String> sii_stringMap = new HashMap<>();

            for (int i = 0; i < sorted_spectrumItem.size(); i++) {
                sii_stringMap.put(sorted_spectrumItem.get(i), estimated_simpleFDR.get(i).toString() + ":"
                        + String.valueOf(estimated_qvalue.get(i)) + ":"
                        + estimated_fdrscore.get(i).toString());

            }

            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {

                SpectrumIdentificationResult sr = iterSpectrumIdentificationResult.next();

                for (SpectrumIdentificationItem sii : sr.getSpectrumIdentificationItem()) {

                    CvParam cvParamestimated_simpleFDR = new CvParam();
                    CvParam cvParamestimated_qvalue = new CvParam();
                    CvParam cvParamfdrscore = new CvParam();
                    String sii_temp = sii_stringMap.get(sii.getId());
                    if (sii_temp != null) {
                        String[] sii_arr = sii_temp.split(":");
                        Cv cv = mzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);

                        cvParamestimated_simpleFDR.setAccession("MS:1002351");
                        cvParamestimated_simpleFDR.setName("PSM-level local FDR");
                        cvParamestimated_simpleFDR.setValue(sii_arr[0]);
                        cvParamestimated_simpleFDR.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_simpleFDR);

                        cvParamestimated_qvalue.setAccession("MS:1002354");
                        cvParamestimated_qvalue.setName("PSM-level q-value");
                        cvParamestimated_qvalue.setValue(sii_arr[1]);
                        cvParamestimated_qvalue.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_qvalue);

                        cvParamfdrscore.setAccession("MS:1002355");
                        cvParamfdrscore.setName("PSM-level FDRScore");
                        cvParamfdrscore.setValue(sii_arr[2]);
                        cvParamfdrscore.setCv(cv);
                        sii.getCvParam().add(cvParamfdrscore);

                    }

                }
                siList.getSpectrumIdentificationResult().add(sr);
            }

            marshaller.marshal(siList, writer);
            writer.write("\n");

            writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");

            writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

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

    // Write the sorted data into a file
    public void writeToMzIdentMLFile(String fileName) {

        writeMzidFile(fileName);

    }

    // Write the sorted data into a file
    public void writeToCsvFile(String fileName) throws Exception {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        String outStrHead = "sorted_spectrumResult.get(i),sorted_peptideNames.get(i) , sorted_decoyOrNot.get(i) ,  sorted_evalues.get(i).toString() , + sorted_scores.get(i).toString() , estimated_simpleFDR.get(i) , estimated_qvalue.get(i) , estimated_fdrscore.get(i) \n";

        for (int i = 0; i < sorted_peptideNames.size(); i++) {

            String outStr = sorted_spectrumResult.get(i) + ","
                    + sorted_peptideNames.get(i) + "," + sorted_decoyOrNot.get(i) + ","
                    + sorted_evalues.get(i).toString() + "," //+ sorted_scores.get(i).toString() + ","
                    + estimated_simpleFDR.get(i) + "," + estimated_qvalue.get(i) + ","
                    + estimated_fdrscore.get(i) + "\n";

            out.write(outStr);
        }
        out.close();
    }

    // Write the sorted data into a file
    public void writeToTsvFile(String fileName) throws Exception {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        String outStrHead = "sorted_spectrumResult.get(i)\tsorted_peptideNames.get(i) \t sorted_decoyOrNot.get(i) \t  sorted_evalues.get(i).toString() \t + sorted_scores.get(i).toString() \t estimated_simpleFDR.get(i) \t estimated_qvalue.get(i) \t estimated_fdrscore.get(i) \n";
        String outStr;
        for (int i = 0; i < sorted_peptideNames.size(); i++) {

            outStr = sorted_spectrumResult.get(i) + "\t"
                    + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i) + "\t"
                    //                    + sorted_evalues.get(i).toString() + "\t" + sorted_scores.get(i).toString() + "\t"
                    + estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.get(i) + "\t"
                    + estimated_fdrscore.get(i) + "\n";

            out.write(outStr);
        }
        out.close();
    }

    /**
     * Compute FDR score and q value etc using the method described in Jones et
     * al. Proteomics, 2009,9, 1220-1229
     */
    public void computeFDRusingJonesMethod() {
        try {
            long startTime = System.currentTimeMillis();

            // Sort the data using evalue
            getEvalueSortedPeptideList();

            long stopTime = System.currentTimeMillis();
            long elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("getEvalueSortedPeptideList time " + elapsedTime / 1000 + " Seconds");
            }

            //insertFakeDecoyInTheEnd(); // RK 15-08-12
            // Force the first sorted e-value to be zero.
            //sorted_evalues.set(0, 0f); // RK 21-10-10
            startTime = System.currentTimeMillis();
            computeSimpleFDR(); // Step 1
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("computeSimpleFDR time " + elapsedTime / 1000 + " Seconds");
            }
            //System.out.println(cvTerm);

            if (fp.size() == 0) {
                System.out.println("No decoys found for search engine Mascot|Omssa|tandem - likely caused by: wrong decoy regex, database doesn't contain decoys or the search reported only identifications for stringent e-values - please allow identifications up to e-value = 10 for omssa and tandem. See omssa and tandem documentation for how to do change this setting");
            } else if (fp.get(fp.size() - 1) == 0) {
                System.out.println("No decoys found for search engine Mascot|Omssa|tandem - likely caused by: wrong decoy regex, database doesn't contain decoys or the search reported only identifications for stringent e-values - please allow identifications up to e-value = 10 for omssa and tandem. See omssa and tandem documentation for how to do change this setting");
            }
            startTime = System.currentTimeMillis();
            computeQValues();   // Step 2
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("computeQValues time " + elapsedTime / 1000 + " Seconds");
            }
            startTime = System.currentTimeMillis();
            computeFDRScore();  // Step 3
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                System.out.println("computeFDRScore time " + elapsedTime / 1000 + " Seconds");
            }

            //removeFakeDecoyFromTheEnd(); // RK 15-08-12
            //        for (int i = 0; i < sorted_peptideNames.size(); i++) {
            //            String outStr = sorted_spectrumResult.get(i) + "\t"
            //                    + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i) + "\t"
            //                    + sorted_evalues.get(i).toString() + "\t" + sorted_scores.get(i).toString() + "\t"
            //                    + estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.get(i) + "\t"
            //                    + estimated_fdrscore.get(i) + "\n";
            //
            //            System.out.println(outStr);
            //        }
            //writeToTsvFile("C:\\Users\\fghali\\Desktop\\mzidlib_test_data_20130821\\_"+startTime+"_.txt");
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Compute simple FDR
     *
     */
    private void computeSimpleFDR() {

        int falsePositiveCount = 0;
        int allTargets = 0;
        double falsePositiveDivRatio;
        double simpleFDR;
        for (int i = 0; i < sorted_peptideNames.size(); i++) {
            if (sorted_decoyOrNot.get(i).equals("true")) {
                falsePositiveCount++;
            } else {
                allTargets++;
            }

            falsePositiveDivRatio = (double) falsePositiveCount / (double) decoyRatio; 	//FP = decoy_count / decoy_ratio

            simpleFDR = falsePositiveDivRatio / (double) allTargets;
            estimated_simpleFDR.add(simpleFDR);					//FDR = FP / (FP + TP)		 TP = all_targets - FP	FP + TP = all_targets - FP  + FP		FP + TP = all_targets		FDR = FP / all_targets
            estimated_qvalue.add(0d);

            // fghali
            fp.add(falsePositiveDivRatio);
            double tpValue = allTargets - falsePositiveDivRatio;
            tp.add(tpValue);

            //System.out.println(tpValue + " " + falsePositiveDivRatio + " " + simpleFDR);
        }
    }

    /**
     * Compute q-value
     *
     */
    private void computeQValues() {

        double immediateMinFdr = estimated_simpleFDR.get(estimated_simpleFDR.size() - 1);
        estimated_qvalue.set(estimated_qvalue.size() - 1, immediateMinFdr);
        double currentFDR;
        for (int i = estimated_simpleFDR.size() - 1; i > 0; i--) {
            currentFDR = estimated_simpleFDR.get(i - 1);

            if (currentFDR < immediateMinFdr) {
                immediateMinFdr = currentFDR;
            }

            estimated_qvalue.set(i - 1, immediateMinFdr);
        }
    }

    /**
     * Compute FDR Score
     *
     */
    private void computeFDRScore() {

        // Initialize the estimated_fdrscore by adding elements containing 0. This needs to be
        // done because of back tracking involved in the algorithm, so simple .add() wouldn't work
        for (int i = 0; i < sorted_peptideNames.size(); i++) {
            estimated_fdrscore.add(0d);
        }

        //float prev_evalue = sorted_evalues.get(0); 			// RK 21-10-10
        //float prev_qvalue = estimated_qvalue.get(0);
        //float prev_prev_evalue =  sorted_evalues.get(0); // previous to previous evalue in case of straight vertical rise in q value without any change in evalue
        double prev_evalue = 0f;	    // RK 21-10-10
        double prev_qvalue = 0f;     // RK 21-10-10
        double prev_prev_evalue = 0f;// RK 21-10-10

        int counterBackwardStep = 0; // RK 21-10-10
        //int counter_backwardStep = 1; // First time, the counter is 1 to account for the origin

        //int i = 1;// RK 21-10-10
        int i = 0; // RK 21-10-10
        for (; i < sorted_peptideNames.size(); i++) {
            double current_evalue = sorted_evalues.get(i);
            double current_qvalue = estimated_qvalue.get(i);

            if (current_qvalue > prev_qvalue) {

                double slope;
                double intercept;

                int id = 0;
                if (current_evalue != prev_evalue) {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue - prev_evalue);
                    id = 1; //RK 22-10-10
                } else {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue - prev_prev_evalue);
                    id = 2; //RK 22-10-10
                }

                //intercept = current_qvalue - slope * current_evalue; //RK 22-10-10
                //RK 22-10-10 - Tells us which co-ordinates to use for calculating intercepts.
                if (id == 1) {
                    intercept = prev_qvalue - slope * prev_evalue;
                } else {
                    intercept = prev_qvalue - slope * prev_prev_evalue;
                }

                if (counterBackwardStep > 0) { // compute the FDR score for flat q-value region
                    for (int k = 0; k <= counterBackwardStep; k++) {
                        int index = i - counterBackwardStep + k;
                        //System.out.println("i = " + i + "Count = " + counter_backwardStep +  " k = " + k + " index = " + index);
                        double fdrScore = slope * sorted_evalues.get(index) + intercept;
                        estimated_fdrscore.set(index, fdrScore);

                        //System.out.println("i = " + i + "Count = " + counter_backwardStep +  " k = " + k + " index::1 = " + index + " slope = " + slope + " intercept = " + intercept + " e-val = " + sorted_evalues.get(index) + " FDR = " + fdrScore);
                    }
                } else { 							// In case an immediate increment in q value is found
                    double fdrScore = slope * current_evalue + intercept;
                    estimated_fdrscore.set(i, fdrScore);

                    //System.out.println("index::2 i = " + i + " FDR = " + fdrScore);
                }

                // Re-initialise the variables
                counterBackwardStep = 0;
                if (current_evalue > prev_evalue) { // the previous e-value will change only if the current e-value is different
                    prev_prev_evalue = prev_evalue;
                    prev_evalue = current_evalue;
                }
                prev_qvalue = current_qvalue;

            } else {
                counterBackwardStep++;
            }
        }

        // In case if we miss to update the very last values in estimated_fdrscore because
        // current_qvalue == prev_qvalue
        if (estimated_fdrscore.get(i - 1) == 0) {
            double lastFdrValue = 0;
            i--;
            while (estimated_fdrscore.get(i) == 0 && i > 0) {
                i--;
            }
            if (i == 0) {
                System.out.println("\n Can't compute FDR. Likely that all the PSM hits are corrects, so nothing to do.");
            } else {
                lastFdrValue = estimated_fdrscore.get(i);
                while (i < estimated_fdrscore.size()) {
                    estimated_fdrscore.set(i, lastFdrValue);

                    i++;
                }
            }
        }

    }
    // Added by FG  for sending the specific CV term to use as an alternative mode, also we need a Boolean for -betterScoresAreLower = true or false

    public FalseDiscoveryRate(File mzid, int decoyRatio, String decoy, String cvTerm, boolean betterScoresAreLower) {

        this.decoyRatio = decoyRatio;
        this.decoy = decoy;
        this.cvTerm = cvTerm;

        if (!cvTerm.equals("") && cvTerm != null) {
            allowedEvalues = cvTerm;
        }

        this.betterScoresAreLower = betterScoresAreLower;
        if (decoy == null || decoy.equals("")) {
            usingFileDecoyAttribute = true;
        }

        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(mzid);
            //mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
            System.out.println("Reading: " + mzid);
            System.out.println("CvTerm: " + allowedEvalues);
            readMzIdentML();
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    // Requiring data of a certain type should be performed at the calling code
    @Deprecated
    public FalseDiscoveryRate(String mzid, String decoyRatio, String decoy, String cvTerm, boolean betterScoresAreLower) {
        this(new File(mzid), Integer.parseInt(decoyRatio), decoy, cvTerm, betterScoresAreLower);
    }

    public List<String> getSorted_spectrumResult() {
        return sorted_spectrumResult;
    }

    public List<String> getSorted_peptideNames() {
        return sorted_peptideNames;
    }

    public List<Double> getSorted_evalues() {
        return sorted_evalues;
    }

    public List<Double> getSorted_scores() {
        return sorted_scores;
    }

    public List<String> getSorted_decoyOrNot() {
        return sorted_decoyOrNot;
    }

    public List<Double> getSorted_simpleFDR() {
        return estimated_simpleFDR;
    }

    public List<Double> getSorted_qValues() {
        return estimated_qvalue;
    }

    public List<Double> getSorted_estimatedFDR() {
        return estimated_fdrscore;
    }

    // F Ghali
    public List<Double> getTP() {
        return tp;
    }

    public List<Double> getFP() {
        return fp;
    }

    /**
     * Clear all the data structures in the class to release all the memory
     */
    public void clearAllData() {
        sorted_spectrumResult.clear();
        sorted_peptideNames.clear();
        sorted_evalues.clear();
        sorted_scores.clear();
        sorted_decoyOrNot.clear();
        estimated_simpleFDR.clear();
        estimated_qvalue.clear();
        estimated_fdrscore.clear();
        // fghali
        fp.clear();
        tp.clear();
    }

    public Map<String, DBSequence> getdBSequenceHashMap() {
        return dBSequenceHashMap;
    }

    public void setdBSequenceHashMap(Map<String, DBSequence> dBSequenceHashMap) {
        this.dBSequenceHashMap = dBSequenceHashMap;
    }

    @Deprecated
    public Map<String, PeptideEvidence> getPeptideEvidenceHashMap() {
        return peptideEvidenceHashMap;
    }

    public Map<String, PeptideEvidence> getPeptideEvidenceMap() {
        return peptideEvidenceHashMap;
    }

    public void setPeptideEvidenceHashMap(Map<String, PeptideEvidence> peptideEvidenceHashMap) {
        this.peptideEvidenceHashMap = peptideEvidenceHashMap;
    }

    public Map<String, Peptide> getPeptideHashMap() {
        return peptideHashMap;
    }

    public void setPeptideHashMap(Map<String, Peptide> peptideHashMap) {
        this.peptideHashMap = peptideHashMap;
    }

    public Map<String, SpectrumIdentificationItem> getSpectrumIdentificationItemHashMap() {
        return spectrumIdentificationItemHashMap;
    }

    public void setSpectrumIdentificationItemHashMap(Map<String, SpectrumIdentificationItem> spectrumIdentificationItemHashMap) {
        this.spectrumIdentificationItemHashMap = spectrumIdentificationItemHashMap;
    }

//    public HashMap<String, SpectrumIdentificationResult> getSpectrumIdentificationResultHashMap() {
//        return spectrumIdentificationResultHashMap;
//    }
//    
//    public void setSpectrumIdentificationResultHashMap(HashMap<String, SpectrumIdentificationResult> spectrumIdentificationResultHashMap) {
//        this.spectrumIdentificationResultHashMap = spectrumIdentificationResultHashMap;
//    }
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        return analysisSoftwareList;
    }

    public AuditCollection getAuditCollection() {
        return auditCollection;
    }

    public Provider getProvider() {
        return provider;
    }

    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        return analysisProtocolCollection;
    }

    public CvList getCvList() {
        return cvList;
    }

    public AnalysisCollection getAnalysisCollection() {
        return analysisCollection;
    }

    public Inputs getInputs() {
        return inputs;
    }

    public String getSearchDatabase_Ref() {
        return searchDatabase_Ref;
    }

    public Map<String, String> getUnimodHashmap() {
        return unimodHashmap;
    }

    public String getSpectraData_ref() {
        return spectraData_ref;
    }

}
