package uk.ac.liv.mzidlib.proteogrouper;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

import org.apache.commons.collections.CollectionUtils;

import uk.ac.liv.mzidlib.util.MzidLibUtils;


/*
 * @author Fawaz Ghali, University of Liverpool, 2012
 */
public class ProteoGrouper {

    private MzIdentMLUnmarshaller unmarshaller;
    private List<String> peptideList = new ArrayList<>();         //ArrayList of all Peptide IDs - 17/01/2013 - ARJ changed this from peptide_refs to peptide sequence since different mod strings break the algorithm
    private List<String> dbSequenceList = new ArrayList<>();
    private List<String> matrix = new ArrayList<>();
    //private List<String> matrix_sii = new ArrayList();            //Removed this matrix since the map from PE to SII is many to many 
    private List<String> matrix_pe = new ArrayList<>();               //This is list of all PeptideEvidence elements
    private Map<String, PeptideEvidence> peptideEvidenceIdHashMap = new HashMap<>();      //HashMap with PE IDs and PE objects
    private Map<String, Peptide> peptideIdHashMap = new HashMap<>();                      //HashMap with Peptide IDs and Peptide objects
    private Map<String, List<SpectrumIdentificationItem>> peptideEvidenceIDToSIIIDsHashMap = new HashMap<>();           //Pointers from PE ID to SII IDs
    private Map<String, SpectrumIdentificationItem> idToSIIHashMap = new HashMap<>();
    private Map<String, List<ProteinDetectionHypothesis>> peptide_pdh_HashMap = new HashMap<>();      //Map from peptide sequence to PDHs containing that peptide
    private Map<ProteinDetectionHypothesis, ProteinAmbiguityGroup> pdh_To_Group_Map = new HashMap<>();      //Map from peptide sequence to the PDHs containing that peptide

    private Map<String, DBSequence> dbSequenceIDToDBSequenceHashMap = new HashMap<>();
    private Map<String, Double> pepSeqToBestScoreMap = new HashMap<>();
    //HashMap<ProteinDetectionHypothesis,ArrayList> pdhToUniquePeps = new HashMap();
    private List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = new ArrayList<>();
    private Map<String, List<String>> protAccessionToAllPepEvidIDs = new HashMap<>();            //ARJ added hashMap so that we only retrieve PepEvids for one protein
    private Map<String, PDHSetMember> setMemberList = new HashMap<>();
    private List<ProteinAmbiguityGroup> proteinAmbiguityGroupList;
    private ProteinDetectionList proteinDetectionList = new ProteinDetectionList();
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("Mascot_MSMS_example.mzid");
    //private String inputFile = "build/classes/resources/orbi_T1_F011188.mzid";
    private String outputFile = "output_of_inference.mzid";
    //private String inputFile = "resources/Toxo_prot_inference_omssa.mzid";
    private String inputFile = "resources/Toxo_prot_inference_omssa_different_mods.mzid";

    private String mzidLibVersion = "1.4-beta"; //TODO - this should really be set at runtime
    //private String inputFile = "build/classes/resources/54merge_mascot.mzid";
    //private String inputFile = "build/classes/resources/F012970_iprg2008.mzid";
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("Toxo_prot_inference_omssa.mzid");
    private String cvAccessionForSIIScore = "MS:1001328";   //<cvParam accession="MS:1001172" name="Mascot:expectation value" cvRef="PSI-MS" value="11.5714944925679" /> used for creating protein score
    private boolean logTransformScore = true;           //e.g. do minus log of evalue to make a protein score
    private boolean includeOnlyBestScorePerPeptide = false;
    private Cv psiCV;
    private Cv unitCV;
    private String repProteinCvAcc = "MS:1001591";          //Used for flagging the group representative
    private String repProteinCvName = "anchor protein";

    private String leadingProteinCvName = "leading protein"; //version 1.2 attributes
    private String leadingProteinCvAcc = "MS:1002401";

    private String repProteinV12CvName = "group representative";
    private String repProteinV12CvAcc = "MS:1002403";

    private String nonLeadingProteinCvName = "non-leading protein"; //version 1.2 attributes
    private String nonLeadingProteinCvAcc = "MS:1002402";

    private String countProteinsCvName = "count of identified protein";
    private String countProteinsCvAcc = "MS:1002404";

    private String cvParamNameForProtScore = "ProteoGrouper:PDH score";
    private String cvParamAccForProtScore = "MS:1002235";
    private String cvParamNameForPAGScore = "ProteoGrouper:PAG score";
    private String cvParamAccForPAGScore = "MS:1002236";
    private String cvParamNameForClusterID = "cluster identifier";
    private String cvParamAccForClusterID = "MS:1002407";
    private boolean verbose = false;
    private boolean requireTwoDistinctSequencesForPassThresholdTrue = false;              //This can be changed as a command line parameter at a later date
    private boolean isV11 = true;
    private Map<String, List<ProteinDetectionHypothesis>> clusterMap = new HashMap<>();
    private final String cvParamForDistinctPeptideSequencs = "MS:1001097";
    private ProteinDetection pd;        //ProteinDetection object which will be added to the file
    private ProteinDetectionProtocol pdProtocol;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private AnalysisSoftwareList analysisSoftwareList;
    private AnalysisCollection anCollection;
    private boolean useProteoAnnotator = false;
    private String performanceFile;
    long startTime, stopTime, elapsedTime;
    StringBuffer bf;
    PrintWriter out1 = null;
    /*
     * For testing only
     */

    public static void main(String[] args) {

        //Boolean requireSIIsToPassThreshold = Boolean.parseBoolean(args[0]);
        //TODO - tidy up command line arguments etc.
        ProteoGrouper proteoGrouper = new ProteoGrouper();

    }

    public ProteoGrouper() {

        boolean requireSIIsToPassThreshold = true;
        this.verbose = false;
        this.logTransformScore = true;
        this.init(requireSIIsToPassThreshold);
        this.algorithm();

    }

    public ProteoGrouper(String inputfile, String outputfile, boolean requireSIIsToPassThreshold, boolean verboseOutput, String cvAccForSIIScore, boolean logTransScore) {

        try {
            this.verbose = verboseOutput;
            this.inputFile = inputfile;
            this.outputFile = outputfile;
            this.cvAccessionForSIIScore = cvAccForSIIScore;
            this.logTransformScore = logTransScore;
            //This has been disabled since it is difficult to implement for PAG-level scoring
            //this.includeOnlyBestScorePerPeptide = includeOnlyBestScorePerPep;
            this.includeOnlyBestScorePerPeptide = true;
            performanceFile = outputfile.substring(0, outputfile.lastIndexOf(".")) + "_ProteoGrouperPerformance.txt";
            new File(performanceFile).delete();
            out1 = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile, true)));
            bf = new StringBuffer();
            startTime = System.currentTimeMillis();
            this.init(requireSIIsToPassThreshold);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;

            bf.append("\n\nInit time ").append(elapsedTime / 1000).append(" Seconds");

            this.algorithm();
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            out1.close();
        }

    }

    /*
     * Optional additional attribute in this constructor for creating mzid version 1.2 files if set to isVersion11 set to false
     * 
     */
    public ProteoGrouper(String inputfile, String outputfile, boolean requireSIIsToPassThreshold, boolean verboseOutput, String cvAccForSIIScore, boolean logTransScore, boolean isVersion11, boolean useProteoAnnotator) {
        try {
            this.verbose = verboseOutput;
            this.inputFile = inputfile;
            this.outputFile = outputfile;
            this.cvAccessionForSIIScore = cvAccForSIIScore;
            this.logTransformScore = logTransScore;

            this.isV11 = isVersion11;
            //This has been disabled since it is difficult to implement for PAG-level scoring
            //this.includeOnlyBestScorePerPeptide = includeOnlyBestScorePerPep;
            this.includeOnlyBestScorePerPeptide = true;
            this.useProteoAnnotator = useProteoAnnotator;
            performanceFile = outputfile.substring(0, outputfile.lastIndexOf(".")) + "_ProteoGrouperPerformance.txt";
            new File(performanceFile).delete();
            out1 = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile, true)));
            bf = new StringBuffer();
            startTime = System.currentTimeMillis();
            init(requireSIIsToPassThreshold);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;

            bf.append("\n\nInit time ").append(elapsedTime / 1000).append(" Seconds");
            algorithm();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    private void init(boolean requireSIIsToPassThreshold) {

        //URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("Mascot_MSMS_example.mzid");
        if (inputFile != null) {

            unmarshaller = new MzIdentMLUnmarshaller(new File(inputFile));
            Iterator<Cv> iterCv = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                } else if (cv.getUri().toLowerCase().contains("unit")) {
                    unitCV = cv;
                }
            }

            Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
                peptideIdHashMap.put(peptide.getId(), peptide);
            }

            Iterator<PeptideEvidence> iterPeptideEvidence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence peptideEvidence = iterPeptideEvidence.next();
                peptideEvidenceIdHashMap.put(peptideEvidence.getId(), peptideEvidence);
            }

            Iterator<DBSequence> iterDBSeq = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSeq.hasNext()) {
                DBSequence dbSeq = iterDBSeq.next();
                dbSequenceIDToDBSequenceHashMap.put(dbSeq.getId(), dbSeq);
            }

            Iterator<SpectrumIdentificationList> iterSpectrumIdentificationList = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
            SpectrumIdentificationList siList = iterSpectrumIdentificationList.next();

            System.out.println("Note: The ProteoGrouper is currently only designed to handle a single list of PSMs, and will take the first list it encounters for building protein groups");

            //Code change by ARJ 21/01/2014 - processes only the first list in the file
            List<SpectrumIdentificationResult> siResList = siList.getSpectrumIdentificationResult();

            for (SpectrumIdentificationResult siResult : siResList) {
                for (SpectrumIdentificationItem spectrumIdentificationItem : siResult.getSpectrumIdentificationItem()) {

                    if (!checkSIIContainsValidCVScore(spectrumIdentificationItem, cvAccessionForSIIScore)) {
                        String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
                        String className = this.getClass().getName();
                        String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because the SpectrumIdentificationItem:" + spectrumIdentificationItem.getId() + "did not contain the required CV term:"
                                + cvAccessionForSIIScore + "\nPlease see the reference guide for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
                        System.out.println(message);
                    }

                    if (spectrumIdentificationItem.isPassThreshold() || !requireSIIsToPassThreshold) {

                        idToSIIHashMap.put(spectrumIdentificationItem.getId(), spectrumIdentificationItem);

                        List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                        for (PeptideEvidenceRef peptideEvidenceRef : peptideEvidenceRefList) {
                            PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

                            Peptide peptide = peptideIdHashMap.get(peptideEvidence.getPeptideRef());

                            String pepSeq = peptide.getPeptideSequence();
                            /*
                             * Changed ARJ 17/01/2013 if
                             * (!peptideList.contains(peptideEvidence.getPeptideRef()))
                             * { peptideList.add(peptideEvidence.getPeptideRef()); }
                             *
                             */

                            if (!peptideList.contains(pepSeq)) {
                                peptideList.add(pepSeq);
                            }

                            List<String> peptideEvidIDList = null;
                            if (!dbSequenceList.contains(peptideEvidence.getDBSequenceRef())) {
                                dbSequenceList.add(peptideEvidence.getDBSequenceRef());

                                peptideEvidIDList = new ArrayList<>();
                                protAccessionToAllPepEvidIDs.put(peptideEvidence.getDBSequenceRef(), peptideEvidIDList);        //Create link from a protein accession to all PEs
                            } else {
                                peptideEvidIDList = protAccessionToAllPepEvidIDs.get(peptideEvidence.getDBSequenceRef());
                            }

                            //New code added since it appeared that there was duplication of PeptideEvidences, if multiple SIIs had same PE
                            if (!peptideEvidIDList.contains(peptideEvidence.getId())) {
                                peptideEvidIDList.add(peptideEvidence.getId());
                            }

                            List<SpectrumIdentificationItem> siiList = null;
                            if (!peptideEvidenceIDToSIIIDsHashMap.containsKey(peptideEvidence.getId())) {
                                siiList = new ArrayList<>();
                                peptideEvidenceIDToSIIIDsHashMap.put(peptideEvidence.getId(), siiList);
                            } else {
                                siiList = peptideEvidenceIDToSIIIDsHashMap.get(peptideEvidence.getId());
                            }

                            //ARJ added this code 17/01/2013 - in case there is erroneous duplication of PEs under an
                            if (!siiList.contains(spectrumIdentificationItem)) {
                                siiList.add(spectrumIdentificationItem);
                            }

                            //ARJ note - the matrix is an ArrayList containing DBsequenceID+Peptideref combinations
                            //matrix_pe is a parallel list of all PeptideEvidences
                            //New note 17/1/2013 - altered matrix to contain DBSequence+PeptideSequence
                            if (verbose) {
                                if (!matrix.contains(peptideEvidence.getDBSequenceRef() + pepSeq)) {
                                    matrix.add(peptideEvidence.getDBSequenceRef() + pepSeq);
                                    //matrix_sii.add(spectrumIdentificationItem.getId());                                 //This adds the SII that corresponds with this ID - Note: this will lose all other SIIs with the same mapping
                                    matrix_pe.add(peptideEvidenceRef.getPeptideEvidenceRef());
                                }
                            }
                        }
                    }
                }

            }

            //Iterator<SpectrumIdentificationItem> iterSpectrumIdentificationItem = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);
            //while (iterSpectrumIdentificationItem.hasNext()) {
            //  SpectrumIdentificationItem spectrumIdentificationItem = iterSpectrumIdentificationItem.next();
            /*
             if (spectrumIdentificationItem.isPassThreshold() || !requireSIIsToPassThreshold) {

             idToSIIHashMap.put(spectrumIdentificationItem.getId(), spectrumIdentificationItem);

             List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
             for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
             PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
             PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

             Peptide peptide = peptideIdHashMap.get(peptideEvidence.getPeptideRef());

             String pepSeq = peptide.getPeptideSequence();

             if (!peptideList.contains(pepSeq)) {
             peptideList.add(pepSeq);
             }

             ArrayList<String> peptideEvidIDList = null;
             if (!dbSequenceList.contains(peptideEvidence.getDBSequenceRef())) {
             dbSequenceList.add(peptideEvidence.getDBSequenceRef());

             peptideEvidIDList = new ArrayList();
             protAccessionToAllPepEvidIDs.put(peptideEvidence.getDBSequenceRef(), peptideEvidIDList);        //Create link from a protein accession to all PEs
             } else {
             peptideEvidIDList = protAccessionToAllPepEvidIDs.get(peptideEvidence.getDBSequenceRef());
             }

             //New code added since it appeared that there was duplication of PeptideEvidences, if multiple SIIs had same PE
             if (!peptideEvidIDList.contains(peptideEvidence.getId())) {
             peptideEvidIDList.add(peptideEvidence.getId());
             }

             List<SpectrumIdentificationItem> siiList = null;
             if (!peptideEvidenceIDToSIIIDsHashMap.containsKey(peptideEvidence.getId())) {
             siiList = new ArrayList();
             peptideEvidenceIDToSIIIDsHashMap.put(peptideEvidence.getId(), siiList);
             } else {
             siiList = peptideEvidenceIDToSIIIDsHashMap.get(peptideEvidence.getId());
             }

             //ARJ added this code 17/01/2013 - in case there is erroneous duplication of PEs under an
             if (!siiList.contains(spectrumIdentificationItem)) {
             siiList.add(spectrumIdentificationItem);
             }

             //ARJ note - the matrix is an ArrayList containing DBsequenceID+Peptideref combinations
             //matrix_pe is a parallel list of all PeptideEvidences
             //New note 17/1/2013 - altered matrix to contain DBSequence+PeptideSequence
             if (verbose) {
             if (!matrix.contains(peptideEvidence.getDBSequenceRef() + pepSeq)) {
             matrix.add(peptideEvidence.getDBSequenceRef() + pepSeq);
             //matrix_sii.add(spectrumIdentificationItem.getId());                                 //This adds the SII that corresponds with this ID - Note: this will lose all other SIIs with the same mapping
             matrix_pe.add(peptideEvidenceRef.getPeptideEvidenceRef());
             }
             }
             }
             }
                
             }
             * 
             */
        } else {
            System.out.println("FILE NOT FOUND:" + inputFile);
        }
    }

    /*
     * Helper method for checking that an SII has the score we will require for ordering
     */
    private boolean checkSIIContainsValidCVScore(SpectrumIdentificationItem siItem, String cvAcc) {
        boolean checked = false;

        for (CvParam cvParam : siItem.getCvParam()) {

            if (cvParam.getAccession().equals(cvAcc)) {
                String value = null;
                value = cvParam.getValue();
                if (value != null) {
                    checked = true;
                }
            }
        }

        return checked;

    }

    private void algorithm() {
        startTime = System.currentTimeMillis();
        buildPDHs();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        bf.append("\n\nbuildPDHs time ").append(elapsedTime / 1000).append(" Seconds");
        startTime = System.currentTimeMillis();
        buildPAGs();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        bf.append("\n\nbuildPAGs time ").append(elapsedTime / 1000).append(" Seconds");
        startTime = System.currentTimeMillis();
        sortPAGs();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        bf.append("\n\nsortPAGs time ").append(elapsedTime / 1000).append(" Seconds");

        if (useProteoAnnotator) {
            applyProteoAnnotator();
        }

        startTime = System.currentTimeMillis();
        buildClusters();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        bf.append("\n\nbuildClusters time ").append(elapsedTime / 1000).append(" Seconds");

        if (verbose) {
            System.out.println("x= Number of peptides:" + peptideList.size());
            System.out.println("y= Number of proteins:" + dbSequenceList.size());

            //Print the headers of the matrix (proteins)
            System.out.print("\t");
            for (String dbSequence : dbSequenceList) {
                System.out.print(dbSequence + "\t");
            }
            System.out.print("\n");

            //Matrix
            for (String peptide : peptideList) {

                System.out.print(peptide + "\t");
                for (String dbSequence : dbSequenceList) {
                    if (matrix.contains(dbSequence + peptide)) {
                        List<SpectrumIdentificationItem> siiList = peptideEvidenceIDToSIIIDsHashMap.get(matrix_pe.get(matrix.indexOf(dbSequence + peptide)));
                        String siiString = "";
                        for (SpectrumIdentificationItem sii : siiList) {
                            //System.out.print(sii + ";");
                            siiString += sii.getId() + ";";
                        }
                        siiString = siiString.substring(0, siiString.length() - 1);        //remove last ;
                        System.out.print(siiString + "\t");
                    } else {
                        System.out.print("-\t");
                    }
                }
                System.out.print("\n");
            }
        }

        // PAGs
        /*
         * for (int i = 0; i < proteinAmbiguityGroupList.size(); i++) {
         * System.out.println("PAG: " + i); ProteinAmbiguityGroup
         * proteinAmbiguityGroup = proteinAmbiguityGroupList.get(i);
         *
         * List<ProteinDetectionHypothesis> proteinDetectionHypothesisList =
         * proteinAmbiguityGroup.getProteinDetectionHypothesis(); for (int j =
         * 0; j < proteinDetectionHypothesisList.size(); j++) {
         * System.out.println("PDH: " +
         * proteinDetectionHypothesisList.get(j).getId());
         * ProteinDetectionHypothesis proteinDetectionHypothesis =
         * proteinDetectionHypothesisList.get(j); List<PeptideHypothesis>
         * peptideHypothesisList =
         * proteinDetectionHypothesis.getPeptideHypothesis(); for (int k = 0; k
         * < peptideHypothesisList.size(); k++) {
         *
         * PeptideHypothesis ph = peptideHypothesisList.get(k); //
         * System.out.println(peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef());
         * //
         * System.out.println(matrix_sii.get(matrix.indexOf(proteinDetectionHypothesis.getDBSequenceRef()
         * +
         * peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef()))
         * + "\t"); } } System.out.println("=========================="); }
         */
        startTime = System.currentTimeMillis();
        writeOutputFile();
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        bf.append("\n\nwriteOutputFile time ").append(elapsedTime / 1000).append(" Seconds");

        bf.append("\n\n");
        bf.append("\n\n");

        out1.println(bf.toString());
        out1.close();

    }

    private void buildPDHs() {
        if (verbose) {
            System.out.println("Building PDHs" + " from " + dbSequenceList.size() + " DBSequence records and " + peptideList.size() + " peptides");
        }

        for (int i = 0; i < dbSequenceList.size(); i++) {
            String dbSequence = dbSequenceList.get(i);
            //ProteinDetectionHypothesis

            if (i % 10 == 0 && verbose) {

                System.out.println("\tProcessed " + i + " proteins into PDHs");
            }

            ProteinDetectionHypothesis proteinDetectionHypothesis = new ProteinDetectionHypothesis();
            proteinDetectionHypothesis.setId("PDH_" + i);
            //proteinDetectionHypothesis.setPassThreshold(true);          

            DBSequence dbSeq = new DBSequence();
            dbSeq.setId(dbSequence);

            proteinDetectionHypothesis.setDBSequence(dbSeq);
            double protScore = 0.0;

            boolean passThreshold = true;

            Map<String, Boolean> tempDistinctPeptideSequence = new HashMap<>();    //Need to count the number of distinct peptide sequence - different mods do not count as distinct
            //Key is peptide sequence, value is always true (not needed)

            List<String> pepEvidIDs = protAccessionToAllPepEvidIDs.get(dbSequence);

            //System.out.println("Protein" + dbSequence + " pepEvid" + pepEvidIDs.size());
            //for (int j = 0; j < peptideList.size(); j++) {     
            //String peptide = peptideList.get(j);
            //if (matrix.contains(dbSequence + peptide)) {
            for (String pepEvidID : pepEvidIDs) {
                //distinctPepCounter++;
                PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(pepEvidID);
                String pepSeq = peptideIdHashMap.get(peptideEvidence.getPeptideRef()).getPeptideSequence();

                if (!tempDistinctPeptideSequence.containsKey(pepSeq)) {
                    tempDistinctPeptideSequence.put(pepSeq, true);
                }

                PeptideHypothesis peptideHypothesis = new PeptideHypothesis();

                //lines removed by ARJ 12th Dec, since I think this was incorrect
                //PeptideEvidence peptideEvidence = new PeptideEvidence();
                //peptideEvidence.setId(pepEvidID);
                peptideHypothesis.setPeptideEvidence(peptideEvidence);
                // SpectrumIdentificationItemRef

                List<SpectrumIdentificationItem> siiList = peptideEvidenceIDToSIIIDsHashMap.get(pepEvidID);

                double maxScore = 0.0;

                //TODO - there is one case not covered by my logic, where a peptide is mapped multiple times to the same protein
                //in this case, the score gets multiplied by the number of mappings
                for (SpectrumIdentificationItem sii : siiList) {
                    SpectrumIdentificationItemRef spectrumIdentificationItemRef = new SpectrumIdentificationItemRef();
                    spectrumIdentificationItemRef.setSpectrumIdentificationItemRef(sii.getId());
                    peptideHypothesis.getSpectrumIdentificationItemRef().add(spectrumIdentificationItemRef);

                    double pepScore = 0.0;

                    for (CvParam cvParam : sii.getCvParam()) {
                        if (cvParam.getAccession().equals(cvAccessionForSIIScore)) {
                            pepScore = Double.parseDouble(cvParam.getValue());
                            if (logTransformScore) {
                                // Updated by FG 11/11/2014 Requested by Andy for Infinity value for ProteoGrouper:PDH score
                                if (pepScore == 0) {
                                    pepScore = 300;
                                }
                                pepScore = -10.0 - (10 * Math.log10(pepScore));     //A score of 0.1 thus scores zero - this is based on the assumption of FDR-type scores rather than e-values
                                if (pepScore < 0) {
                                    pepScore = 0.0;
                                }
                            }

                            if (includeOnlyBestScorePerPeptide) {
                                if (pepScore > maxScore) {
                                    maxScore = pepScore;
                                }
                            } else {
                                protScore += pepScore;
                            }
                        }
                    }
                }
                protScore += maxScore;  //only set if includeBest = true;                   
                proteinDetectionHypothesis.getPeptideHypothesis().add(peptideHypothesis);
                //Add score to a map, needed for protein group scoring
                if (pepSeqToBestScoreMap.get(pepSeq) == null) {
                    pepSeqToBestScoreMap.put(pepSeq, maxScore);
                } else {
                    double currentBestScore = pepSeqToBestScoreMap.get(pepSeq);
                    if (maxScore > currentBestScore) {
                        pepSeqToBestScoreMap.put(pepSeq, maxScore);
                    }
                }

                //String peptide = peptideEvidenceIdHashMap.get(pepEvidID).getPeptideRef();     //not sure if this is ID mapped, may need another hashmap
                //List list = peptide_pdh_HashMap.get(peptide); Altered by ARJ 17/01/2013
                String ilMappedpepSeq = pepSeq.replace('L', 'J').replace('I', 'J');      //Treat I or L differences as the same
                List<ProteinDetectionHypothesis> list = peptide_pdh_HashMap.get(ilMappedpepSeq);
                if (list == null) {
                    list = new ArrayList<ProteinDetectionHypothesis>();
                    //peptide_pdh_HashMap.put(peptide, list);
                    peptide_pdh_HashMap.put(ilMappedpepSeq, list);
                }
                list.add(proteinDetectionHypothesis);
                //}                
            }

            if (tempDistinctPeptideSequence.keySet().size() < 2 && requireTwoDistinctSequencesForPassThresholdTrue) {
                passThreshold = false;
            }

            //The number of distinct peptide sequences is needed by the make PAG algorithm 17/01/2013
            proteinDetectionHypothesis.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001097", "distinct peptide sequences", psiCV, "" + tempDistinctPeptideSequence.keySet().size()));
            proteinDetectionHypothesis.setPassThreshold(passThreshold);
            proteinDetectionHypothesis.getCvParam().add(MzidLibUtils.makeCvParam(cvParamAccForProtScore, cvParamNameForProtScore, psiCV, "" + protScore));
            proteinDetectionHypothesisList.add(proteinDetectionHypothesis);

        }

        if (verbose) {
            System.out.print("...finished building PDHs\n");
        }
    }

    /*
     * New method following this algorithm
     *
     * Start with the proteins with the fewest peptides and always look equal
     * and upwards: (Retrieve all the proteins that could be involved in a set
     * relationship from its peptides) if no unique peptides: Test for being a
     * subset Test for being a sameset Test for being multiply subsumed i.e. Are
     * all peptides contained within proteins that have more overall peptides
     *
     * If none of these cases - then it is a master protein
     *
     */
    private void buildPAGs() {

        //First sort total protein list by number of peptides
        proteinAmbiguityGroupList = proteinDetectionList.getProteinAmbiguityGroup();

        List<String> allMastersAndMembersForTesting = new ArrayList<>();         //Arraylist used to debug which PDHs have not been put into groups
        Map<ProteinDetectionHypothesis, Map<ProteinDetectionHypothesis, Integer>> multiplySubsumedPDHs = new HashMap<>();
        Map<ProteinDetectionHypothesis, ProteinAmbiguityGroup> pdhToPAGMap = new HashMap<>();  //Map from PDH to PAG
        Map<ProteinDetectionHypothesis, Boolean> pdhIsRepProt = new HashMap<>();               //Map for whether a PDH is a rep prot

        int numPAGsPassThreshold = 0;
        int pdhAssignedToGroupCount = 0;

        //First pass of the List to find unique peptides and map IDs to PDHs
        for (ProteinDetectionHypothesis pdh : proteinDetectionHypothesisList) {
            List<String> uniquePeptides = new ArrayList<>();     //list of peptides that can only be matched to this protein

            setMemberList.put(pdh.getId(), new PDHSetMember(pdh));

            //Let's find the unique peptides here
            /*
             * Removed in code update by ARJ 1/03/2013 - unique peptides have no
             * extra special status now - these are found in main loop below
             * for(String pepSeq : getPeptidesFromPDH(pdh)){ String ilModpepSeq
             * = pepSeq.replaceAll("L","J" ).replaceAll("I","J" ); //do look-up
             * replacing I and L with ambiguity code
             * List<ProteinDetectionHypothesis> list =
             * peptide_pdh_HashMap.get(ilModpepSeq); if(list.size()==1){
             * if(verbose){ System.out.println("\tPDH:" + pdh.getDBSequenceRef()
             * + " has unique peptide: " + pepSeq); }
             * uniquePeptides.add(pepSeq); } }
             *
             *
             */
        }

        //Code for sorting pdhList first by number of peptides then alphabetically, smallest then later in the alpabet first
        Collections.sort(proteinDetectionHypothesisList, new Comparator<ProteinDetectionHypothesis>() {

            @Override
            public int compare(ProteinDetectionHypothesis pdh1, ProteinDetectionHypothesis pdh2) {
                int i = comparePDHs(pdh1, pdh2, cvParamForDistinctPeptideSequencs);
                return i;
            }
        });

        if (verbose) {
            System.out.println("About to build PAGs for: " + proteinDetectionHypothesisList.size() + " PDHs");
        }

        //Now reprocess PDHs to assign peptides to best protein candidate
        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis outerPDH = proteinDetectionHypothesisList.get(i);

            List<String> razorPeptides = new ArrayList<>();      //list of peptides that are most likely mapped to this protein i.e. all other mapped proteins have less evidence (fewer peptides), 
            //thus this is a "pseudo-unique" peptides for making groups from these proteins

            List<String> uniquePeptides = new ArrayList<>();
            //ArrayList<String> uniquePeptides = pdhToUniquePeps.get(outerPDH);
            //ArrayList<String> tiedPeptides = new ArrayList();       //List of peptides with no evidence either way, arbitrarily assigned by on alphabetical order
            Map<ProteinDetectionHypothesis, Integer> countOccurencesOfPDHForAllPeptidesInThisProtein = new HashMap<>();  //Needed to work out which protein contributes the most to subsuming this protein

            List<String> outerPeptides = getPeptidesFromPDH(outerPDH);
            int numOuterPeptides = outerPeptides.size();

            if (verbose) {
                System.out.println("Processing " + outerPDH.getDBSequenceRef());
            }
            //Boolean isMaster=false;
            boolean sameSetOrSubset = false;
            //Boolean isMultiplySubsumed = true;

            //if(uniquePeptides.size()==0){
            if (verbose) {
                System.out.println("\tPDH:" + outerPDH.getDBSequenceRef() + " has " + numOuterPeptides + " peptides:" + outerPeptides.toString());
            }

            for (int j = proteinDetectionHypothesisList.size() - 1; j > i; j--) {
                ProteinDetectionHypothesis innerPDH = proteinDetectionHypothesisList.get(j);
                //code = getSameorSubsets(outerPDH,innerPDH);

                List<String> innerPeptides = getPeptidesFromPDH(innerPDH);

                if (outerPeptides.size() == innerPeptides.size()) {
                    if (CollectionUtils.isEqualCollection(outerPeptides, innerPeptides)) {
                        if (verbose) {
                            System.out.println("\tSame set: " + outerPDH.getDBSequenceRef() + " contained by " + innerPDH.getDBSequenceRef());
                        }
                        sameSetOrSubset = true;
                        PDHSetMember candidateMasterMember = setMemberList.get(innerPDH.getId());
                        PDHSetMember samesetMember = setMemberList.get(outerPDH.getId());
                        samesetMember.setIsSameSet(true);
                        addSubMemberRecursively(candidateMasterMember, samesetMember);
                        //candidateMasterMember.getAllContainedMembers().add(samesetMember);          //Now add 

                        j = i;
                    }
                } else if (outerPeptides.size() < innerPeptides.size()) {     //because the PDH are ordered from big to small (num peptides) - subsets will be set before same set
                    if (CollectionUtils.isProperSubCollection(outerPeptides, innerPeptides)) {
                        if (verbose) {
                            System.out.println("\tPDH: " + outerPDH.getDBSequenceRef() + " is a subset of " + innerPDH.getDBSequenceRef());
                        }
                        PDHSetMember candidateMasterMember = setMemberList.get(innerPDH.getId());
                        PDHSetMember subsetMember = setMemberList.get(outerPDH.getId());
                        subsetMember.setIsSubset(true);
                        //candidateMasterMember.getAllContainedMembers().add(subsetMember);          //Now add 
                        addSubMemberRecursively(candidateMasterMember, subsetMember);
                        j = i;
                        sameSetOrSubset = true;
                    }
                }
            }

            //The proteins that are same or subsets are not going to own any peptides. We now assign ownership of peptides to PDHs
            if (!sameSetOrSubset) {

                for (String pepSeq : outerPeptides) {
                    String ilModpepSeq = pepSeq.replace('L', 'J').replace('I', 'J');
                    List<ProteinDetectionHypothesis> listOfProteinsMappedFromSamePeptide = peptide_pdh_HashMap.get(ilModpepSeq);

                    if (listOfProteinsMappedFromSamePeptide.size() == 1) {
                        if (verbose) {
                            System.out.println("\tPDH:" + outerPDH.getDBSequenceRef() + " has unique peptide: " + pepSeq);
                        }
                        uniquePeptides.add(pepSeq);
                    } else {

                        //Check if ALL of these peptides are mapped to a protein with stronger evidence
                        //If any peptide is not mapped to a protein with more peptides, then this can be a group representative
                        //int maxPeptidesMappedForProteinsContainingThisPeptide = numOuterPeptides; //remove in code update by ARJ 1/3/2013
                        ProteinDetectionHypothesis maxPDHForThisPeptide = outerPDH;
                        //TODO - flaw in the code logic, we should only be looking further down the list, not back
                        for (ProteinDetectionHypothesis innerPDH : listOfProteinsMappedFromSamePeptide) {
                            if (innerPDH != outerPDH) {
                                //This code is only used for working out which is the most subsuming protein
                                if (!countOccurencesOfPDHForAllPeptidesInThisProtein.containsKey(innerPDH)) {
                                    countOccurencesOfPDHForAllPeptidesInThisProtein.put(innerPDH, 1);
                                } else {
                                    int count = countOccurencesOfPDHForAllPeptidesInThisProtein.get(innerPDH);
                                    count++;
                                    countOccurencesOfPDHForAllPeptidesInThisProtein.put(innerPDH, count);
                                }

                                /*
                                 * int numPeptidesMappedToInnerProtein =
                                 * getPeptidesFromPDH(innerPDH).size();
                                 * if(numPeptidesMappedToInnerProtein >
                                 * maxPeptidesMappedForProteinsContainingThisPeptide){
                                 * maxPeptidesMappedForProteinsContainingThisPeptide
                                 * = numPeptidesMappedToInnerProtein;
                                 * maxPDHForThisPeptide = innerPDH; } else
                                 * if(numPeptidesMappedToInnerProtein ==
                                 * maxPeptidesMappedForProteinsContainingThisPeptide
                                 * &&
                                 * comparePDHs(innerPDH,maxPDHForThisPeptide,cvParamForDistinctPeptideSequencs)
                                 * > 0){ maxPDHForThisPeptide = innerPDH; //i.e.
                                 * there is a tie for which is the best protein
                                 * for a given peptide, so we take the best
                                 * scoring }
                                 *
                                 */
                                if (comparePDHs(innerPDH, maxPDHForThisPeptide, cvParamForDistinctPeptideSequencs) > 0) {
                                    maxPDHForThisPeptide = innerPDH;      //Inner PDH has stronger evidence that outerPDH, so will claim this peptide
                                }
                            }
                        }

                        if (verbose && maxPDHForThisPeptide != outerPDH) {
                            System.out.println(" This protein:" + maxPDHForThisPeptide.getDBSequenceRef() + " has more evidence for owning pep seq:" + pepSeq + " hence this is not a razor peptide");
                        } else if (verbose) {
                            System.out.println(" no better candidate protein for owning pep seq:" + pepSeq + " so this assigned as a razor for this protein");
                        }

                        if (maxPDHForThisPeptide == outerPDH) {
                            razorPeptides.add(pepSeq);
                            if (verbose) {
                                System.out.println("Razor peptide:" + pepSeq + " for:" + outerPDH.getDBSequenceRef() + " based on having more peptides than any other, better score or alphabetical order");
                            }
                        }
                    }
                }
            }

            //}  //end of if for uniquePeptide.size==0
            if (verbose) {

                System.out.println("Unique:" + uniquePeptides.size() + " razor:" + razorPeptides.size());
            }

            if (uniquePeptides.size() > 0 || razorPeptides.size() > 0) {
                ProteinAmbiguityGroup pag = new ProteinAmbiguityGroup();
                numPAGsPassThreshold++;
                //Needed for v1.2 files
              //  if (!isV11) {
                    pag.getCvParam().add(MzidLibUtils.makeCvParam("MS:1002415", "protein group passes threshold", psiCV, "true"));
                //}

                List<String> allScoringPeps = new ArrayList<>();
                allScoringPeps.addAll(uniquePeptides);
                allScoringPeps.addAll(razorPeptides);
                addPAGScore(allScoringPeps, pag);

                if (verbose) {
                    System.out.println("Creating PAG with PDH " + outerPDH.getDBSequenceRef() + " as representative protein");
                }
                PDHSetMember masterMember = setMemberList.get(outerPDH.getId());
                List<ProteinDetectionHypothesis> pdhList = pag.getProteinDetectionHypothesis();
                List<PDHSetMember> containedMembers = masterMember.getAllContainedMembers();
                //setAllContainedMembersRecursive(masterMember,containedMembers);    //Should not be needed but safer to make a call to this to check if any contained members need pulling up

                pdhList.add(outerPDH);
                pdh_To_Group_Map.put(outerPDH, pag); //Needed for cluster formation

                pdhToPAGMap.put(outerPDH, pag);
                pdhIsRepProt.put(outerPDH, true);

                if (isV11) {
                    outerPDH.getCvParam().add(MzidLibUtils.makeCvParam(repProteinCvAcc, repProteinCvName, psiCV));
                } else {
                    outerPDH.getCvParam().add(MzidLibUtils.makeCvParam(leadingProteinCvAcc, leadingProteinCvName, psiCV));
                    outerPDH.getCvParam().add(MzidLibUtils.makeCvParam(repProteinV12CvAcc, repProteinV12CvName, psiCV));
                }

                if (uniquePeptides.size() > 0) {
                    //outerPDH.getCvParam().add(makeCvParam("MS:99999","uniquePeptides",psiCV,uniquePeptides.toString().replaceAll(",", ";").replaceAll("\\[", "").replaceAll("\\]", "")));
                    outerPDH.getUserParam().add(MzidLibUtils.makeUserParam("unique peptides", uniquePeptides.toString().replace(',', ';').replaceAll("\\[", "").replaceAll("\\]", "")));
                }

                if (razorPeptides.size() > 0) {
                    outerPDH.getUserParam().add(MzidLibUtils.makeUserParam("razor peptides", razorPeptides.toString().replace(',', ';').replaceAll("\\[", "").replaceAll("\\]", "")));

                    //outerPDH.getCvParam().add(makeCvParam("MS:99999","razorPeptides",psiCV,razorPeptides.toString().replaceAll(",", ";").replaceAll("\\[", "").replaceAll("\\]", "")));
                }


                /*
                 * ARJ 28/02/2013
                 *
                 * May want to add another score in here for uniqueness score or
                 * razor score, based only on those peptides
                 *
                 */
                /*
                 * if(tiedPeptides.size()>0){
                 * outerPDH.getCvParam().add(makeCvParam("MS:99999","tiedPeptides",psiCV,tiedPeptides.toString().replaceAll(",",
                 * ";").replaceAll("\\[", "").replaceAll("\\]", ""))); }
                 *
                 */
                allMastersAndMembersForTesting.add(outerPDH.getId());
                pdhAssignedToGroupCount++;

                for (PDHSetMember member : containedMembers) {

                    ProteinDetectionHypothesis pdh = member.getPdh();
                    allMastersAndMembersForTesting.add(pdh.getId());
                    pdhList.add(pdh);
                    pdh_To_Group_Map.put(pdh, pag);
                    pdhAssignedToGroupCount++;

                    if (member.getIsSameSet()) {
                        //ARJ Jan 2014: New recommendations state that reference should be given to superset PDH ID
                        //pdh.getCvParam().add(utils.makeCvParam("MS:1001594", "sequence same-set protein", psiCV, dbSequenceIDToDBSequenceHashMap.get(outerPDH.getDBSequenceRef()).getAccession()));
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001594", "sequence same-set protein", psiCV, outerPDH.getId()));
                        if (!isV11) {//In version 1.2, we can flag that these are equally good - also leading proteins
                            pdh.getCvParam().add(MzidLibUtils.makeCvParam(leadingProteinCvAcc, leadingProteinCvName, psiCV));
                        }

                    } else if (member.getIsSubset()) {
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001596", "sequence sub-set protein", psiCV, outerPDH.getId()));
                        if (!isV11) {
                            pdh.getCvParam().add(MzidLibUtils.makeCvParam(nonLeadingProteinCvAcc, nonLeadingProteinCvName, psiCV));
                        }
                    } else if (member.getIsMultiplySubsumed()) {
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001598", "sequence subsumable protein", psiCV, outerPDH.getId()));
                        if (!isV11) {
                            pdh.getCvParam().add(MzidLibUtils.makeCvParam(nonLeadingProteinCvAcc, nonLeadingProteinCvName, psiCV));
                        }
                    } else {
                        System.out.println("Error - protein is not a master and not sameset, subset or multiply subsumed?");
                    }
                }

                //pag.setId("PAG_" + pagCount);
                proteinAmbiguityGroupList.add(pag);
                //pagCount++;
            } //else if(tiedPeptides.size()<numOuterPeptides && !sameSetOrSubset){
            else if (!sameSetOrSubset) {
                if (verbose) {
                    System.out.println("\tProt: " + outerPDH.getDBSequenceRef() + "is multiply sub-sumed since it has no unique, no razor peptides and number of tied peptides is less than full set");
                }
                multiplySubsumedPDHs.put(outerPDH, countOccurencesOfPDHForAllPeptidesInThisProtein);

                /*
                 * int maxPeptides = 0; ProteinDetectionHypothesis maxPDH=null;
                 * for(ProteinDetectionHypothesis pdh :
                 * countOccurencesOfPDHForAllPeptidesInThisProtein.keySet()){
                 * if(countOccurencesOfPDHForAllPeptidesInThisProtein.get(pdh) >
                 * maxPeptides &&
                 * getPeptidesFromPDH(pdh).size()>=numOuterPeptides){
                 * maxPeptides =
                 * countOccurencesOfPDHForAllPeptidesInThisProtein.get(pdh);
                 * maxPDH = pdh; } }
                 *
                 * if(maxPDH!=null){ if(verbose){ System.out.println("\tPDH:" +
                 * outerPDH.getDBSequenceRef() + " is multiply subsumed by" +
                 * maxPDH.getDBSequenceRef() + " since it shares " + maxPeptides
                 * + " peptides " + " this protein has: " + numOuterPeptides); }
                 * PDHSetMember candidateMasterMember =
                 * setMemberList.get(maxPDH.getId()); PDHSetMember
                 * multiplySubsumedMember = setMemberList.get(outerPDH.getId());
                 * multiplySubsumedMember.setIsMultiplySubsumed(true);
                 * //candidateMasterMember.getAllContainedMembers().add(multiplySubsumedMember);
                 * //Now add
                 * addSubMemberRecursively(candidateMasterMember,multiplySubsumedMember);
                 * } else{ System.out.println("No subsuming protein found for: "
                 * + outerPDH.getDBSequenceRef()); }
                 *
                 */
            }

            if (verbose) {
                System.out.println("\tFinal assignment: " + outerPDH.getDBSequenceRef() + " uniquePepides: " + uniquePeptides.toString() + " sameSet or subset: " + sameSetOrSubset + " razors:  " + razorPeptides.toString());
            }
        }

        /*
         * Now handle all multiply sub-sumed and their children Current
         * algorithm implementation is rather long-winded and certainly not
         * perfect, needs work
         *
         */
        for (ProteinDetectionHypothesis msPDH : multiplySubsumedPDHs.keySet()) {

            if (verbose) {

                System.out.println("\tProcessing mult sub prot: " + msPDH.getDBSequenceRef());
            }

            Map<ProteinDetectionHypothesis, Integer> countOccurencesOfPDHForAllPeptidesInThisProtein = multiplySubsumedPDHs.get(msPDH);

            int maxPeptides = 0;
            int numPeptides = getPeptidesFromPDH(msPDH).size();

            ProteinDetectionHypothesis repProtein = null;

            String subsumingProteins = "";

            for (ProteinDetectionHypothesis pdh : countOccurencesOfPDHForAllPeptidesInThisProtein.keySet()) {

                //System.out.println("\t\tpdh:" + pdh.getDBSequenceRef() +" is rep?" + pdhIsRepProt.get(pdh));
                if (pdhIsRepProt.get(pdh) != null) {

                    if (getPeptidesFromPDH(pdh).size() >= numPeptides) {
                        //subsumingProteins += dbSequenceIDToDBSequenceHashMap.get(pdh.getDBSequenceRef()).getAccession() + ";";
                        subsumingProteins += pdh.getId() + ";";
                        if (countOccurencesOfPDHForAllPeptidesInThisProtein.get(pdh) > maxPeptides) {
                            maxPeptides = countOccurencesOfPDHForAllPeptidesInThisProtein.get(pdh);
                            repProtein = pdh;
                        }
                    }
                }
            }

            if (subsumingProteins.length() > 0) {
                subsumingProteins = subsumingProteins.substring(0, subsumingProteins.length() - 1);  //remove final ;
            }

            if (repProtein != null) {
                ProteinAmbiguityGroup pag = pdhToPAGMap.get(repProtein);
                List<ProteinDetectionHypothesis> pdhList = pag.getProteinDetectionHypothesis();
                pdhList.add(msPDH);
                pdh_To_Group_Map.put(msPDH, pag);
                msPDH.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001598", "sequence subsumable protein", psiCV, subsumingProteins));

                //In version 1.2, need to flag it is non-leading
                if (!isV11) {
                    msPDH.getCvParam().add(MzidLibUtils.makeCvParam(nonLeadingProteinCvAcc, nonLeadingProteinCvName, psiCV));
                }
                allMastersAndMembersForTesting.add(msPDH.getId());
                if (verbose) {
                    System.out.println("Adding " + msPDH.getDBSequenceRef() + " as sub-sumed by:" + repProtein.getDBSequenceRef());
                }
                pdhAssignedToGroupCount++;

                PDHSetMember msMember = setMemberList.get(msPDH.getId());
                List<PDHSetMember> containedMembers = msMember.getAllContainedMembers();
                for (PDHSetMember member : containedMembers) {

                    ProteinDetectionHypothesis pdh = member.getPdh();
                    allMastersAndMembersForTesting.add(pdh.getId());
                    pdhList.add(pdh);
                    pdhAssignedToGroupCount++;
                    pdh_To_Group_Map.put(pdh, pag); //ARJ 02072014 - seems this case was missed

                    if (!isV11) { //rest of proteins contained by this one, must be non-leading in version 1.2 specs
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam(nonLeadingProteinCvAcc, nonLeadingProteinCvName, psiCV));
                    }

                    if (member.getIsSameSet()) {
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001594", "sequence same-set protein", psiCV, msPDH.getId()));
                    } else if (member.getIsSubset()) {
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001596", "sequence sub-set protein", psiCV, msPDH.getId()));
                    } else if (member.getIsMultiplySubsumed()) {
                        pdh.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001598", "sequence subsumable protein", psiCV, subsumingProteins));
                    } else {
                        System.out.println("Error - protein is not a master and not sameset, subset or multiply subsumed?");
                    }
                }

            } else {
                if (verbose) {
                    System.out.println("Should add " + msPDH.getDBSequenceRef() + " as sub-sumed but no rep protein found");
                }
            }

            //countOccurencesOfPDHForAllPeptidesInThisProtein
        }

//        System.out.println("Number of PDHs assigned to groups:" + pdhAssignedToGroupCount);
//        System.out.println("Total PDH:" + proteinDetectionHypothesisList.size());
        //if (!isV11) {
            proteinDetectionList.getCvParam().add(MzidLibUtils.makeCvParam(countProteinsCvAcc, countProteinsCvName, psiCV, "" + numPAGsPassThreshold));
        //}
        //System.out.println("Total PAGs:" + pagCount);

//        for (ProteinDetectionHypothesis pdh : proteinDetectionHypothesisList) {
//            if (!allMastersAndMembersForTesting.contains(pdh.getId())) {
//                System.out.println("PDH: " + pdh.getDBSequenceRef() + " not assigned to a group");
//            }
//        }
    }

    //Takes array of peptides (uniques and razors combined), sum scores and adds it to the PAG
    private double addPAGScore(List<String> peps, ProteinAmbiguityGroup pag) {
        double pagScore = 0.0;

        for (String pep : peps) {
            pagScore += pepSeqToBestScoreMap.get(pep);
        }
        pag.getCvParam().add(MzidLibUtils.makeCvParam(cvParamAccForPAGScore, cvParamNameForPAGScore, psiCV, "" + pagScore));

        return pagScore;
    }

    /*
     * This assumes scores are ordered high to low (high being better)
     */
    private void sortPAGs() {
        Collections.sort(proteinAmbiguityGroupList, new Comparator<ProteinAmbiguityGroup>() {

            @Override
            public int compare(ProteinAmbiguityGroup pag1, ProteinAmbiguityGroup pag2) {
                int i = 0;

                String cvAcc = "";
                if (isV11) {
                    cvAcc = repProteinCvAcc;
                } else {
                    cvAcc = repProteinV12CvAcc;
                }
                ProteinDetectionHypothesis pdh1 = getRepresentativePDH(pag1, cvAcc);
                ProteinDetectionHypothesis pdh2 = getRepresentativePDH(pag2, cvAcc);

                if (getScoreFromPAG(pag1, cvParamAccForPAGScore) > getScoreFromPAG(pag2, cvParamAccForPAGScore)) {
                    i = -1;

                } else if (getScoreFromPAG(pag1, cvParamAccForProtScore) < getScoreFromPAG(pag2, cvParamAccForPAGScore)) {
                    i = +1;
                } else {
                    i = pdh2.getDBSequenceRef().compareTo(pdh1.getDBSequenceRef());
                }

                System.out.println(pdh1.getId() + "," + pdh1.getDBSequenceRef() + "," + getScoreFromPAG(pag1, cvParamAccForPAGScore));
                System.out.println(pdh2.getId() + "," + pdh2.getDBSequenceRef() + "," + getScoreFromPAG(pag2, cvParamAccForPAGScore));

                /*
                 * ProteinDetectionHypothesis pdh1 = getRepresentativePDH(pag1,
                 * repProteinCvAcc); ProteinDetectionHypothesis pdh2 =
                 * getRepresentativePDH(pag2, repProteinCvAcc);
                 *
                 * if (getScoreFromPDH(pdh1, cvParamAccForProtScore) <
                 * getScoreFromPDH(pdh2, cvParamAccForProtScore)) { i = -1; }
                 * else if (getScoreFromPDH(pdh1, cvParamAccForProtScore) >
                 * getScoreFromPDH(pdh2, cvParamAccForProtScore)) { i = +1; }
                 * else{ i =
                 * pdh2.getDBSequenceRef().compareTo(pdh1.getDBSequenceRef()); }
                 *
                 */
                return i;
            }
        });

        int pagCount = 0;
        for (ProteinAmbiguityGroup pag : proteinAmbiguityGroupList) {
            pag.setId("PAG_" + pagCount);
            pagCount++;
        }
        System.out.println("Total PAGs:" + pagCount);
    }

    public int comparePDHs(ProteinDetectionHypothesis pdh1, ProteinDetectionHypothesis pdh2, String cvParamForDistinctPeptideSequencs) {
        int i = 0;

        //Insert code for checking if only one of the two has a unique peptide - then put pdh1 or 2 first
        //ArrayList<String> pdh1UniquePeptides = pdhToUniquePeps.get(pdh1);
        //ArrayList<String> pdh2UniquePeptides = pdhToUniquePeps.get(pdh2);
        if (getCountFromPDH(pdh1, cvParamForDistinctPeptideSequencs) < getCountFromPDH(pdh2, cvParamForDistinctPeptideSequencs)) {
            i--;
        } else if (getCountFromPDH(pdh1, cvParamForDistinctPeptideSequencs) > getCountFromPDH(pdh2, cvParamForDistinctPeptideSequencs)) {
            i++;
        } else {   //tied number of distinct peptide

            if (getScoreFromPDH(pdh1, cvParamAccForProtScore) < getScoreFromPDH(pdh2, cvParamAccForProtScore)) {
                i--;
            } else if (getScoreFromPDH(pdh1, cvParamAccForProtScore) > getScoreFromPDH(pdh2, cvParamAccForProtScore)) {
                i++;
            } else {
                i = pdh2.getDBSequenceRef().compareTo(pdh1.getDBSequenceRef());
            }
        }

        /*
         * Changed by ARJ 1/03/2013 if(pdh1UniquePeptides.size()> 0 &&
         * pdh2UniquePeptides.size()==0){ i=+1; } else
         * if(pdh2UniquePeptides.size()> 0 && pdh1UniquePeptides.size()==0){
         * i=-1; } else
         * if(getCountFromPDH(pdh1,cvParamForDistinctPeptideSequencs) <
         * getCountFromPDH(pdh2,cvParamForDistinctPeptideSequencs)){ i=-1; }
         * else if(getCountFromPDH(pdh1,cvParamForDistinctPeptideSequencs) >
         * getCountFromPDH(pdh2,cvParamForDistinctPeptideSequencs)){ i=+1; }
         * else{ //tied number of distinct peptide
         *
         * if (getScoreFromPDH(pdh1, cvParamAccForProtScore) <
         * getScoreFromPDH(pdh2, cvParamAccForProtScore)) { i = -1; } else if
         * (getScoreFromPDH(pdh1, cvParamAccForProtScore) >
         * getScoreFromPDH(pdh2, cvParamAccForProtScore)) { i = +1; } else{ i =
         * pdh2.getDBSequenceRef().compareTo(pdh1.getDBSequenceRef()); } }
         *
         *
         */
        //write code to compare next by alphabetical order ... complete putting PDHs into correct PAGs etc...PAGs
        return i;
    }

    private ProteinDetectionHypothesis getRepresentativePDH(ProteinAmbiguityGroup pag, String cvAccForRep) {

        ProteinDetectionHypothesis repPDH = null;
        for (ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()) {
            for (CvParam cvParam : pdh.getCvParam()) {
                if (cvParam.getAccession().equals(cvAccForRep)) {
                    repPDH = pdh;
                    break;
                }
            }
        }
        return repPDH;
    }

    /*
     * Helper method to retrieve a double score value from a given cvParam
     */
    private double getScoreFromPDH(ProteinDetectionHypothesis pdh, String cvAccForScore) {

        double score = 0.0;

        for (CvParam cvParam : pdh.getCvParam()) {
            if (cvParam.getAccession().equals(cvAccForScore)) {
                score = Double.parseDouble(cvParam.getValue());
                break;
            }
        }
        return score;
    }

    /*
     * Helper method to retrieve a double score value from a given cvParam
     */
    private double getScoreFromPAG(ProteinAmbiguityGroup pag, String cvAccForScore) {

        double score = 0.0;

        for (CvParam cvParam : pag.getCvParam()) {
            if (cvParam.getAccession().equals(cvAccForScore)) {
                score = Double.parseDouble(cvParam.getValue());
                break;
            }
        }
        return score;
    }

    /*
     * Helper method to retrieve a double score value from a given cvParam
     */
    private int getCountFromPDH(ProteinDetectionHypothesis pdh, String cvAccForScore) {

        int count = -1;

        for (CvParam cvParam : pdh.getCvParam()) {
            if (cvParam.getAccession().equals(cvAccForScore)) {
                count = Integer.parseInt(cvParam.getValue());
                break;
            }
        }
        return count;
    }

    /*
     * Method to assign cluster identifiers to all PDHs
     *
     */
    private void buildClusters() {

        int clusterID = 0;
        boolean newClusterAssigned = false;

        Map<ProteinAmbiguityGroup, Boolean> groupPlacedInCluster = new HashMap<>();

        //
        for (int i = 0; i < peptideList.size(); i++) {
            String pepSeq = peptideList.get(i);
            String ilModpepSeq = pepSeq.replace('L', 'J').replace('I', 'J');
            List<ProteinDetectionHypothesis> pdhList = peptide_pdh_HashMap.get(ilModpepSeq);
            //System.out.println("Peptide: " + pEvid);

            Map<String, List<PDHSetMember>> clusterIDMap = new HashMap<>();
            for (ProteinDetectionHypothesis pdh : pdhList) {
                //System.out.println("Protein: " + pdh.getDBSequenceRef());

                PDHSetMember pdhMember = setMemberList.get(pdh.getId());
                String memberClusterID = pdhMember.getClusterID();
                if (memberClusterID == null) {
                    memberClusterID = "" + clusterID;
                    pdhMember.setClusterID(memberClusterID);
                    newClusterAssigned = true;
                }
                List<PDHSetMember> clusterMembers = clusterIDMap.get(memberClusterID);
                if (clusterMembers == null) {
                    clusterMembers = new ArrayList<>();
                }
                clusterMembers.add(pdhMember);
                clusterIDMap.put(memberClusterID, clusterMembers);
            }

            if (clusterIDMap.size() > 1) {
                //System.out.println("Multiple cluster to combine:" + clusterIDMap.keySet());

                int lowestID = 999999;
                //first loop to establish which is lowest ID
                for (String id : clusterIDMap.keySet()) {
                    int intID = Integer.parseInt(id);
                    if (intID < lowestID) {
                        lowestID = intID;
                    }
                }
                //second loop to set all PDHs to have the lowest cluster ID
                for (String id : clusterIDMap.keySet()) {
                    List<PDHSetMember> clusterMembers = clusterIDMap.get(id);
                    for (PDHSetMember member : clusterMembers) {
                        member.setClusterID("" + lowestID);
                    }
                }

            }
            if (newClusterAssigned) {
                clusterID++;
            }
            newClusterAssigned = false;
        }

        //member.getPdh().getCvParam().add(makeCvParam(cvParamAccForClusterID,cvParamNameForClusterID,psiCV,""+lowestID));
        for (String memberID : setMemberList.keySet()) {
            PDHSetMember member = setMemberList.get(memberID);
            //member.getPdh().getCvParam().add(makeCvParam(cvParamAccForClusterID, cvParamNameForClusterID, psiCV, member.getClusterID()));
            if (!isV11) {
                //member.getPdh().getCvParam().add(utils.makeCvParam(cvParamAccForClusterID,cvParamNameForClusterID, psiCV, member.getClusterID()));
                //Now add the CV term to the parent group if it hasn't already been set
                ProteinAmbiguityGroup parentPAG = pdh_To_Group_Map.get(member.getPdh());

                if (groupPlacedInCluster.get(parentPAG) == null) {

                    parentPAG.getCvParam().add(MzidLibUtils.makeCvParam(cvParamAccForClusterID, cvParamNameForClusterID, psiCV, member.getClusterID()));
                    groupPlacedInCluster.put(parentPAG, true);
                }
            } else {
                member.getPdh().getUserParam().add(MzidLibUtils.makeUserParam(cvParamNameForClusterID, member.getClusterID()));
            }
        }
        //
    }

    /*
     * Recursive call to pull up any contained members underneath one top level
     * member
     */
    private void addSubMemberRecursively(PDHSetMember masterMember, PDHSetMember subMember) {

        if (verbose) {
            System.out.println("Recursive call on master: " + masterMember.getPdh().getDBSequenceRef());
            System.out.println("\tsubMember: " + subMember.getPdh().getDBSequenceRef());
        }

        List<PDHSetMember> masterContainedMembers = masterMember.getAllContainedMembers();

        if (!masterContainedMembers.contains(subMember)) {
            masterContainedMembers.add(subMember);
            if (verbose) {
                System.out.println("\tAdd subMember: " + subMember.getPdh().getDBSequenceRef());
            }

            for (PDHSetMember subSubMember : subMember.getAllContainedMembers()) {
                addSubMemberRecursively(masterMember, subSubMember);
            }
        }

        //subMember.getAllContainedMembers().clear();
    }

    /*
     * Method adapted with algorithm changes 17/01/2013 to return an array of
     * peptide sequences, previously returning array of Peptide refs
     *
     */
    private List<String> getPeptidesFromPDH(ProteinDetectionHypothesis pdh) {

        Map<String, Boolean> distinctPeptideSequence = new HashMap<>();    //Need to count the number of distinct peptide sequence - different mods do not count as distinct
        for (PeptideHypothesis ph : pdh.getPeptideHypothesis()) {

            String pepSeq = peptideIdHashMap.get(peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef()).getPeptideSequence();        //Get peptide sequence, via PE and Pep ID hash maps

            if (!distinctPeptideSequence.containsKey(pepSeq)) {
                distinctPeptideSequence.put(pepSeq, true);
            }

        }

        //Convert HashMap keys to an array
        List<String> peptides = new ArrayList<>();
        for (String pepSeq : distinctPeptideSequence.keySet()) {
            peptides.add(pepSeq);
        }

        //peptides.add(peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef());
        return peptides;
    }

    private void handlePDProtocol() {

        analysisProtocolCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
        analysisSoftwareList = unmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
        anCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);

        pdProtocol = new ProteinDetectionProtocol();
        pdProtocol.setId("PD_Protocol1");

        ParamList threshold = new ParamList();
        threshold.getCvParam().add(MzidLibUtils.makeCvParam("MS:1001494", "no threshold", psiCV));
        pdProtocol.setThreshold(threshold);

        analysisProtocolCollection.setProteinDetectionProtocol(pdProtocol);

        List<AnalysisSoftware> analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();

        AnalysisSoftware analysisSoftware = new AnalysisSoftware();
        analysisSoftware.setName("ProteoGrouper");
        Param tempParam = new Param();
        //tempParam.setParamGroup(makeCvParam("MS:1001475","OMSSA",psiCV));
        tempParam.setParam(MzidLibUtils.makeCvParam("MS:1002241", "ProteoGrouper", psiCV));
        analysisSoftware.setSoftwareName(tempParam);
        analysisSoftware.setId("ProteoGrouper");
        analysisSoftware.setVersion(mzidLibVersion);
        pdProtocol.setAnalysisSoftware(analysisSoftware);

        analysisSoftwares.add(analysisSoftware);

        pd = new ProteinDetection();
        pd.setId("PD1");
        proteinDetectionList.setId("PDL_1");
        pd.setProteinDetectionList(proteinDetectionList);
        pd.setProteinDetectionProtocol(pdProtocol);
        List<InputSpectrumIdentifications> specIdentList = pd.getInputSpectrumIdentifications();
        InputSpectrumIdentifications inputSI = new InputSpectrumIdentifications();
        Iterator<SpectrumIdentificationList> iterSiList = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
        inputSI.setSpectrumIdentificationList(iterSiList.next());
        specIdentList.add(inputSI);
        anCollection.setProteinDetection(pd);

    }

    // Write values
    // Updated by FG 31/01/2013 Solving out of memory error
    private void writeOutputFile() {
        Writer writer = null;
        try {
            writer = new FileWriter(outputFile);
            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();
            writer.write(marshaller.createXmlHeader() + "\n");
            String mzID = unmarshaller.getMzIdentMLId();

            //Creating version 1.2 output           
            if (isV11) {
                if (mzID != null) {
                    writer.write(marshaller.createMzIdentMLStartTag(mzID) + "\n");

                } else {
                    writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
                }
            } else {
                String mzidHeader = "";
                if (mzID != null) {
                    mzidHeader = marshaller.createMzIdentMLStartTag(mzID) + "\n";

                } else {
                    mzidHeader = marshaller.createMzIdentMLStartTag("12345") + "\n";
                }
                mzidHeader = mzidHeader.replace("version=\"1.1.0\"", "version=\"1.2.0\"");    //Replace version number

                writer.write(mzidHeader);
            }

            CvList cvList = unmarshaller.unmarshal(MzIdentMLElement.CvList);
            if (cvList != null) {
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");
            handlePDProtocol(); //Method for setting analysisSoftwareList, analysisCollection and analysisProtocolCollection

            /* Need bespoke code for this, since ProteinDetectionProtocol has to reference this
             if (analysisSoftwareList != null) {
             AnalysisSoftware analysisSoftware = new AnalysisSoftware();
             analysisSoftware.setName(this.getClass().getSimpleName()); 
             analysisSoftware.setId(mzidLibVersion);
             analysisSoftware.setVersion(mzID);
             analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);

             marshaller.marshal(analysisSoftwareList, writer);
             }
             * 
             */
            marshaller.marshal(analysisSoftwareList, writer);
            writer.write("\n");
            Provider provider = unmarshaller.unmarshal(MzIdentMLElement.Provider);
            if (provider != null) {
                marshaller.marshal(provider, writer);
            }
            writer.write("\n");
            AuditCollection auditCollection = unmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            if (auditCollection != null) {
                marshaller.marshal(auditCollection, writer);
            }
            writer.write("\n");
            SequenceCollection sequenceCollection = unmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");
            //AnalysisCollection analysisCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            if (anCollection != null) {
                marshaller.marshal(anCollection, writer);
            }
            writer.write("\n");
            //AnalysisProtocolCollection analysisProtocolCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            if (analysisProtocolCollection != null) {
                marshaller.marshal(analysisProtocolCollection, writer);
            }
            writer.write("\n");
            writer.write(marshaller.createDataCollectionStartTag() + "\n");
            writer.write("\n");
            Inputs inputs = unmarshaller.unmarshal(MzIdentMLElement.Inputs);
            if (inputs != null) {
                marshaller.marshal(inputs, writer);
            }
            writer.write("\n");
            writer.write(marshaller.createAnalysisDataStartTag() + "\n");
            String spectrumIdentificationListRef = "";
            if (anCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = anCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();
            siList.setId(spectrumIdentificationListRef);
            Iterator<FragmentationTable> iterFragmentationTable = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {

                SpectrumIdentificationResult sr = iterSpectrumIdentificationResult.next();

                siList.getSpectrumIdentificationResult().add(sr);
            }
            //            if (siListList != null) {
            //
            //            }
            marshaller.marshal(siList, writer);
            writer.write("\n");
            //writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");
            if (proteinDetectionList != null) {
                marshaller.marshal(proteinDetectionList, writer);
            }
            //writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");
            writer.write(marshaller.createMzIdentMLClosingTag());
            writer.close();
            System.out.println("Output written to " + outputFile);
        } /*
         * public class CustomComparator implements
         * Comparator<ProteinDetectionHypothesis> { //@Override public int
         * compare(ProteinDetectionHypothesis pdh1, ProteinDetectionHypothesis
         * pdh2) { return
         * pdh1.getPeptideHypothesis().size().compareTo(pdh2.getPeptideHypothesis().size());
         *
         * }
         * }
         */ catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } finally {
            try {
                writer.close();
            } catch (IOException ex) {
                String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
                String className = this.getClass().getName();
                String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                        + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
                System.out.println(message);
            }
        }

    }
    /*
     * public class CustomComparator implements
     * Comparator<ProteinDetectionHypothesis> { //@Override public int
     * compare(ProteinDetectionHypothesis pdh1, ProteinDetectionHypothesis pdh2)
     * { return
     * pdh1.getPeptideHypothesis().size().compareTo(pdh2.getPeptideHypothesis().size());
     *
     * }
     * }
     */

    private void applyProteoAnnotator() {
        Map<String, Set<String>> distinctPeptides = new HashMap<>();
        for (ProteinAmbiguityGroup proteinAmbiguityGroup : proteinAmbiguityGroupList) {
            List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
            for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
                List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesisList.get(i).getPeptideHypothesis();
                for (PeptideHypothesis peptideHypothesis : peptideHypothesisList) {
                    String peptideEvidenceRef = peptideHypothesis.getPeptideEvidenceRef();
                    PeptideEvidence pe = peptideEvidenceIdHashMap.get(peptideEvidenceRef);
                    String peptideRef = pe.getPeptideRef();
                    Peptide peptide = peptideIdHashMap.get(peptideRef);
                    String dbSeqRef = pe.getDBSequenceRef();
                    DBSequence dBSequence = dbSequenceIDToDBSequenceHashMap.get(dbSeqRef);
                    String accession = dBSequence.getAccession();
                    accession = accession.split("\\|")[1];
                    String newValue = accession.substring(0, 1);
                    if (!pe.isIsDecoy()) {
                        newValue += "_NotDecoy";
                    } else {
                        newValue += "_Decoy";
                    }
                    Set<String> value = distinctPeptides.get(peptide.getPeptideSequence());
                    if (value != null) {
                        if (!value.contains(newValue)) {
                            value.add(newValue);
                        }
                    } else {
                        Set<String> hs = new HashSet<>();
                        hs.add(newValue);
                        distinctPeptides.put(peptide.getPeptideSequence(), hs);
                    }

                }
            }
        }

        for (ProteinAmbiguityGroup proteinAmbiguityGroup : proteinAmbiguityGroupList) {

            List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
            ProteinDetectionHypothesis anchorProteinDetectionHypothesis = null;
            for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
                ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
                List<CvParam> cvParamListproteinDetectionHypothesis = proteinDetectionHypothesis.getCvParam();
                for (int s = 0; s < cvParamListproteinDetectionHypothesis.size(); s++) {
                    CvParam cvParam = cvParamListproteinDetectionHypothesis.get(s);
                    String accession = cvParam.getAccession();
                    if (accession.equals("MS:1001591")) {
                        anchorProteinDetectionHypothesis = proteinDetectionHypothesis;
                        break;
                    }
                }
            }
            if (anchorProteinDetectionHypothesis != null) {
                List<UserParam> userParamList = anchorProteinDetectionHypothesis.getUserParam();
                Set<String> razorAndUniquePeptide = new HashSet();
                String razor = "";
                String unique = "";
                for (int i = 0; i < userParamList.size(); i++) {
                    UserParam userParam = userParamList.get(i);
                    if (userParam.getName().equals("unique peptides")) {
                        unique = userParam.getValue();
                    }
                    if (userParam.getName().equals("razor peptides")) {
                        razor = userParam.getValue();
                    }
                }
                if (!razor.equals("")) {
                    razorAndUniquePeptide.addAll(Arrays.asList(razor.split("; ")));
                }
                if (!unique.equals("")) {
                    razorAndUniquePeptide.addAll(Arrays.asList(unique.split("; ")));
                }
                int countNonA = 0;
                double scoreNonA = 0;
                String nonAPeptide = "";

                for (Entry<String, Set<String>> entry : distinctPeptides.entrySet()) {
                    String peptideSequence = entry.getKey();
                    Set<String> value = entry.getValue();
                    if (!value.contains("A_NotDecoy") && razorAndUniquePeptide.contains(peptideSequence)) {
                        countNonA = countNonA + value.size();
                        scoreNonA = scoreNonA + pepSeqToBestScoreMap.get(peptideSequence).doubleValue();
                        nonAPeptide = nonAPeptide + peptideSequence + ";";
                    }
                }
                proteinAmbiguityGroup.getUserParam().add(MzidLibUtils.makeUserParam("nonAPeptide", nonAPeptide));
                proteinAmbiguityGroup.getCvParam().add(MzidLibUtils.makeCvParam("MS:1002474", "ProteoAnnotator:non-canonical gene model score", psiCV, String.valueOf(scoreNonA)));
                proteinAmbiguityGroup.getCvParam().add(MzidLibUtils.makeCvParam("MS:1002475", "ProteoAnnotator:count alternative peptides", psiCV, String.valueOf(countNonA)));

            }

        }
    }
}
