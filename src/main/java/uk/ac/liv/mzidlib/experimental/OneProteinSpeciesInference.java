package uk.ac.liv.mzidlib.experimental;

import java.io.File;
import java.io.FileOutputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBException;

import uk.ac.ebi.JmzIdentMLParser;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;


/*
 * @author Fawaz Ghali, University of Liverpool, 2012
 */
public class OneProteinSpeciesInference {

    private MzIdentMLUnmarshaller unmarshaller;
    private List<String> peptideList = new ArrayList<>();
    private List<String> dbSequenceList = new ArrayList<>();
    private List<String> matrix = new ArrayList<>();
    private List<String> matrix_sii = new ArrayList<>();
    private List<String> matrix_pe = new ArrayList<>();
    private Map<String, PeptideEvidence> peptideEvidenceIdHashMap = new HashMap<>();
    private Map<String, Peptide> peptideIdHashMap = new HashMap<>();
    private Map<String, List<ProteinDetectionHypothesis>> peptide_pdh_HashMap = new HashMap<>();
    private Map<String, String> pdhId_description_HashMap = new HashMap<>();
    private Map<String, String> dbseqId_taxon_HashMap = new HashMap<>();
    private List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = new ArrayList<>();
    private List<ProteinAmbiguityGroup> proteinAmbiguityGroupList = new ArrayList<>();
    ProteinDetectionList proteinDetectionList = new ProteinDetectionList();
    private String score = "X\\!Tandem:hyperscore";
    private String transform = "false";
    private double proteinScore = 0;

    private static double TEMP_SCORE_CUTOFF = 40;           //Testing thresholds for removing low quality identifications

    private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("F012436.mzid");

    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("F012436_LM_NERC.mzid");
    private int speciesNameCount = 0;

    private void init() {

        try {

            Iterator<PeptideEvidence> iterPeptideEvidence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence peptideEvidence = iterPeptideEvidence.next();
                peptideEvidenceIdHashMap.put(peptideEvidence.getId(), peptideEvidence);

            }

            Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
                peptideIdHashMap.put(peptide.getId(), peptide);

            }

            Iterator<SpectrumIdentificationItem> iterSpectrumIdentificationItem = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);
            while (iterSpectrumIdentificationItem.hasNext()) {
                SpectrumIdentificationItem spectrumIdentificationItem = iterSpectrumIdentificationItem.next();
                if (spectrumIdentificationItem.isPassThreshold()) {
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                    for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
                        PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());
                        if (!peptideList.contains(peptideEvidence.getPeptideRef())) {
                            peptideList.add(peptideEvidence.getPeptideRef());
                        }
                        if (!dbSequenceList.contains(peptideEvidence.getDBSequenceRef())) {
                            dbSequenceList.add(peptideEvidence.getDBSequenceRef());
                        }
                        if (!matrix.contains(peptideEvidence.getDBSequenceRef() + peptideEvidence.getPeptideRef())) {
                            matrix.add(peptideEvidence.getDBSequenceRef() + peptideEvidence.getPeptideRef());
                            matrix_sii.add(spectrumIdentificationItem.getId());
                            matrix_pe.add(peptideEvidenceRef.getPeptideEvidenceRef());

                        }
                    }
                }
            }
            extractSpeciesName();

        } catch (Exception e) {
            System.out.println("At init()" + e.getMessage());
            e.printStackTrace();
        }
    }

    private void algorithm() {
        buildPDHs();
        buildScore();
        //buildScoreUniquePepSequences();
        buildPAGs();
        writeFdrValues();
        // test Matrix
//        for (int i = 0; i < peptideList.size(); i++) {
//            String peptide = peptideList.get(i);
//            for (int j = 0; j < dbSequenceList.size(); j++) {
//                String dbSequence = dbSequenceList.get(j);
//                if (matrix.contains(dbSequence + peptide)) {
//                    System.out.print(matrix_sii.get(matrix.indexOf(dbSequence + peptide)) + "\t");
//                } else {
//                    System.out.print("-\t");
//                }
//            }
//            System.out.print("\n");
//        }
        // test PAGs
//        for (int i = 0; i < proteinAmbiguityGroupList.size(); i++) {
//            System.out.println("PAG: " + i);
//            ProteinAmbiguityGroup proteinAmbiguityGroup = proteinAmbiguityGroupList.get(i);
//
//            List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = proteinAmbiguityGroup.getProteinDetectionHypothesis();
//            for (int j = 0; j < proteinDetectionHypothesisList.size(); j++) {
//                System.out.println("PDH: " + proteinDetectionHypothesisList.get(j).getId());
//                ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(j);
//                List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesis.getPeptideHypothesis();
//                for (int k = 0; k < peptideHypothesisList.size(); k++) {

//                   PeptideHypothesis ph = peptideHypothesisList.get(k);
//                    System.out.println(peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef());
//                    System.out.println(matrix_sii.get(matrix.indexOf(proteinDetectionHypothesis.getDBSequenceRef() + peptideEvidenceIdHashMap.get(ph.getPeptideEvidenceRef()).getPeptideRef())) + "\t");
//                }
//            }
//            System.out.println("==========================");
//        }
    }

    public static void main(String[] args) {
        try {
            OneProteinSpeciesInference proteinInference = new OneProteinSpeciesInference();
            File myFile = new File("C:\\Work\\PSI\\mzIdentML\\mzidLib\\svn\\trunk\\resources\\sheep_identifier_fix.mzid");

            if (args != null && args.length == 3) {
                //File myFile = new File("C:\\uni\\BIOL 703\\Project Work\\Species Identifier\\sheep_identifier.mzid");
                proteinInference.unmarshaller = new MzIdentMLUnmarshaller(myFile);
                // proteinInference.unmarshaller = new MzIdentMLUnmarshaller(file);
                proteinInference.score = args[1];
                proteinInference.transform = args[2];
            } else {
                proteinInference.unmarshaller = new MzIdentMLUnmarshaller(myFile);
            }
            proteinInference.init();
            proteinInference.algorithm();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void buildPDHs() {
        for (int i = 0; i < dbSequenceList.size(); i++) {
            String dbSequence = dbSequenceList.get(i);
            //ProteinDetectionHypothesis
            ProteinDetectionHypothesis proteinDetectionHypothesis = new ProteinDetectionHypothesis();
            proteinDetectionHypothesis.setId("PDH_" + i);
            proteinDetectionHypothesis.setPassThreshold(true);
            DBSequence dbSeq = new DBSequence();
            dbSeq.setId(dbSequence);
            proteinDetectionHypothesis.setDBSequence(dbSeq);
            for (int j = 0; j < peptideList.size(); j++) {
                String peptide = peptideList.get(j);
                if (matrix.contains(dbSequence + peptide)) {

                    // SpectrumIdentificationItemRef
                    SpectrumIdentificationItemRef spectrumIdentificationItemRef = new SpectrumIdentificationItemRef();
                    spectrumIdentificationItemRef.setSpectrumIdentificationItemRef(matrix_sii.get(matrix.indexOf(dbSequence + peptide)));
                    //PeptideHypothesis
                    PeptideHypothesis peptideHypothesis = new PeptideHypothesis();
                    PeptideEvidence peptideEvidence = new PeptideEvidence();
                    peptideEvidence.setId(matrix_pe.get(matrix.indexOf(dbSequence + peptide)));
                    peptideHypothesis.setPeptideEvidence(peptideEvidence);

                    peptideHypothesis.getSpectrumIdentificationItemRef().add(spectrumIdentificationItemRef);
                    proteinDetectionHypothesis.getPeptideHypothesis().add(peptideHypothesis);
                    List<ProteinDetectionHypothesis> list = peptide_pdh_HashMap.get(peptide);
                    if (list == null) {
                        list = new ArrayList<>();
                        peptide_pdh_HashMap.put(peptide, list);
                    }
                    list.add(proteinDetectionHypothesis);

                }

            }

            proteinDetectionHypothesisList.add(proteinDetectionHypothesis);

        }

    }

    private void buildPAGs() {
        ProteinAmbiguityGroup proteinAmbiguityGroup = new ProteinAmbiguityGroup();
        proteinAmbiguityGroup.setId("PAG_1");
        proteinAmbiguityGroup.getProteinDetectionHypothesis().addAll(proteinDetectionHypothesisList);
        proteinDetectionList.getProteinAmbiguityGroup().add(proteinAmbiguityGroup);
    }

    private void buildScore() {

        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
            System.out.println("Prot" + proteinDetectionHypothesis.getDBSequenceRef());
            proteinDetectionHypothesis.setPassThreshold(false);
            List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesis.getPeptideHypothesis();
            //  Loop through all PeptideHypothesis elements
            for (int j = 0; j < peptideHypothesisList.size(); j++) {
                // Retrieve the best SII score
                double bestScorePerPeptideHypothesis = 0;
                SpectrumIdentificationItem bestSii = null;
                PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(j);
                List<SpectrumIdentificationItemRef> spectrumIdentificationItemRefList = peptideHypothesis.getSpectrumIdentificationItemRef();
                for (int k = 0; k < spectrumIdentificationItemRefList.size(); k++) {
                    try {
                        SpectrumIdentificationItemRef spectrumIdentificationItemRef = spectrumIdentificationItemRefList.get(k);
                        SpectrumIdentificationItem spectrumIdentificationItem = unmarshaller.unmarshal(SpectrumIdentificationItem.class, spectrumIdentificationItemRef.getSpectrumIdentificationItemRef());
                        List<CvParam> cvParamList = spectrumIdentificationItem.getCvParam();
                        for (int l = 0; l < cvParamList.size(); l++) {
                            CvParam cvParam = cvParamList.get(l);
                            if (cvParam.getName().equals(score)) {
                                double score = Double.valueOf(cvParam.getValue()).doubleValue();
                                if (transform.equals("true")) {
                                    score = Math.log10(score);
                                }
                                if (score > bestScorePerPeptideHypothesis) {
                                    bestScorePerPeptideHypothesis = score;
                                    bestSii = spectrumIdentificationItem;
//                                    topRankingPDH = proteinDetectionHypothesis;
                                }

                            }

                        }
                    } catch (JAXBException ex) {

                        ex.printStackTrace();
                    }
                }
                if (bestSii != null) // Divide this score by the number of proteins to which it is matched
                {
//                    bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / bestSii.getPeptideEvidenceRef().size();

                    //bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / speciesNameCount;
                    //Now we need to divide the score by the number of different taxa to which it is matched
                    List<PeptideEvidenceRef> peptideEvidenceRefList = bestSii.getPeptideEvidenceRef();
                    Map<String, Boolean> countTaxaMap = new HashMap<>();
                    for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);
                        PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

                        String taxon = dbseqId_taxon_HashMap.get(peptideEvidence.getDBSequenceRef());
                        if (taxon != null) {
                            countTaxaMap.put(taxon, true);
                        } else {
                            System.out.println("Taxon not recognized");
                        }
                    }

                    //System.out.println("Best SII: " + bestSii.getId() + " score: " + bestScorePerPeptideHypothesis + " / " + countTaxaMap.size() + " = " + bestScorePerPeptideHypothesis / countTaxaMap.size());
                    //if(countTaxaMap.size() == 1){
                    //bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis;          //Nothing to do here
                    System.out.println("\tScoring pep " + bestSii.getId() + " " + bestSii.getPeptideRef() + " score: " + bestScorePerPeptideHypothesis);

                    if (bestSii.getPeptideRef().contains("x")) {

                        bestScorePerPeptideHypothesis = 0;  //remove peps containing x's
                    }
                    //}
                    //else{
                    //    bestScorePerPeptideHypothesis = 0; 
                    //}

                    // bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / countTaxaMap.size();   
                }
                // Add up the PHscores from all PeptideHypothesis elements to arrive at a final â€œprotein scoreâ€�
                proteinScore = proteinScore + bestScorePerPeptideHypothesis;

            }
            CvParam cvParam = new CvParam();
            cvParam.setAccession("MS:1001171");
            cvParam.setName(score);
            cvParam.setValue("" + proteinScore);
            proteinDetectionHypothesis.getCvParam().add(cvParam);

//            System.out.println("PDH "+proteinDetectionHypothesis.getId() + " has score: "+ proteinScore);
            proteinScore = 0;

        }
        // set only this one to â€œpassThreshold = trueâ€�
        // Find the top-ranking PDH
        ProteinDetectionHypothesis topRankingPDH = null;
        double topPDH = 0;
//        System.out.println(proteinDetectionHypothesisList.size());
        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
            double score = Double.valueOf(proteinDetectionHypothesis.getCvParam().get(0).getValue());
            if (score > topPDH) {
                topPDH = score;
                topRankingPDH = proteinDetectionHypothesis;
            }
        }
        topRankingPDH.setPassThreshold(true);

        Iterator<SpectraData> spectraDataList = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectraData);
        while (spectraDataList.hasNext()) {
            SpectraData specData = spectraDataList.next();
            int filePos = specData.getLocation().lastIndexOf("/");
            String file = specData.getLocation().substring(filePos + 1);
            System.out.print(file + "\t");
        }

        //int filePos = xmlFileURL.toString().lastIndexOf("/");
        //String file = xmlFileURL.toString().substring(filePos+1);
        System.out.println(topRankingPDH.getDBSequenceRef() + "\t" + topRankingPDH.getCvParam().get(0).getValue() + "\t" + pdhId_description_HashMap.get(topRankingPDH.getDBSequenceRef()));

    }

    /*
     * 
     * ARJ 25/06/2012 Adapted this method, since peptides with modifications in multiple positions were being scored separately
     */
    private void buildScoreUniquePepSequences() {

        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);

            proteinDetectionHypothesis.setPassThreshold(false);
            List<PeptideHypothesis> peptideHypothesisList = proteinDetectionHypothesis.getPeptideHypothesis();
            //  Loop through all PeptideHypothesis elements

            Map<String, Double> pepSeqToBestScore = new HashMap<>();

            for (int j = 0; j < peptideHypothesisList.size(); j++) {
                // Retrieve the best SII score
                double bestScorePerPeptideHypothesis = 0;
                SpectrumIdentificationItem bestSii = null;
                PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(j);
                List<SpectrumIdentificationItemRef> spectrumIdentificationItemRefList = peptideHypothesis.getSpectrumIdentificationItemRef();
                for (int k = 0; k < spectrumIdentificationItemRefList.size(); k++) {
                    try {
                        SpectrumIdentificationItemRef spectrumIdentificationItemRef = spectrumIdentificationItemRefList.get(k);
                        SpectrumIdentificationItem spectrumIdentificationItem = unmarshaller.unmarshal(SpectrumIdentificationItem.class, spectrumIdentificationItemRef.getSpectrumIdentificationItemRef());
                        List<CvParam> cvParamList = spectrumIdentificationItem.getCvParam();
                        for (int l = 0; l < cvParamList.size(); l++) {
                            CvParam cvParam = cvParamList.get(l);
                            if (cvParam.getName().equals(score)) {
                                double score = Double.valueOf(cvParam.getValue());
                                if (transform.equals("true")) {
                                    score = Math.log10(score);
                                }
                                if (score > bestScorePerPeptideHypothesis) {
                                    bestScorePerPeptideHypothesis = score;
                                    bestSii = spectrumIdentificationItem;
//                                    topRankingPDH = proteinDetectionHypothesis;
                                }

                            }

                        }
                    } catch (JAXBException ex) {

                        ex.printStackTrace();
                    }
                }
                if (bestSii != null) // Divide this score by the number of proteins to which it is matched
                {
//                    bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / bestSii.getPeptideEvidenceRef().size();

                    //bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / speciesNameCount;
                    //Now we need to divide the score by the number of different taxa to which it is matched
                    List<PeptideEvidenceRef> peptideEvidenceRefList = bestSii.getPeptideEvidenceRef();
                    Map<String, Boolean> countTaxaMap = new HashMap<>();
                    for (int k = 0; k < peptideEvidenceRefList.size(); k++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(k);
                        PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

                        String taxon = dbseqId_taxon_HashMap.get(peptideEvidence.getDBSequenceRef());
                        if (taxon != null) {
                            countTaxaMap.put(taxon, true);
                        } else {
                            System.out.println("Taxon not recognized");
                        }
                    }

                    //System.out.println("Best SII: " + bestSii.getId() + " score: " + bestScorePerPeptideHypothesis + " / " + countTaxaMap.size() + " = " + bestScorePerPeptideHypothesis / countTaxaMap.size());
                    if (countTaxaMap.size() == 1) {
                        //bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis;          //Nothing to do here
                        //System.out.println("\tScoring pep " + bestSii.getId() + " " + bestSii.getPeptideRef() + " score: " + bestScorePerPeptideHypothesis);

                        if (bestSii.getPeptideRef().contains("X")) {

                            bestScorePerPeptideHypothesis = 0;  //remove peps containing x's
                        }
                    } else {
                        bestScorePerPeptideHypothesis = 0;
                    }

                    // bestScorePerPeptideHypothesis = bestScorePerPeptideHypothesis / countTaxaMap.size();   
                    Peptide pep = peptideIdHashMap.get(bestSii.getPeptideRef());
                    String pepSeq = pep.getPeptideSequence();
                    //System.out.println("\tScoring pep:" + pep.getPeptideSequence());

                    double currentBestScoreForPep = 0.0;
                    if (pepSeqToBestScore.get(pepSeq) != null) {
                        currentBestScoreForPep = pepSeqToBestScore.get(pepSeq);
                    }

                    if (bestScorePerPeptideHypothesis > currentBestScoreForPep) {
                        pepSeqToBestScore.put(pep.getPeptideSequence(), bestScorePerPeptideHypothesis);
                    }

                }
                // Add up the PHscores from all PeptideHypothesis elements to arrive at a final â€œprotein scoreâ€�
                //proteinScore = proteinScore + bestScorePerPeptideHypothesis;

            }

            for (String pepSeq : pepSeqToBestScore.keySet()) {
                if (pepSeqToBestScore.get(pepSeq) > TEMP_SCORE_CUTOFF) {
                    proteinScore = proteinScore + pepSeqToBestScore.get(pepSeq);
                }
            }

            //Print some output for debugging purposes
            if (proteinScore > 0) {
                //System.out.println("Prot" + proteinDetectionHypothesis.getDBSequenceRef());
                for (String pepSeq : pepSeqToBestScore.keySet()) {
                    if (pepSeqToBestScore.get(pepSeq) > TEMP_SCORE_CUTOFF) {
                        System.out.println(proteinDetectionHypothesis.getDBSequenceRef() + "\t" + pepSeq + "\t" + pepSeqToBestScore.get(pepSeq) + "\t" + proteinScore);
                    }
                }
            }

            CvParam cvParam = new CvParam();
            cvParam.setAccession("MS:1001171");
            cvParam.setName(score);
            cvParam.setValue("" + proteinScore);
            proteinDetectionHypothesis.getCvParam().add(cvParam);

//            System.out.println("PDH "+proteinDetectionHypothesis.getId() + " has score: "+ proteinScore);
            proteinScore = 0;

        }
        // set only this one to â€œpassThreshold = trueâ€�
        // Find the top-ranking PDH
        ProteinDetectionHypothesis topRankingPDH = null;
        double topPDH = 0;
//        System.out.println(proteinDetectionHypothesisList.size());
        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
            double score = Double.valueOf(proteinDetectionHypothesis.getCvParam().get(0).getValue());
            if (score > topPDH) {
                topPDH = score;
                topRankingPDH = proteinDetectionHypothesis;
            }
        }
        topRankingPDH.setPassThreshold(true);

        Iterator<SpectraData> spectraDataList = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectraData);
        while (spectraDataList.hasNext()) {
            SpectraData specData = spectraDataList.next();
            int filePos = specData.getLocation().lastIndexOf("/");
            String file = specData.getLocation().substring(filePos + 1);
            System.out.print(file + "\t");
        }

        //int filePos = xmlFileURL.toString().lastIndexOf("/");
        //String file = xmlFileURL.toString().substring(filePos+1);
        System.out.println(topRankingPDH.getDBSequenceRef() + "\t" + topRankingPDH.getCvParam().get(0).getValue() + "\t" + pdhId_description_HashMap.get(topRankingPDH.getDBSequenceRef()));

    }

    // Write values
    private void writeFdrValues() {
        try {
            MzIdentML mzIdentML = unmarshaller.unmarshal(MzIdentMLElement.MzIdentML);
            proteinDetectionList.setId("PDL_1");
//            System.out.println("Number of PDHs before: "+proteinDetectionList.getProteinAmbiguityGroup().get(0).getProteinDetectionHypothesis().size());
            mzIdentML.getDataCollection().getAnalysisData().setProteinDetectionList(proteinDetectionList);
//            System.out.println("Number of PDHs: "+mzIdentML.getDataCollection().getAnalysisData().getProteinDetectionList().getProteinAmbiguityGroup().get(0).getProteinDetectionHypothesis().size());
            MzIdentMLMarshaller m = new MzIdentMLMarshaller();
            File myFile = new File("C:\\Work\\PSI\\mzIdentML\\mzidLib\\svn\\trunk\\resources\\sheep_identifier_out.mzid");
            m.marshal(mzIdentML, new FileOutputStream(myFile));

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
            ex.printStackTrace();
        }
    }

    private void extractSpeciesName() {
        Iterator<DBSequence> iter = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        List<String> speciesNameList = new ArrayList<>();
        while (iter.hasNext()) {
            DBSequence result = iter.next();
            System.out.println("DBSequence: " + result.getId());
            String id = result.getId();
            String speciesName = id.substring(id.indexOf("[") + 1, id.indexOf("]"));
            System.out.println("taxon: " + speciesName);

            dbseqId_taxon_HashMap.put(result.getId(), speciesName);

            /*
             List<CvParam> cvParamList = result.getCvParam();
             for (int i = 0; i < cvParamList.size(); i++) {
             CvParam cvParam = cvParamList.get(i);
             if(cvParam.getName().equals("protein description")){
             String value = cvParam.getValue();

             pdhId_description_HashMap.put(result.getId(),value);
             if(value!=null){
             int start = value.lastIndexOf("taxon=");
             if(start>0){
             String speciesName = value.substring(start);                        
             dbseqId_taxon_HashMap.put(result.getId(),speciesName);
             if(!speciesNameList.contains(speciesName)){
             speciesNameList.add(speciesName);
                            
             //System.out.println("speciesName: " + speciesName);
             }
             }
             }
             } 
             }
             */
        }
        speciesNameCount = speciesNameList.size();
        //System.out.println("species count: " + speciesNameCount);

    }
}
