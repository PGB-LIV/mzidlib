package uk.ac.liv.mzidlib.converters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author jonesar
 */
public class MzIdentMLToCSV {

    private MzIdentMLUnmarshaller unmarshaller;
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("CPTAC_Progenesis_Identifications.mzid");
    //private URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("55merge_mascot_full.mzid");
    private String xmlFile = "resources/55merge_mascot_full.mzid";
    private List<String> peptideList = new ArrayList<>();
    private List<String> dbSequenceList = new ArrayList<>();
    private List<String> matrix = new ArrayList<>();
    private List<String> matrix_sii = new ArrayList<>();
    private List<String> matrix_pe = new ArrayList<>();
    private Map<String, PeptideEvidence> peptideEvidenceIdHashMap = new HashMap<>();
    private Map<String, Peptide> peptideIdHashMap = new HashMap<>();
    private Map<String, SpectraData> spectraDataIdHashMap = new HashMap<>();
    private Map<String, SpectrumIdentificationItem> siiIdHashMap = new HashMap<>();
    private Map<String, SpectrumIdentificationResult> siiIdToSirHashMap = new HashMap<>();
    private Map<Integer, String> columnToScoreMap = new HashMap<>();
    private Map<Integer, String> columnToProtScoreMap = new HashMap<>();
    private Map<String, DBSequence> dbSequenceIdHashMap = new HashMap<>();
    private Map<String, ProteinDetectionHypothesis> pdhIdHashMap = new HashMap<>();
    private Map<String, List<ProteinDetectionHypothesis>> peptide_pdh_HashMap = new HashMap<>();
    private List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = new ArrayList<>();
    private List<ProteinAmbiguityGroup> proteinAmbiguityGroupList = new ArrayList<>();
    private ProteinDetectionList proteinDetectionList = new ProteinDetectionList();
    private String sep = ",";
    private String pagHeader = "PAG ID" + sep + "PAG score" + sep + "protein accession" + sep + "Pass Threshold (Protein)" + sep + "description" + sep + "group membership" + sep;
    private String spectrumHeader = "Raw data location" + sep + "Spectrum ID" + sep + "Spectrum Title" + sep + "Retention Time (s)" + sep;
    private String psmHeader = "PSM_ID" + sep + "rank" + sep + "Pass Threshold" + sep + "Calc m/z" + sep + "Exp m/z" + sep + "Charge" + sep + "Sequence" + sep + "Modifications";
    private String pScoreHeader = "";     //Protein score header will be set only after reading the file
    private String scoreHeader = "";     //This will be set only after reading the file
    // Added by Fawaz Ghali 13/05/2014 exportProteoAnnotator
    private String exportProteoAnnotatorHeader = "countNonA" + sep + "scoreNonA" + sep + "nonAPeptide" + sep + "A genes" + sep + "protein group-level q-value";
    private String endPsmHeader = sep + "proteinacc_start_stop_pre_post_;" + sep + "Is decoy";
    private String representativeProteinAcc = "MS:1001591";     //Used to identify the representative of each group - only used for the one line export of PAGS
    private boolean isVerbose = true;

    public static void main(String[] args) {
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();

        //TODO - Undecided which if any command line arguments to include - minimally need to know whether to export Peptides or PAGs
        if (args != null && args.length == 3) {
            mzidToCsv.unmarshaller = new MzIdentMLUnmarshaller(new File(args[0]));
            mzidToCsv.init(args[1], args[2]);

        } else {

            mzidToCsv.unmarshaller = new MzIdentMLUnmarshaller(new File(mzidToCsv.xmlFile));
            mzidToCsv.init("out.csv", "exportPSMs");
            //mzidToCsv.init("out.csv","exportProteinGroups");

            //System.out.println("Error - correct usage MzIdentMLToCSV inputFile.mzid outputFile.csv [exportProteinGroups|exportPSMs|exportProteinsOnly]");
        }

    }

    // FG: this function is used to call this class from outside the library
    public void useMzIdentMLToCSV(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, String out, String option, boolean verbose) {
        isVerbose = verbose;
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        mzidToCsv.unmarshaller = mzIdentMLUnmarshaller;
        mzidToCsv.init(out, option);

    }

    // FG: this function is used to call this class from outside the library
    public void useMzIdentMLToCSV(String mzid, String out, String option, boolean verbose) {
        isVerbose = verbose;
        MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid));
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        mzidToCsv.unmarshaller = mzIdentMLUnmarshaller;
        mzidToCsv.init(out, option);

    }

    private void init(String outputFile, String exportOption) {
        Writer out = null;
        try {
            out = new BufferedWriter(new FileWriter(outputFile));
            //Read all the objects we will need into hashes that are not automatically resolved by object reference
            if (isVerbose) {
                System.out.print("About to iterate over PepEvid...");
            }
            Iterator<PeptideEvidence> iterPeptideEvidence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence peptideEvidence = iterPeptideEvidence.next();
                peptideEvidenceIdHashMap.put(peptideEvidence.getId(), peptideEvidence);
            }
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over Peptide");
            }
            Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
                peptideIdHashMap.put(peptide.getId(), peptide);
            }
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over Spectra Data");
            }
            Iterator<SpectraData> iterSpectraData = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectraData);
            while (iterSpectraData.hasNext()) {
                SpectraData spectraData = iterSpectraData.next();
                spectraDataIdHashMap.put(spectraData.getId(), spectraData);
            }
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over DBsequence");
            }
            Iterator<DBSequence> iterDBSequence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {
                DBSequence dbSequence = iterDBSequence.next();
                dbSequenceIdHashMap.put(dbSequence.getId(), dbSequence);
            }
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over PDH");
            }
            Iterator<ProteinDetectionHypothesis> iterPDH = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinDetectionHypothesis);
            Integer pCounter = 0;
            while (iterPDH.hasNext()) {
                ProteinDetectionHypothesis pdh = iterPDH.next();
                pdhIdHashMap.put(pdh.getId(), pdh);

                for (CvParam cvParam : pdh.getCvParam()) {
                    if (cvParam.getAccession().equals("MS:1001591") || cvParam.getAccession().equals("MS:1001592") || cvParam.getAccession().equals("MS:1001593")
                            || cvParam.getAccession().equals("MS:1001594") || cvParam.getAccession().equals("MS:1001595") || cvParam.getAccession().equals("MS:1001596")
                            || cvParam.getAccession().equals("MS:1001597")
                            || cvParam.getAccession().equals("MS:1001598")
                            || cvParam.getAccession().equals("MS:1001599")) {        //do nothing - these are specifically handled
                        //ToDO this code could be improved using an array of values...
                    } else if (cvParam.getValue() != null) {
                        if (!columnToProtScoreMap.containsValue(cvParam.getName())) {
                            columnToProtScoreMap.put(pCounter, cvParam.getName());
                            pCounter++;
                        }
                    }
                }

                for (UserParam userParam : pdh.getUserParam()) {
                    if (!columnToProtScoreMap.containsValue(userParam.getName())) {
                        columnToProtScoreMap.put(pCounter, userParam.getName());
                        pCounter++;
                    }

                }

            }
            for (int i = 0; i < pCounter; i++) {
                pScoreHeader += columnToProtScoreMap.get(i) + sep;
            }
            //Now let's see what scores we have in the file
            //TODO - I'm not sure this is the fastest way to parse the files; these are unmarshalled again below - inefficient?
            //Iterator<SpectrumIdentificationItem> iterSII = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);
            Integer counter = 0;
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to iterate over SIR");
            }
            Iterator<SpectrumIdentificationResult> iterSIR = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            List<SpectrumIdentificationResult> sirList = new ArrayList<>();
            while (iterSIR.hasNext()) {
                SpectrumIdentificationResult sir = iterSIR.next();
                sirList.add(sir);

                List<SpectrumIdentificationItem> listSII = sir.getSpectrumIdentificationItem();

                for (SpectrumIdentificationItem sii : listSII) {
                    siiIdHashMap.put(sii.getId(), sii);
                    siiIdToSirHashMap.put(sii.getId(), sir);
                    for (CvParam cvParam : sii.getCvParam()) {
                        if (cvParam.getValue() != null) {
                            if (!columnToScoreMap.containsValue(cvParam.getName())) {
                                columnToScoreMap.put(counter, cvParam.getName());
                                counter++;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < counter; i++) {
                scoreHeader += sep + columnToScoreMap.get(i);
            }
            if (isVerbose) {
                System.out.println("...done");
                System.out.print("About to create output");
            }
            if (exportOption.equals("exportPSMs")) {

                out.write(spectrumHeader + psmHeader + scoreHeader);
                out.write(endPsmHeader + "\n");

                //Iterator<SpectrumIdentificationResult> iterSIR = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);                    
                for (SpectrumIdentificationResult sir : sirList) {

                    String sirLine = sirToString(sir);

                    List<SpectrumIdentificationItem> listSII = sir.getSpectrumIdentificationItem();

                    for (SpectrumIdentificationItem sii : listSII) {
                        out.write(sirLine + sep + siiToString(sii) + "\n");
                    }
                }

            } // Added by Fawaz Ghali 23/9/2015 to export peptides
            else if (exportOption.equals("exportPeptides")) {

                out.write(psmHeader + scoreHeader);
                out.write(endPsmHeader + "\n");

                Iterator<ProteinAmbiguityGroup> iterPAG = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                List<String> pepRefGroupIDList = new ArrayList();
                while (iterPAG.hasNext()) {
                    ProteinAmbiguityGroup pag = iterPAG.next();
                    List<CvParam> cvParamList = pag.getCvParam();
                    double qValue = 0;
                    for (int i = 0; i < cvParamList.size(); i++) {
                        CvParam cvParam = cvParamList.get(i);
                        if (cvParam.getAccession().equals("MS:1002373")) {
                            qValue = Double.valueOf(cvParam.getValue()).doubleValue();
                            break;
                        }

                    }
                    List<SpectrumIdentificationItem> siiList = new ArrayList();
                    if (qValue < 0.01) {
                        ProteinDetectionHypothesis repPdh = getRepresentativePDH(pag, representativeProteinAcc);
                        List<PeptideHypothesis> peptideHypothesisList = repPdh.getPeptideHypothesis();
                        for (int i = 0; i < peptideHypothesisList.size(); i++) {
                            PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(i);
                            List<SpectrumIdentificationItemRef> siiRefList = peptideHypothesis.getSpectrumIdentificationItemRef();
                            for (int j = 0; j < siiRefList.size(); j++) {
                                SpectrumIdentificationItemRef spectrumIdentificationItemRef = siiRefList.get(j);
                                SpectrumIdentificationItem sii = siiIdHashMap.get(spectrumIdentificationItemRef.getSpectrumIdentificationItemRef());
                                siiList.add(sii);

                            }
                        }
                        double bestScore = 0;
                        double siiqValue = 0;
                        String peptideRefGroupID = "";
                        SpectrumIdentificationItem bestSII;

                        for (int i = 0; i < siiList.size(); i++) {
                            SpectrumIdentificationItem sii = siiList.get(i);
                            double combinedFDR = 0;
                            List<CvParam> cvParamList1 = sii.getCvParam();
                            for (CvParam cvParam : cvParamList1) {
                                String accession = cvParam.getAccession();
                                if (accession.equals("MS:1002356")) {
                                    combinedFDR = Double.valueOf(cvParam.getValue()).doubleValue();
                                }
                                if (accession.equals("MS:1001868")) {
                                    siiqValue = Double.valueOf(cvParam.getValue()).doubleValue();
                                }
                                if (accession.equals("MS:1002520")) {
                                    peptideRefGroupID = cvParam.getValue();
                                }

                            }
                            if (combinedFDR < bestScore) {
                                bestScore = combinedFDR;
                                bestSII = sii;
                            }
                            if (siiqValue < 0.01 && !pepRefGroupIDList.contains(peptideRefGroupID)) {
                                pepRefGroupIDList.add(peptideRefGroupID);
                                out.write(siiToString(sii) + "\n");
                            }

                        }

                    }

                }

            } else if (exportOption.equals("exportProteinGroups")) {

                out.write(pagHeader);
                out.write(pScoreHeader);
                out.write(spectrumHeader + psmHeader + scoreHeader);
                out.write(endPsmHeader + "\n");

                Iterator<ProteinAmbiguityGroup> iterPAG = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterPAG.hasNext()) {
                    ProteinAmbiguityGroup pag = iterPAG.next();

                    String pagLine = pagToString(pag);

                    //handle PDHs
                    for (ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()) {
                        String pdhLine = pagLine;
                        pdhLine += pdhToString(pdh);

                        for (PeptideHypothesis pepH : pdh.getPeptideHypothesis()) {

                            List<SpectrumIdentificationItemRef> siiRefList = pepH.getSpectrumIdentificationItemRef();
                            for (SpectrumIdentificationItemRef siiRef : siiRefList) {
                                SpectrumIdentificationResult sir = siiIdToSirHashMap.get(siiRef.getSpectrumIdentificationItemRef());

                                SpectrumIdentificationItem sii = siiIdHashMap.get(siiRef.getSpectrumIdentificationItemRef());
                                out.write(pdhLine + sirToString(sir) + sep + siiToString(sii) + "\n");
                            }
                        }
                    }
                }

            } else if (exportOption.equals("exportRepProteinPerPAGOnly")) {
                out.write(pagHeader);
                out.write(pScoreHeader);
                out.write("\n");

                Iterator<ProteinAmbiguityGroup> iterPAG = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterPAG.hasNext()) {
                    ProteinAmbiguityGroup pag = iterPAG.next();
                    String pagLine = pagToString(pag);

                    ProteinDetectionHypothesis repPdh = getRepresentativePDH(pag, representativeProteinAcc);

                    String pdhLine = pagLine;

                    if (repPdh != null) {
                        pdhLine += pdhToString(repPdh);
                    }
                    out.write(pdhLine + "\n");
                }

            } else if (exportOption.equals("exportProteoAnnotator")) {     // Added by Fawaz Ghali 13/05/2014 exportProteoAnnotator

                out.write(pagHeader);
                out.write(pScoreHeader);
                // Added by Fawaz Ghali 13/05/2014 exportProteoAnnotator

                out.write(exportProteoAnnotatorHeader);
                out.write("\n");

                Iterator<ProteinAmbiguityGroup> iterPAG = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterPAG.hasNext()) {
                    ProteinAmbiguityGroup pag = iterPAG.next();
                    String pagLine = pagToString(pag);

                    ProteinDetectionHypothesis repPdh = getRepresentativePDH(pag, representativeProteinAcc);

                    String pdhLine = pagLine;

                    if (repPdh != null) {
                        pdhLine += pdhToString(repPdh);
                    }
                    // Added by Fawaz Ghali 13/05/2014 exportProteoAnnotator
                    String proteoAnnotatorLine = pdhLine;
                    proteoAnnotatorLine = proteoAnnotatorLine + proteoAnnotatorLineToString(pag);
                    out.write(proteoAnnotatorLine + "\n");
                }
            } else if (exportOption.equals("exportProteinsOnly")) {
                out.write(pagHeader);
                out.write(pScoreHeader);
                out.write("\n");

                Iterator<ProteinAmbiguityGroup> iterPAG = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterPAG.hasNext()) {
                    ProteinAmbiguityGroup pag = iterPAG.next();
                    String pagLine = pagToString(pag);

                    //handle PDHs
                    for (ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()) {
                        String pdhLine = pagLine;
                        pdhLine += pdhToString(pdh);
                        out.write(pdhLine + "\n");
                    }
                }
            } else {
                System.out.println("Error - correct usage MzIdentMLToCSV inputFile outputFile -exportType [exportProteinGroups|exportPSMs|exportProteinsOnly]");
            }
            out.close();
            System.out.println("Output written to " + outputFile);
        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } finally {
            try {
                out.close();
            } catch (IOException ex) {
                String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
                String className = this.getClass().getName();
                String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                        + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
                System.out.println(message);
            }
        }

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

    private String pagToString(ProteinAmbiguityGroup pag) {
        String pagString = pag.getId() + sep;

        if (pag.getCvParam().isEmpty()) {
            CvParam scoreParam = pag.getCvParam().get(0);
            pagString += scoreParam.getValue() + sep;
        } else {
            pagString += sep;
        }

        return pagString;

    }

    private String pdhToString(ProteinDetectionHypothesis pdh) {
        String pdhString = "\"" + dbSequenceIdHashMap.get(pdh.getDBSequenceRef()).getAccession() + "\"" + sep + pdh.isPassThreshold() + sep;

        DBSequence dbSeq = dbSequenceIdHashMap.get(pdh.getDBSequenceRef());
        String protDesc = "";

        String protGroupMembership = "";

        if (dbSeq != null) {
            for (CvParam cvParam : dbSeq.getCvParam()) {
                if (cvParam.getAccession().equals("MS:1001088")) {        //Protein description

                    String description = cvParam.getValue().replaceAll("\"", "");    //remove internal "
                    protDesc = "\"" + description + "\"";
                }
            }
        }

        Map<String, String> mapNameToValue = new HashMap<>();

        for (CvParam cvParam : pdh.getCvParam()) {
            if (cvParam.getAccession().equals("MS:1001591") || cvParam.getAccession().equals("MS:1001592") || cvParam.getAccession().equals("MS:1001593")
                    || cvParam.getAccession().equals("MS:1001594") || cvParam.getAccession().equals("MS:1001595") || cvParam.getAccession().equals("MS:1001596")
                    || cvParam.getAccession().equals("MS:1001597")
                    || cvParam.getAccession().equals("MS:1001598")
                    || cvParam.getAccession().equals("MS:1001599")) {        //Protein description
                protGroupMembership = "\"" + cvParam.getName();
                if (cvParam.getValue() != null) {
                    protGroupMembership += ":" + cvParam.getValue();
                }
                protGroupMembership += "\"";
            } else {
                mapNameToValue.put(cvParam.getName(), cvParam.getValue());
            }
        }

        for (UserParam userParam : pdh.getUserParam()) {

            mapNameToValue.put(userParam.getName(), userParam.getValue());

        }

        pdhString += protDesc + sep + protGroupMembership + sep;

        //Handle scores
        for (int i = 0; i < columnToProtScoreMap.size(); i++) {
            String score = columnToProtScoreMap.get(i);
            //System.out.println("test2" + score);
            if (mapNameToValue.containsKey(score)) {
                String scoreValue = mapNameToValue.get(score);
                //System.out.println("test3" + scoreValue);
                pdhString += scoreValue + sep;
            } else {
                pdhString += sep;
            }
        }

        return pdhString;

    }

    private String siiToString(SpectrumIdentificationItem sii) {
        String siiString = "";

        siiString = "\"" + sii.getId() + "\"" + sep + sii.getRank() + sep + sii.isPassThreshold() + sep + sii.getCalculatedMassToCharge() + sep + sii.getExperimentalMassToCharge() + sep + sii.getChargeState();
        Peptide pep = peptideIdHashMap.get(sii.getPeptideRef());    //get Peptide via the hash for this object
        siiString += sep + "\"" + pep.getPeptideSequence() + "\"";

        //Handle Mods
        siiString += sep;

        String modString = "";

        if (pep.getModification() != null) {
            int i = 0;
            for (Modification mod : pep.getModification()) {
                if (i > 0) {
                    modString += ";"; //Add an extra separator between mods                                        
                }
                modString += modToString(mod);
                i++;
            }
        }

        if (pep.getSubstitutionModification() != null) {
            int i = 0;
            for (SubstitutionModification subMod : pep.getSubstitutionModification()) {
                if (i > 0 || !modString.equals("")) {
                    modString += ";"; //Add an extra separator between mods                                        
                }
                modString += subModToString(subMod);
                i++;
            }
        }

        siiString += modString;

        Map<String, String> mapNameToValue = new HashMap<>();
        for (AbstractParam param : sii.getParamGroup()) {
            mapNameToValue.put(param.getName(), param.getValue());
            //System.out.println("test1" + param.getName() + "-> " + param.getValue());
        }

        //Handle scores
        for (int i = 0; i < columnToScoreMap.size(); i++) {
            String score = columnToScoreMap.get(i);
            //System.out.println("test2" + score);
            if (mapNameToValue.containsKey(score)) {
                String scoreValue = mapNameToValue.get(score);
                //System.out.println("test3" + scoreValue);
                siiString += sep + scoreValue;
            } else {
                siiString += sep;
            }
        }

        //Handle all protein maps
        siiString += sep + "\"";
        List<PeptideEvidenceRef> peptideEvidenceRefList = sii.getPeptideEvidenceRef();
        Boolean isDecoy = false;
        for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
            PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
            PeptideEvidence peptideEvidence = peptideEvidenceIdHashMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

            DBSequence dbSeq = dbSequenceIdHashMap.get(peptideEvidence.getDBSequenceRef());
            if (i > 0) {
                siiString += ";"; //Add an extra separator between mods                                        
            }
            siiString += dbSeq.getAccession() + "_" + peptideEvidence.getStart() + "_" + peptideEvidence.getEnd() + "_" + peptideEvidence.getPre() + "_" + peptideEvidence.getPost();
            if (peptideEvidence.isIsDecoy()) {
                isDecoy = true;
            }
        }
        siiString += "\"";

        siiString += sep + isDecoy;

        return siiString;

    }

    private String sirToString(SpectrumIdentificationResult sir) {
        String sirString = "";

        SpectraData spectraData = spectraDataIdHashMap.get(sir.getSpectraDataRef());
        sirString += spectraData.getLocation() + sep + "\"" + sir.getSpectrumID() + "\"";

        Double rtInSeconds = -1.0;
        String spectrumTitle = "";

        // <cvParam accession="MS:1001114" name="retention time(s)"  cvRef="PSI-MS" value="3488.676" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
        //  <cvParam accession="MS:1000796" name="spectrum title"  cvRef="PSI-MS" value="mam_050108o_CPTAC_study6_6E004.6805.6805.1" />
        //
        for (CvParam cvParam : sir.getCvParam()) {
            // Updated by FG: checking for old CV param 1114 or newer correct CV term 16.
            if (cvParam.getAccession().equals("MS:1001114") || cvParam.getAccession().equals("MS:1000016")) {
                if (cvParam.getUnitAccession().equals("UO:0000010")) {
                    rtInSeconds = Double.parseDouble(cvParam.getValue());
                } else if (cvParam.getUnitAccession().equals("UO:0000031")) {
                    rtInSeconds = Double.parseDouble(cvParam.getValue()) / 60;    //Convert minutes to seconds
                } else {
                    System.out.println("Error parsing RT - unit not recognised");
                }
            }

            if (cvParam.getAccession().equals("MS:1000796")) {
                spectrumTitle = cvParam.getValue();
            }
        }

        sirString += sep + "\"" + spectrumTitle + "\"" + sep + rtInSeconds;

        return sirString;

    }

    /*
     * Method to convert an mzid Mod element into a string of type
     * ModName@location
     */
    private String modToString(Modification mod) {

        String modString = "";

        if (mod.getCvParam() != null) {

            String modName = "";
            for (CvParam cvParam : mod.getCvParam()) {
                modName = cvParam.getName();
            }

            modString += modName;
        } else {

            if (mod.getMonoisotopicMassDelta() != null) {
                modString += mod.getMonoisotopicMassDelta();
            } else if (mod.getAvgMassDelta() != null) {
                modString += mod.getAvgMassDelta();
            }
        }

        if (mod.getLocation() != null) {
            modString += ":" + mod.getLocation();
        }
        return modString;

    }

    /*
     * Method to create and return a string representation of a substitution
     * modification
     */
    private String subModToString(SubstitutionModification subMod) {
        return subMod.getOriginalResidue() + "=>" + subMod.getReplacementResidue() + "@" + subMod.getLocation();

    }

    // Added by Fawaz Ghali 13/05/2014 exportProteoAnnotator
    private String proteoAnnotatorLineToString(ProteinAmbiguityGroup pag) {

        String line = "";
        List<UserParam> userParams = pag.getUserParam();
        String countNonA = "";
        String scoreNonA = "";
        String nonAPeptide = "";
        String aGenes = "";
        String qValue = "";
        for (int i = 0; i < userParams.size(); i++) {
            UserParam userParam = userParams.get(i);

            if (userParam.getName().equals("nonAPeptide")) {
                nonAPeptide = userParam.getValue();
            }
        }
        List<CvParam> cvParamList = pag.getCvParam();
        for (int i = 0; i < cvParamList.size(); i++) {
            CvParam cvParam = cvParamList.get(i);
            if (cvParam.getAccession().equals("MS:1002474")) {
                scoreNonA = cvParam.getValue();
            }

            if (cvParam.getAccession().equals("MS:1002475")) {
                countNonA = cvParam.getValue();
            }
            if (cvParam.getAccession().equals("MS:1002373")) {
                qValue = cvParam.getValue();
            }

        }

        List<ProteinDetectionHypothesis> proteinDetectionHypothesisList = pag.getProteinDetectionHypothesis();
        for (int i = 0; i < proteinDetectionHypothesisList.size(); i++) {
            ProteinDetectionHypothesis proteinDetectionHypothesis = proteinDetectionHypothesisList.get(i);
            if (proteinDetectionHypothesis.getDBSequenceRef().startsWith("dbseq_generic|A_")) {
                aGenes = aGenes + proteinDetectionHypothesis.getDBSequenceRef() + ";";
            }
        }

        line = countNonA + sep + scoreNonA + sep + nonAPeptide + sep + aGenes + sep + qValue;
        return line;
    }
}
