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
import java.util.Map.Entry;


import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.MzIdentMLObject;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.multiplesearch.TreeSortForIndices;
import uk.ac.liv.mzidlib.util.FileHandler;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 * @author Fawaz Ghali, University of Liverpool, 2013
 * Last updated 2016.
 */
public class FalseDiscoveryRateGlobal {

    // Determine FDR level 1) PSM, 2) Peptide, 3)ProteinGroup
    private String fdrLevel = "PSM";
    // PAG or PDH
    private String proteinLevel = "PDH";
    private double decoyRatio = 1;
    private String decoy;
    private boolean usingFileDecoyAttribute = false;
    //TODO: Make this as an array of possible values than a long string
    private String allowedEvalues = "MS:1001330" + ";" + "MS:1001172" + ";" + "MS:1001159" + ";" + "MS:1001328" + ";" + "MS:1002045" + ";" + "MS:1002053";
    private boolean betterScoresAreLower = true;
    /*
     * For each Spectrum Result id we need to store :: - peptides associated-
     * scores - evalue - decoy = true/false
     */
    private List<String> spectrumResult = new ArrayList<>();
    private List<String> spectrumItem = new ArrayList<>();
    private List<String> peptideNames = new ArrayList<>();
    private List<Double> evalues = new ArrayList<>();
    private List<String> decoyOrNot = new ArrayList<>();

    /*
     * the above information in sorted order
     */
    private List<String> sorted_spectrumResult = new ArrayList<>();
    private List<String> sorted_spectrumItem = new ArrayList<>();
    private List<String> sorted_peptideNames = new ArrayList<>();
    private List<Double> sorted_evalues = new ArrayList<>();
    private List<String> sorted_decoyOrNot = new ArrayList<>();

    /*
     * Store the estimated FDR, q-Value, FDR Score here
     */
    private List<Double> estimated_simpleFDR = new ArrayList<>();
    private List<Double> estimated_qvalue = new ArrayList<>();
    private List<Double> estimated_fdrscore = new ArrayList<>();
    private List<Double> tp = new ArrayList<>();
    private List<Double> fp = new ArrayList<>();
    /*
     * MzIdentML elements
     */
    // private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private Map<String, DBSequence> dBSequenceMap = new HashMap<>();
    private Map<String, PeptideEvidence> peptideEvidenceMap = new HashMap<>();
    private Map<String, Peptide> peptideMap = new HashMap<>();

    /*
     * The data read from MzIdentML file after parsing
     */
    private Map<String, String> peptideIdAndSequence = new HashMap<>();

    // metadata
    private AnalysisSoftwareList analysisSoftwareList;
    private AuditCollection auditCollection;
    private Provider provider;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private CvList cvList;
    private AnalysisCollection analysisCollection;
    private Inputs inputs;
    private String searchDatabase_Ref;

    // Store all PSM related information to a specific peptide
    private Map<String, List<String>> peptidePSMMap = new HashMap<>();
    private Map<String, Integer> peptidePSMCount = new HashMap<>();
    private Map<String, String> bestScorePSM = new HashMap<>();

    private MzidLibUtils mzidLibUtils;



    public static void main(String args[])
    {
        if (args.length < 8) {
            System.err.println("Missing parameters. Parameters are: \n" +
                    " InputFile,\n decoyRatio,\n decoyValue,\n cvTerm,\n Beter Scores are Lower of Not (true/false),\n fdrLevel,\n ProteinLevel,\n OutputFile");
            return;
        }
        FalseDiscoveryRateGlobal fdrGlobal = new FalseDiscoveryRateGlobal(args[0], args[1], args[2], args[3], Boolean.valueOf(args[4]), args[5], args[6]);
        fdrGlobal.computeFDRusingJonesMethod();
        fdrGlobal.writeToMzIdentMLFile(args[7]);
    }

    public FalseDiscoveryRateGlobal(String mzid, String decoyRatio, String decoy, String cvTerm, boolean betterScoresAreLower, String fdrLevel, String proteinLevel)
    {
        this.decoyRatio = Double.valueOf(decoyRatio);
        this.decoy = decoy;
        this.fdrLevel = fdrLevel;
        this.proteinLevel = proteinLevel;

        if (!cvTerm.isEmpty()) {
            allowedEvalues = cvTerm;
        }

        this.betterScoresAreLower = betterScoresAreLower;
        if (decoy == null || decoy.isEmpty()) {
            usingFileDecoyAttribute = true;
        }
        
        try {
            //File mzidFile = FileHandler.handleFile(mzid, true, true);
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid));
            readMzIdentML();
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    /**
     * Compute FDR score and q value etc using the method described in Jones et
     * al. Proteomics, 2009,9, 1220-1229
     */
    public void computeFDRusingJonesMethod()
    {
        getEvalueSortedPeptideList();
        computeSimpleFDR();
        computeQValues();
        computeFDRScore();
    }


    // Write the sorted data into a file
    public void writeToMzIdentMLFile(String fileName)
    {
        writeMzidFile(fileName);
    }


    private void readMetaData()
    {
        mzidLibUtils = new MzidLibUtils();
        cvList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
        analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
        auditCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
        provider = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);
        analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
        analysisCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
        inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
        searchDatabase_Ref = inputs.getSearchDatabase().get(0).getId();
    }

    private void readSequenceCollection()
    {
        Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        while (iterDBSequence.hasNext()) {
            DBSequence dBSequence = iterDBSequence.next();
            dBSequenceMap.put(dBSequence.getId(), dBSequence);
        }
        getPeptideEvidenceMap().clear();
        Iterator<PeptideEvidence> iterPeptideEvidence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
        while (iterPeptideEvidence.hasNext()) {
            PeptideEvidence pe = iterPeptideEvidence.next();
            peptideEvidenceMap.put(pe.getId(), pe);
        }
        Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
        while (iterPeptide.hasNext()) {
            Peptide peptide = iterPeptide.next();
            peptideMap.put(peptide.getId(), peptide);
            String pepId = peptide.getId();
            String pepSeq = peptide.getPeptideSequence();
            peptideIdAndSequence.put(pepId, pepSeq);
        }
    }

    private void readMzIdentML()
    {
        readMetaData();
        readSequenceCollection();
        switch (fdrLevel) {
            case "Peptide":
                readMzIdentMLPeptide();
                break;
            case "ProteinGroup":
                readMzIdentMLProteinGroup();
                break;
            case "PSM":
                readMzIdentMLPSM();
            default:
                System.err.println("Unknown FDR Level:" + fdrLevel);
                break;
        }
    }

    private void readMzIdentMLPeptide()
    {
        try {

            Iterator<MzIdentMLObject> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.next();
                String spectrumResultId = spectrumIdentificationResult.getSpectrumID();

                for (SpectrumIdentificationItem spectrumIdentItem : spectrumIdentificationResult.getSpectrumIdentificationItem()) {
                    List<CvParam> cvParamListSpectrumIdentificationItem = spectrumIdentItem.getCvParam();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentItem.getPeptideEvidenceRef();
                    boolean isdecoy = false;
                    for (PeptideEvidenceRef peptideEvidenceRef1 : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence1 = getPeptideEvidenceMap().get(peptideEvidenceRef1.getPeptideEvidenceRef());
                        if (peptideEvidence1 != null) {
                            isdecoy = peptideEvidence1.isIsDecoy();
                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence1.isIsDecoy()) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                if (decoy != null && !decoy.equals("") && (decoy.length() > 1)) {
                                    if (dBSequenceMap.get(peptideEvidence1.getDBSequenceRef()).getAccession().contains(decoy)) {
                                        isdecoy = true;
                                        break;
                                    }
                                } else {
                                    System.out.println("Error - no decoy value set, need to use alternative constructor");
                                }
                            }
                        }
                    }
                    double combinedFDR = 0;
                    for (CvParam cvParam : cvParamListSpectrumIdentificationItem) {
                        String accession = cvParam.getAccession();

                        if (allowedEvalues.contains(accession)) {
                            combinedFDR = Double.valueOf(cvParam.getValue());

                        }
                    }

                    String peptideRef = spectrumIdentItem.getPeptideRef();

                    if (!peptidePSMMap.containsKey(peptideRef)) {
                        List<String> tmp = new ArrayList<>();
                        tmp.add(spectrumIdentItem.getId());
                        peptidePSMMap.put(peptideRef, tmp);
                        peptidePSMCount.put(peptideRef, 1);
                    } else if (peptidePSMMap.containsKey(peptideRef) && !peptidePSMMap.get(peptideRef).contains(spectrumIdentItem.getId())) {
                        peptidePSMMap.get(peptideRef).add(spectrumIdentItem.getId());
                        int oldCount = peptidePSMCount.get(peptideRef);
                        oldCount = oldCount + 1;
                        peptidePSMCount.put(peptideRef, oldCount);
                    }

                    if (!bestScorePSM.containsKey(peptideRef)) {
                        String newKey = spectrumIdentItem.getId() + ":_:" + combinedFDR + ":_:" + Boolean.valueOf(isdecoy).toString() + ":_:" + spectrumResultId;
                        bestScorePSM.put(peptideRef, newKey);
                    } else {
                        String value = bestScorePSM.get(peptideRef);
                        if (value != null && value.contains(":_:")) {
                            String[] psm_fdr = value.split(":_:");
                            if (combinedFDR < Double.valueOf(psm_fdr[1])) {
                                bestScorePSM.remove(peptideRef);
                                String newKey = spectrumIdentItem.getId() + ":_:" + combinedFDR + ":_:" + Boolean.valueOf(isdecoy).toString() + ":_:" + spectrumResultId;
                                bestScorePSM.put(peptideRef, newKey);

                            }
                        }
                    }
                }
            }

            spectrumResult.clear();
            spectrumItem.clear();
            peptideNames.clear();
            evalues.clear();
            decoyOrNot.clear();

            for (Entry<String, String> entry : bestScorePSM.entrySet()) {
                String peptideRef = entry.getKey();
                String value = entry.getValue();
                String[] values = value.split(":_:");
                peptideNames.add(peptideRef);
                spectrumItem.add(values[0]);
                evalues.add(Double.valueOf(values[1]));
                decoyOrNot.add(values[2]);
                spectrumResult.add(values[3]);

            }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    private void readMzIdentMLProteinGroup()
    {
        try {

            Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller
                    .unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {
                ProteinAmbiguityGroup proteinAmbiguityGroup = iterProteinAmbiguityGroup.next();
                boolean anchorProtein = false;
                ProteinDetectionHypothesis anchorProteinDetectionHypothesis = null;

                String proteinAmbiguityGroupId = proteinAmbiguityGroup.getId();
                double combinedFDR = 0;
                boolean isdecoy = false;
                //
                String dbSeqRef;

                for (int j = 0; j < proteinAmbiguityGroup.getProteinDetectionHypothesis().size(); j++) {
                    ProteinDetectionHypothesis proteinDetectionHypothesis = proteinAmbiguityGroup
                            .getProteinDetectionHypothesis().get(j);
                    List<CvParam> cvParamListproteinDetectionHypothesis = proteinDetectionHypothesis.getCvParam();
                    for (CvParam cvParam : cvParamListproteinDetectionHypothesis) {
                        String accession = cvParam.getAccession();
                        if (accession.equals("MS:1001591")) {
                            anchorProteinDetectionHypothesis = proteinDetectionHypothesis;
                            anchorProtein = true;
                            break;
                        }
                    }
                }
                if (!anchorProtein) {
                    anchorProteinDetectionHypothesis = proteinAmbiguityGroup.getProteinDetectionHypothesis().get(0);

                }
                List<PeptideHypothesis> peptideHypothesisList = anchorProteinDetectionHypothesis.getPeptideHypothesis();
                dbSeqRef = anchorProteinDetectionHypothesis.getDBSequenceRef();

                for (PeptideHypothesis peptideHypothesis : peptideHypothesisList) {
                    PeptideEvidence peptideEvidence1 = getPeptideEvidenceMap()
                            .get(peptideHypothesis.getPeptideEvidenceRef());
                    if (peptideEvidence1 != null) {
                        isdecoy = peptideEvidence1.isIsDecoy();
                        if (usingFileDecoyAttribute) {
                            if (peptideEvidence1.isIsDecoy()) {
                                isdecoy = true;
                                break;
                            }
                        } else {
                            if (decoy != null && !decoy.equals("") && (decoy.length() > 1)) {
                                if (dBSequenceMap.get(peptideEvidence1.getDBSequenceRef()).getAccession()
                                        .contains(decoy)) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                System.out.println("Error - no decoy value set, need to use alternative constructor");
                            }
                        }
                    }
                }


                if (proteinLevel.equals("PAG")) {
                    for (CvParam cvParam : proteinAmbiguityGroup.getCvParam()) {
                        String accession = cvParam.getAccession();
                        if (allowedEvalues.contains(accession)) {
                            combinedFDR = Double.valueOf(cvParam.getValue());
                        }
                    }
                }
                if (proteinLevel.equals("PDH")) {
                    for (CvParam cvParam : anchorProteinDetectionHypothesis.getCvParam()) {
                        String accession = cvParam.getAccession();
                        if (allowedEvalues.contains(accession)) {
                            combinedFDR = Double.valueOf(cvParam.getValue());
                        }
                    }

                }

                spectrumResult.add(proteinAmbiguityGroupId);
                spectrumItem.add(anchorProteinDetectionHypothesis.getId());
                peptideNames.add(dbSeqRef);
                evalues.add(combinedFDR);
                decoyOrNot.add(String.valueOf(isdecoy));

            }

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    private void readMzIdentMLPSM()
    {
        try {
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult = iterSpectrumIdentificationResult.next();
                String spectrumResultId = spectrumIdentificationResult.getSpectrumID();
                for (SpectrumIdentificationItem spectrumIdentItem : spectrumIdentificationResult.getSpectrumIdentificationItem()) {
                    String peptideRef = spectrumIdentItem.getPeptideRef();
                    List<CvParam> cvParamListSpectrumIdentificationItem = spectrumIdentItem.getCvParam();
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentItem.getPeptideEvidenceRef();
                    boolean isdecoy = false;
                    for (PeptideEvidenceRef peptideEvidenceRef1 : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence1 = getPeptideEvidenceMap().get(peptideEvidenceRef1.getPeptideEvidenceRef());
                        if (peptideEvidence1 != null) {
                            isdecoy = peptideEvidence1.isIsDecoy();
                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence1.isIsDecoy()) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                if (decoy != null && !decoy.equals("") && (decoy.length() > 1)) {
                                    if (dBSequenceMap.get(peptideEvidence1.getDBSequenceRef()).getAccession().contains(decoy)) {
                                        isdecoy = true;
                                        break;
                                    }
                                } else {
                                    System.out.println("Error - no decoy value set, need to use alternative constructor");
                                }
                            }
                        }
                    } // end of for each

                    if (isdecoy) {
                        decoyOrNot.add("true");
                    } else {
                        decoyOrNot.add("false");

                    }
                    spectrumResult.add(spectrumResultId);
                    peptideNames.add(peptideRef);
                    spectrumItem.add(spectrumIdentItem.getId());
                    for (CvParam cvParam : cvParamListSpectrumIdentificationItem) {
                        String accession = cvParam.getAccession();
                        if (allowedEvalues.contains(accession)) {
                            evalues.add(new Double(cvParam.getValue()));
                        }

                    }
                }
            }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    /**
     * Sort the data using evalue
     * Make simple arrayLists from complicated data structures returned after
     * reading mzIdentML file. These simple arrayLists will be used for further
     * calculations of various scores.
     */
    private void getEvalueSortedPeptideList()
    {
        try {
            // Call the sorting routine to find the indices of sorted evalues
            TreeSortForIndices sortClass = new TreeSortForIndices();
            Integer[] sortOrderForEvalues = sortClass.sortTheValueColumn(evalues.toArray(), betterScoresAreLower);

            // Arrange the values in the order determined by the sorting
            // operation, such that, each index an arraylist can be map to
            // the entries in other arraylists for the "same" index
            int index;
            for (Integer sortOrderForEvalue : sortOrderForEvalues) {
                index = sortOrderForEvalue;
                sorted_spectrumResult.add(spectrumResult.get(index));
                sorted_spectrumItem.add(spectrumItem.get(index));
                sorted_peptideNames.add(peptideNames.get(index));
                sorted_evalues.add(evalues.get(index));
                sorted_decoyOrNot.add(decoyOrNot.get(index));

            }

            // Clear the memory for the items no more needed
            spectrumResult.clear();
            peptideNames.clear();
            spectrumItem.clear();
            evalues.clear();
            decoyOrNot.clear();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }


    private void writeUnchangedPart(Writer writer, MzIdentMLMarshaller marshaller)
    {
        try {

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
                if (fdrLevel.equals("Peptide")) {
                    String peptideRef = pe.getPeptideRef();
                    int psmCount = 0;
                    if (peptidePSMCount.get(peptideRef) != null) {
                        psmCount = peptidePSMCount.get(peptideRef);
                    }
                    pe.getUserParam().add(makeUserParam("psm_count", String.valueOf(psmCount)));
                }
                sequenceCollection.getPeptideEvidence().add(pe);
            }

            marshaller.marshal(sequenceCollection, writer);
            writer.write("\n");

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
                Cv cv = mzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);
                CvParam peptideLevelScoring = new CvParam();
                peptideLevelScoring.setAccession("MS:1002490");
                peptideLevelScoring.setName("peptide-level scoring");
                peptideLevelScoring.setCv(cv);
                analysisProtocolCollection.getSpectrumIdentificationProtocol().get(0).getAdditionalSearchParams()
                        .getCvParam().add(peptideLevelScoring);
                CvParam groupPSMsBySequence = new CvParam();
                groupPSMsBySequence.setAccession("MS:1002496");
                groupPSMsBySequence.setName("group PSMs by sequence");
                groupPSMsBySequence.setCv(cv);
                analysisProtocolCollection.getSpectrumIdentificationProtocol().get(0).getAdditionalSearchParams()
                        .getCvParam().add(groupPSMsBySequence);
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

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    //TODO: This method has to be broken into PSM, Peptide or ProteinGroup-specific functions
    //Not only to increase the readability, but to avoid if-conditions inside loops. Grrr.
    //Can't justify a function with 200 lines of code!
    private void writeMzidFile(String csvFileName)
    {
        try {
            Writer writer = new FileWriter(csvFileName);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            writeUnchangedPart(writer, marshaller);

            SpectrumIdentificationList siList = new SpectrumIdentificationList();
            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            siList.setId(spectrumIdentificationListRef);

            Iterator<SpectrumIdentificationList> iterSpectrumIdentificationList = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationList);
            while (iterSpectrumIdentificationList.hasNext()) {
                SpectrumIdentificationList sr = iterSpectrumIdentificationList.next();
                siList.getCvParam().addAll(sr.getCvParam());
            }

            Iterator<FragmentationTable> iterFragmentationTable = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {
                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Map<String, String> sorted_stringMap = new HashMap<>();
            switch (fdrLevel) {
                case "PSM":
                    for (int i = 0; i < sorted_spectrumItem.size(); i++) {
                        sorted_stringMap.put(sorted_spectrumItem.get(i), estimated_simpleFDR.get(i).toString() + ":"
                                + String.valueOf(estimated_qvalue.get(i)) + ":" + estimated_fdrscore.get(i).toString());
                    }
                    break;
                case "Peptide":
                    for (int i = 0; i < sorted_peptideNames.size(); i++) {
                        sorted_stringMap.put(sorted_peptideNames.get(i), estimated_simpleFDR.get(i).toString() + ":"
                                + String.valueOf(estimated_qvalue.get(i)) + ":" + estimated_fdrscore.get(i).toString());
                    }
                    break;
                default:
                    for (int i = 0; i < sorted_spectrumItem.size(); i++) {
                        sorted_stringMap.put(sorted_spectrumResult.get(i),
                                estimated_simpleFDR.get(i).toString() + ":" + String.valueOf(estimated_qvalue.get(i)));
                    }
                    break;
            }

            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult sr = iterSpectrumIdentificationResult.next();
                for (SpectrumIdentificationItem sii : sr.getSpectrumIdentificationItem()) {
                    String sortedValue = null;
                    if (fdrLevel.equals("PSM")) {
                        sortedValue = sorted_stringMap.get(sii.getId());
                    } else if (fdrLevel.equals("Peptide")) {
                        sortedValue = sorted_stringMap.get(sii.getPeptideRef());
                    }

                    if (sortedValue != null && fdrLevel.equals("PSM")) {
                        CvParam cvParamestimated_simpleFDR = new CvParam();
                        CvParam cvParamestimated_qvalue = new CvParam();
                        CvParam cvParamfdrscore = new CvParam();
                        String[] sii_arr = sortedValue.split(":");
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

                    if (sortedValue != null && fdrLevel.equals("Peptide")) {
                         //Check if CVs are already in the file
                        List<CvParam> listParams = sii.getCvParam();
                        Iterator<CvParam> i = listParams.iterator();
                        while (i.hasNext()) {
                            CvParam get = i.next();
                            if (get.getAccession().equals("MS:1002359")
                                    || get.getAccession().equals("MS:1001868")
                                    || get.getAccession().equals("MS:1002360")
                                    || get.getAccession().equals("MS:1002500")
                                    || get.getAccession().equals("MS:1002520")) {
                                i.remove();
                            }
                        }
                        CvParam cvParamestimated_simpleFDR = new CvParam();
                        CvParam cvParamestimated_qvalue = new CvParam();
                        CvParam cvParamfdrscore = new CvParam();
                        String[] sii_arr = sortedValue.split(":");
                        Cv cv = mzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);

                        cvParamestimated_simpleFDR.setAccession("MS:1002359");
                        cvParamestimated_simpleFDR.setName("distinct peptide-level local FDR");
                        cvParamestimated_simpleFDR.setValue(sii_arr[0]);
                        cvParamestimated_simpleFDR.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_simpleFDR);

                        cvParamestimated_qvalue.setAccession("MS:1001868");
                        cvParamestimated_qvalue.setName("distinct peptide-level q-value");
                        cvParamestimated_qvalue.setValue(sii_arr[1]);
                        cvParamestimated_qvalue.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_qvalue);

                        cvParamfdrscore.setAccession("MS:1002360");
                        cvParamfdrscore.setName("distinct peptide-level FDRScore");
                        cvParamfdrscore.setValue(sii_arr[2]);
                        cvParamfdrscore.setCv(cv);
                        sii.getCvParam().add(cvParamfdrscore);

                        // Added by Fawaz Ghali: update to mzid 1.2 20/01/2015
                        CvParam peptidePassesThreshold = new CvParam();
                        peptidePassesThreshold.setAccession("MS:1002500");
                        peptidePassesThreshold.setName("peptide passes threshold");
                        peptidePassesThreshold.setValue("true");
                        peptidePassesThreshold.setCv(cv);
                        sii.getCvParam().add(peptidePassesThreshold);
                        CvParam peptideGroupID = new CvParam();
                        peptideGroupID.setAccession("MS:1002520");
                        peptideGroupID.setName("peptide group ID");
                        peptideGroupID.setValue(sii.getPeptideRef());
                        peptideGroupID.setCv(cv);
                        sii.getCvParam().add(peptideGroupID);

                    }

                }
                siList.getSpectrumIdentificationResult().add(sr);
            }

            marshaller.marshal(siList, writer);
            writer.write("\n");

            if (fdrLevel.equals("PSM") || fdrLevel.equals("Peptide")) {
                writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
                writer.write(marshaller.createDataCollectionClosingTag() + "\n");

                writer.write(marshaller.createMzIdentMLClosingTag());
            } else {
                ProteinDetectionList pdList = new ProteinDetectionList();
                pdList.setId("PDL_1");
                Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller
                        .unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterProteinAmbiguityGroup.hasNext()) {
                    ProteinAmbiguityGroup proteinAmbiguityGroup = iterProteinAmbiguityGroup.next();

                    for (ProteinDetectionHypothesis proteinDetectionHypothesis : proteinAmbiguityGroup
                            .getProteinDetectionHypothesis()) {
                        List<CvParam> cvParamListproteinDetectionHypothesis = proteinDetectionHypothesis.getCvParam();
                        for (CvParam cvParam : cvParamListproteinDetectionHypothesis) {
                            String accession = cvParam.getAccession();

                            if (accession.equals("MS:1001591")) {
                                String sortedValue = sorted_stringMap.get(proteinAmbiguityGroup.getId());
                                CvParam cvParamestimated_simpleFDR = new CvParam();
                                CvParam cvParamestimated_qvalue = new CvParam();
                                String[] sii_arr = sortedValue.split(":");
                                Cv cv = mzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);
                                cvParamestimated_simpleFDR.setAccession("MS:1002370");
                                cvParamestimated_simpleFDR.setName("protein group-level local FDR");
                                cvParamestimated_simpleFDR.setValue(sii_arr[0]);
                                cvParamestimated_simpleFDR.setCv(cv);

                                proteinAmbiguityGroup.getCvParam().add(cvParamestimated_simpleFDR);

                                cvParamestimated_qvalue.setAccession("MS:1002373");
                                cvParamestimated_qvalue.setName("protein group-level q-value");
                                cvParamestimated_qvalue.setValue(sii_arr[1]);
                                cvParamestimated_qvalue.setCv(cv);

                                proteinAmbiguityGroup.getCvParam().add(cvParamestimated_qvalue);
                            }
                        }
                    }

                    pdList.getProteinAmbiguityGroup().add(proteinAmbiguityGroup);
                }
                marshaller.marshal(pdList, writer);
                writer.write("\n");
                writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
                writer.write(marshaller.createDataCollectionClosingTag() + "\n");

                writer.write(marshaller.createMzIdentMLClosingTag());
                writer.write("\n");

            }
            writer.close();

            System.out.println("Output written to " + csvFileName);

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className
                    + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }





    /**
     * Compute simple FDR
     */
    private void computeSimpleFDR()
    {
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

            falsePositiveDivRatio = (double) falsePositiveCount / decoyRatio;

            simpleFDR = falsePositiveDivRatio / (double) allTargets;
            estimated_simpleFDR.add(simpleFDR);
            estimated_qvalue.add(0d);

            fp.add(falsePositiveDivRatio);
            double tpValue = allTargets - falsePositiveDivRatio;
            tp.add(tpValue);
        }
        if (fp.size() == 0 || (fp.get(fp.size() - 1) == 0)) {
            System.out.println("No decoys found for search engine Mascot|Omssa|tandem - likely caused by: wrong decoy regex, database doesn't contain decoys or the search reported only identifications for stringent e-values - please allow identifications up to e-value = 10 for omssa and tandem. See omssa and tandem documentation for how to do change this setting");
        }
    }

    /**
     * Compute q-value
     */
    private void computeQValues()
    {
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
     */
    private void computeFDRScore()
    {
        // This method cannot be used for ProteinGroup FDRLevels
        if (fdrLevel.equals("ProteinGroup")) {
            return;
        }
        // Initialize the estimated_fdrscore by adding elements containing 0.
        // This needs to be
        // done because of back tracking involved in the algorithm, so simple
        // .add() wouldn't work
        for (int i = 0; i < sorted_peptideNames.size(); i++) {
            estimated_fdrscore.add(0d);

        }

        // previous evalue in case of straight vertical rise in q value without
        // any change in evalue
        double prev_evalue = 0f;
        double prev_qvalue = 0f;
        double prev_prev_evalue = 0f;

        int counter_backwardStep = 0;

        int lastNonZeroScoreIdx = -1;
        int i=0;
        for (; i < sorted_peptideNames.size(); i++) {
            double current_evalue = sorted_evalues.get(i);
            double current_qvalue = estimated_qvalue.get(i);

            if (current_qvalue > prev_qvalue) {

                double slope;
                double intercept;

                // Work out the slope and the intercept
                // Tells us which co-ordinates to use for
                if (current_evalue != prev_evalue) {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue - prev_evalue);
                    intercept = prev_qvalue - slope * prev_evalue;
                } else {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue - prev_prev_evalue);
                    intercept = prev_qvalue - slope * prev_prev_evalue;
                }

                if (counter_backwardStep > 0) { // compute the FDR score for flat q-value region
                    for (int k = 0; k <= counter_backwardStep; k++) {
                        int index = i - counter_backwardStep + k;
                        double fdrScore = slope * sorted_evalues.get(index) + intercept;
                        estimated_fdrscore.set(index, fdrScore);
                        lastNonZeroScoreIdx = fdrScore>0 && index > lastNonZeroScoreIdx?index:lastNonZeroScoreIdx;
                    }
                } else { // In case an immediate increment in q value is found
                    double fdrScore = slope * current_evalue + intercept;
                    estimated_fdrscore.set(i, fdrScore);
                    lastNonZeroScoreIdx = fdrScore>0 && i > lastNonZeroScoreIdx?i:lastNonZeroScoreIdx;
                }

                // Re-initialise the variables
                counter_backwardStep = 0;
                if (current_evalue > prev_evalue) { // the previous e-value will
                    // change only if the current e-value is different
                    prev_prev_evalue = prev_evalue;
                    prev_evalue = current_evalue;
                }
                prev_qvalue = current_qvalue;

            } else {
                counter_backwardStep++;
            }
        }

        // In case if we miss to update the very last values in
        // estimated_fdrscore because
        // current_qvalue == prev_qvalue
        //TODO: Pretty dangerous tactic here in terms of manipulating an expired loop trip count.
        //TODO: Fix this please! Seriously who on the earth would write a code like this? Must be sacked!
        if (estimated_fdrscore.get(i - 1) == 0) {
            double lastFdrValue;
            i = i - 1;
            // TODO: Appears totally novel tactic to locate a last occuring, non-zero estimated score.
            // TODO: Record the index while putting values! IOOOOO!
            while (estimated_fdrscore.get(i) == 0 && i > 0) {
                i--;
            }
            if (i == 0) {
                throw new RuntimeException(
                        "\n Can't compute FDR. Likely that all the hits are corrects, so nothing to do.");
            } else {
                System.out.println("Last FDR Value Found at " + i + " and pre-recorded version is: " + lastNonZeroScoreIdx);
                lastFdrValue = estimated_fdrscore.get(i);
                while (i < estimated_fdrscore.size()) {
                    estimated_fdrscore.set(i, lastFdrValue);
                    i++;
                }
            }
        }


//        double lastFdrValue = estimated_fdrscore.get(lastNonZeroScoreIdx);
//        for (int i=lastNonZeroScoreIdx; i < estimated_fdrscore.size(); i++) {
//            estimated_fdrscore.set(i, lastFdrValue);
//        }

    }

    @SuppressWarnings("unused")
    public List<String> getSorted_spectrumResult()
    {
        return sorted_spectrumResult;
    }

    @SuppressWarnings("unused")
    public List<String> getSorted_peptideNames()
    {
        return sorted_peptideNames;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_evalues()
    {
        return sorted_evalues;
    }

    @SuppressWarnings("unused")
    public List<String> getSorted_decoyOrNot()
    {
        return sorted_decoyOrNot;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_simpleFDR()
    {
        return estimated_simpleFDR;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_qValues()
    {
        return estimated_qvalue;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_estimatedFDR()
    {
        return estimated_fdrscore;
    }

    @SuppressWarnings("unused")
    public List<Double> getTP()
    {
        return tp;
    }

    @SuppressWarnings("unused")
    public List<Double> getFP()
    {
        return fp;
    }

    /**
     * Clear all the data structures in the class to release all the memory
     */
    @SuppressWarnings("unused")
    public void clearAllData()
    {
        sorted_spectrumResult.clear();
        sorted_peptideNames.clear();
        sorted_evalues.clear();
        sorted_decoyOrNot.clear();
        estimated_simpleFDR.clear();
        estimated_qvalue.clear();
        estimated_fdrscore.clear();
        fp.clear();
        tp.clear();
    }

    @SuppressWarnings("unused")
    public Map<String, DBSequence> getdBSequenceHashMap()
    {
        return dBSequenceMap;
    }

    @SuppressWarnings("unused")
    public void setdBSequenceHashMap(Map<String, DBSequence> dBSequenceMap)
    {
        this.dBSequenceMap = dBSequenceMap;
    }

    @Deprecated
    @SuppressWarnings("unused")
    public Map<String, PeptideEvidence> getPeptideEvidenceHashMap()
    {
        return peptideEvidenceMap;
    }

    @SuppressWarnings("unused")
    public void setPeptideEvidenceHashMap(Map<String, PeptideEvidence> peptideEvidenceMap)
    {
        this.peptideEvidenceMap = peptideEvidenceMap;
    }

    @SuppressWarnings("unused")
    private Map<String, PeptideEvidence> getPeptideEvidenceMap()
    {
        return peptideEvidenceMap;
    }

    @SuppressWarnings("unused")
    public Map<String, Peptide> getPeptideHashMap()
    {
        return peptideMap;
    }

    @SuppressWarnings("unused")
    public void setPeptideHashMap(Map<String, Peptide> peptideMap)
    {
        this.peptideMap = peptideMap;
    }

    public AnalysisSoftwareList getAnalysisSoftwareList()
    {
        return analysisSoftwareList;
    }

    public AuditCollection getAuditCollection()
    {
        return auditCollection;
    }

    public Provider getProvider()
    {
        return provider;
    }

    public AnalysisProtocolCollection getAnalysisProtocolCollection()
    {
        return analysisProtocolCollection;
    }

    public CvList getCvList()
    {
        return cvList;
    }

    public AnalysisCollection getAnalysisCollection()
    {
        return analysisCollection;
    }

    public Inputs getInputs()
    {
        return inputs;
    }

    @SuppressWarnings("unused")
    public Map<String, String> getFromXMLPeptideSequenceHash()
    {
        return peptideIdAndSequence;
    }

    @SuppressWarnings("unused")
    public String getSearchDatabase_Ref()
    {
        return searchDatabase_Ref;
    }


    @SuppressWarnings("unused")
    // Write the sorted data into a file
    public void writeToTsvFile(String fileName) throws Exception
    {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        String outStr;
        for (int i = 0; i < sorted_peptideNames.size(); i++) {

            outStr = sorted_spectrumResult.get(i) + "\t" + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i)
                    + "\t"
                    + estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.get(i) + "\t" + estimated_fdrscore.get(i)
                    + "\n";

            out.write(outStr);
        }
        out.close();
    }

    // Write the sorted data into a file
    @SuppressWarnings("unused")
    public void writeToCsvFile(String fileName) throws Exception
    {

        Writer out = new BufferedWriter(new FileWriter(fileName));
        String outStrHead = "sorted_spectrumResult.get(i),sorted_peptideNames.get(i) , sorted_decoyOrNot.get(i) ,  sorted_evalues.get(i).toString() , + sorted_scores.get(i).toString() , estimated_simpleFDR.get(i) , estimated_qvalue.get(i) , estimated_fdrscore.get(i) \n";

        for (int i = 0; i < sorted_peptideNames.size(); i++) {

            String outStr = sorted_spectrumResult.get(i) + "," + sorted_peptideNames.get(i) + ","
                    + sorted_decoyOrNot.get(i) + "," + sorted_evalues.get(i).toString() + ","
                    + estimated_simpleFDR.get(i) + "," + estimated_qvalue.get(i) + "," + estimated_fdrscore.get(i)
                    + "\n";

            out.write(outStr);
        }
        out.close();
    }

    //  Write the sorted data into a file
    @SuppressWarnings("unused")
    public void writeTheSortedDataToFile(String fileName) throws Exception
    {

        // Basic check to see that each peptide has a evalue
        if (sorted_peptideNames.size() != sorted_evalues.size()) {
            throw (new Exception("Number of entries = " + sorted_peptideNames.size()
                    + "in sorted_peptideNames don't match with the number of entries = " + sorted_evalues.size()
                    + "in sorted_evalues"));
        }

        Writer out = new BufferedWriter(new FileWriter(fileName));

        for (int i = 0; i < sorted_peptideNames.size(); i++) {
            String outStr = sorted_spectrumResult.get(i) + "\t" + sorted_spectrumItem.get(i) + "\t"
                    + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i) + "\t"
                    + sorted_evalues.get(i).toString() + "\n";

            out.write(outStr);
        }

        out.close();
    }


    private UserParam makeUserParam(String name, String value)
    {
        UserParam userParam = new UserParam();
        userParam.setName(name);
        userParam.setValue(value);
        return userParam;
    }
}
