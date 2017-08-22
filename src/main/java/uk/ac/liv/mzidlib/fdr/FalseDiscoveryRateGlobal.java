
package uk.ac.liv.mzidlib.fdr;

import static uk.ac.liv.mzidlib.constants.CvConstants.*;

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
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.mutable.MutableInt;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentationTable;
import uk.ac.ebi.jmzidml.model.mzidml.Param;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.multiplesearch.TreeSortForIndices;
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

    //TODO: Make this as an array of possible values than a long string
    private String allowedEvalues = XTANDEM_EXPECT.getAccession() + ";"
            + MASCOT_EXPECTATION_VALUE.getAccession() + ";"
            + SEQUEST_EXPECTATION_VALUE.getAccession() + ";"
            + OMSSA_EVALUE.getAccession() + ";"
            + PROTEINPROSPECTOR_EXPECTATION_VALUE.getAccession() + ";"
            + MSGF_EVALUE.getAccession();
    private boolean betterScoresAreLower = true;

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

    private MzIdentMLVersion version;

    private FalseDiscoveryRateGlobalReader fdrgReader;

    public FalseDiscoveryRateGlobal(String mzid, String decoyRatio, String decoy,
                                    String cvTerm, boolean betterScoresAreLower,
                                    String fdrLevel, String proteinLevel,
                                    String mzidVer) {
        this.decoyRatio = Double.valueOf(decoyRatio);
        this.decoy = decoy;
        this.fdrLevel = fdrLevel;
        this.proteinLevel = proteinLevel;

        if (!cvTerm.isEmpty()) {
            allowedEvalues = cvTerm;
        }

        this.betterScoresAreLower = betterScoresAreLower;

        if (mzidVer.equals("1.1")) {
            this.version = MzIdentMLVersion.Version_1_1;
        } else if (mzidVer.equals("1.2")) {
            this.version = MzIdentMLVersion.Version_1_2;
        } else {
            System.out.println(
                    "The input mzIdentML version is not recognizable. Using the default version 1.2 for output files.");
            this.version = MzIdentMLVersion.Version_1_2;
        }

        fdrgReader = new FalseDiscoveryRateGlobalReader(new File(mzid),
                                                        this.decoy,
                                                        this.fdrLevel,
                                                        this.proteinLevel);
    }

    /**
     * Compute FDR score and q value etc using the method described in Jones et
     * al. Proteomics, 2009,9, 1220-1229
     */
    public void computeFDRusingJonesMethod() {
        getEvalueSortedPeptideList();
        computeSimpleFDR();
        computeQValues();
        computeFDRScore();
    }

    /**
     * Sort the data using evalue
     * Make simple arrayLists from complicated data structures returned after
     * reading mzIdentML file. These simple arrayLists will be used for further
     * calculations of various scores.
     */
    private void getEvalueSortedPeptideList() {
        try {
            // Call the sorting routine to find the indices of sorted evalues
            TreeSortForIndices sortClass = new TreeSortForIndices();
            Integer[] sortOrderForEvalues = sortClass.sortTheValueColumn(
                    fdrgReader.getEvalues().toArray(), betterScoresAreLower);

            // Arrange the values in the order determined by the sorting
            // operation, such that, each index an arraylist can be map to
            // the entries in other arraylists for the "same" index
            int index;
            for (Integer sortOrderForEvalue : sortOrderForEvalues) {
                index = sortOrderForEvalue;
                sorted_spectrumResult.add(fdrgReader.getSpectrumResult().get(
                        index));
                sorted_spectrumItem.add(fdrgReader.getSpectrumItem().get(index));
                sorted_peptideNames.add(fdrgReader.getPeptideNames().get(index));
                sorted_evalues.add(fdrgReader.getEvalues().get(index));
                sorted_decoyOrNot.add(fdrgReader.getDecoyOrNot().get(index));

            }

            // Clear the memory for the items no more needed
            fdrgReader.getSpectrumResult().clear();
            fdrgReader.getPeptideNames().clear();
            fdrgReader.getSpectrumItem().clear();
            fdrgReader.getEvalues().clear();
            fdrgReader.getDecoyOrNot().clear();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Compute simple FDR
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

            falsePositiveDivRatio = (double) falsePositiveCount / decoyRatio;

            simpleFDR = falsePositiveDivRatio / (double) allTargets;
            estimated_simpleFDR.add(simpleFDR);
            estimated_qvalue.add(0d);

            fp.add(falsePositiveDivRatio);
            double tpValue = allTargets - falsePositiveDivRatio;
            tp.add(tpValue);
        }
        if (fp.isEmpty() || (fp.get(fp.size() - 1) == 0)) {
            System.out.println(
                    "No decoys found for search engine Mascot|Omssa|tandem - likely caused by: wrong decoy regex, database doesn't contain decoys or the search reported only identifications for stringent e-values - please allow identifications up to e-value = 10 for omssa and tandem. See omssa and tandem documentation for how to do change this setting");
        }
    }

    /**
     * Compute q-value
     */
    private void computeQValues() {
        double immediateMinFdr = estimated_simpleFDR.get(estimated_simpleFDR.
                size() - 1);
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
    private void computeFDRScore() {
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
        int i = 0;
        for (; i < sorted_peptideNames.size(); i++) {
            double current_evalue = sorted_evalues.get(i);
            double current_qvalue = estimated_qvalue.get(i);

            if (current_qvalue > prev_qvalue) {

                double slope;
                double intercept;

                // Work out the slope and the intercept
                // Tells us which co-ordinates to use for
                if (current_evalue != prev_evalue) {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue
                            - prev_evalue);
                    intercept = prev_qvalue - slope * prev_evalue;
                } else {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue
                            - prev_prev_evalue);
                    intercept = prev_qvalue - slope * prev_prev_evalue;
                }

                if (counter_backwardStep > 0) { // compute the FDR score for flat q-value region
                    for (int k = 0; k <= counter_backwardStep; k++) {
                        int index = i - counter_backwardStep + k;
                        double fdrScore = slope * sorted_evalues.get(index)
                                + intercept;
                        estimated_fdrscore.set(index, fdrScore);
                        lastNonZeroScoreIdx = fdrScore > 0 && index
                                > lastNonZeroScoreIdx ? index
                                        : lastNonZeroScoreIdx;
                    }
                } else { // In case an immediate increment in q value is found
                    double fdrScore = slope * current_evalue + intercept;
                    estimated_fdrscore.set(i, fdrScore);
                    lastNonZeroScoreIdx = fdrScore > 0 && i
                            > lastNonZeroScoreIdx ? i : lastNonZeroScoreIdx;
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
                System.out.println("Last FDR Value Found at " + i
                        + " and pre-recorded version is: " + lastNonZeroScoreIdx);
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

    // Write the sorted data into a file
    public void writeToMzIdentMLFile(String fileName) {
        writeMzidFile(fileName);
    }

    //TODO: This method has to be broken into PSM, Peptide or ProteinGroup-specific functions
    //Not only to increase the readability, but to avoid if-conditions inside loops. Grrr.
    //Can't justify a function with 200 lines of code!
    private void writeMzidFile(String csvFileName) {

        try (Writer writer = new FileWriter(csvFileName)) {
            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller(version);

            writer.write(marshaller.createXmlHeader() + "\n");

            MzIdentMLUnmarshaller mzIdentMLUnmarshaller = fdrgReader.
                    getMzIdentMLUnmarshaller();

            String mzID = mzIdentMLUnmarshaller.getMzIdentMLId();
            if (mzID != null) {
                writer.write(marshaller.createMzIdentMLStartTag(mzID) + "\n");
            } else {
                writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            }

            if (fdrgReader.getCvList() != null) {
                marshaller.marshal(fdrgReader.getCvList(), writer);
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
            Cv psiCV = MzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);
            param.setParam(MzidLibUtils.makeCvParam("MS:1002237", "mzidLib",
                                                    psiCV));
            analysisSoftware.setSoftwareName(param);
            fdrgReader.getAnalysisSoftwareList().getAnalysisSoftware().add(
                    analysisSoftware);
            marshaller.marshal(fdrgReader.getAnalysisSoftwareList(), writer);
            writer.write("\n");

            if (fdrgReader.getProvider() != null) {
                marshaller.marshal(fdrgReader.getProvider(), writer);
            }
            writer.write("\n");

            if (fdrgReader.getAuditCollection() != null) {
                marshaller.marshal(fdrgReader.getAuditCollection(), writer);
            }
            writer.write("\n");
            SequenceCollection sequenceCollection = new SequenceCollection();
            Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.
                    unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {

                sequenceCollection.getDBSequence().add(iterDBSequence.next());

            }

            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.
                    unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide pe = iterPeptide.next();
                sequenceCollection.getPeptide().add(pe);
            }
            Iterator<PeptideEvidence> iterPeptideEvidence
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence pe = iterPeptideEvidence.next();
                if (fdrLevel.equals("Peptide")) {
                    String peptideRef = pe.getPeptideRef();
                    int psmCount = 0;
                    if (fdrgReader.getPeptidePSMCount().get(peptideRef) != null) {
                        psmCount = fdrgReader.getPeptidePSMCount().get(
                                peptideRef);
                    }
                    pe.getUserParam().add(makeUserParam("psm_count", String.
                                                        valueOf(psmCount)));
                }
                sequenceCollection.getPeptideEvidence().add(pe);
            }

            marshaller.marshal(sequenceCollection, writer);
            writer.write("\n");

            if (fdrgReader.getAnalysisCollection() != null) {
                marshaller.marshal(fdrgReader.getAnalysisCollection(), writer);
            }
            writer.write("\n");

            if (fdrgReader.getAnalysisProtocolCollection() != null) {
                List<SpectrumIdentificationProtocol> sipList
                        = fdrgReader.getAnalysisProtocolCollection().
                        getSpectrumIdentificationProtocol();
                for (SpectrumIdentificationProtocol sip : sipList) {
                    if (sip.getAdditionalSearchParams().getCvParam().contains(
                            CvConstants.NO_SPECIAL_PROCESSING)) {
                        sip.getAdditionalSearchParams().getCvParam().remove(
                                CvConstants.NO_SPECIAL_PROCESSING);
                    }

                    if (!sip.getAdditionalSearchParams().getCvParam().contains(
                            CvConstants.PEPTIDE_LEVEL_SCORING)) {
                        sip.getAdditionalSearchParams().getCvParam().add(
                                CvConstants.PEPTIDE_LEVEL_SCORING);
                    }

                    if (!sip.getAdditionalSearchParams().getCvParam().contains(
                            CvConstants.GROUP_PSMS_BY_SEQUENCE)) {
                        sip.getAdditionalSearchParams().getCvParam().add(
                                CvConstants.GROUP_PSMS_BY_SEQUENCE);
                    }
                }
                marshaller.marshal(fdrgReader.getAnalysisProtocolCollection(),
                                   writer);

            }
            writer.write("\n");
            writer.write(marshaller.createDataCollectionStartTag() + "\n");
            writer.write("\n");
            if (fdrgReader.getInputs() != null) {
                marshaller.marshal(fdrgReader.getInputs(), writer);
            }
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");

            SpectrumIdentificationList siList = new SpectrumIdentificationList();
            String spectrumIdentificationListRef = "";
            if (fdrgReader.getAnalysisCollection().getSpectrumIdentification().
                    size() > 0) {
                spectrumIdentificationListRef = fdrgReader.
                        getAnalysisCollection().
                        getSpectrumIdentification().get(0).
                        getSpectrumIdentificationListRef();
            }
            siList.setId(spectrumIdentificationListRef);

            Iterator<SpectrumIdentificationList> iterSpectrumIdentificationList
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationList);
            while (iterSpectrumIdentificationList.hasNext()) {
                SpectrumIdentificationList sr = iterSpectrumIdentificationList.
                        next();
                siList.getCvParam().addAll(sr.getCvParam());
            }

            Iterator<FragmentationTable> iterFragmentationTable
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {
                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Map<String, String> sorted_stringMap = new HashMap<>();
            switch (fdrLevel) {
                case "PSM":
                    for (int i = 0; i < sorted_spectrumItem.size(); i++) {
                        sorted_stringMap.put(sorted_spectrumItem.get(i),
                                             estimated_simpleFDR.get(i).
                                             toString() + ":"
                                             + String.valueOf(estimated_qvalue.
                                                     get(i)) + ":"
                                             + estimated_fdrscore.get(i).
                                             toString());
                    }
                    break;
                case "Peptide":
                    for (int i = 0; i < sorted_peptideNames.size(); i++) {
                        sorted_stringMap.put(sorted_peptideNames.get(i),
                                             estimated_simpleFDR.get(i).
                                             toString() + ":"
                                             + String.valueOf(estimated_qvalue.
                                                     get(i)) + ":"
                                             + estimated_fdrscore.get(i).
                                             toString());
                    }
                    break;
                default:
                    for (int i = 0; i < sorted_spectrumItem.size(); i++) {
                        sorted_stringMap.put(sorted_spectrumResult.get(i),
                                             estimated_simpleFDR.get(i).
                                             toString() + ":" + String.valueOf(
                                                     estimated_qvalue.get(i)));
                    }
                    break;
            }

            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult sr
                        = iterSpectrumIdentificationResult.next();
                for (SpectrumIdentificationItem sii : sr.
                        getSpectrumIdentificationItem()) {
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
                        Cv cv = MzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);

                        cvParamestimated_simpleFDR.setAccession("MS:1002351");
                        cvParamestimated_simpleFDR.
                                setName("PSM-level local FDR");
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
                        Cv cv = MzidLibUtils.getpsiCV(mzIdentMLUnmarshaller);

                        cvParamestimated_simpleFDR.setAccession("MS:1002359");
                        cvParamestimated_simpleFDR.setName(
                                "distinct peptide-level local FDR");
                        cvParamestimated_simpleFDR.setValue(sii_arr[0]);
                        cvParamestimated_simpleFDR.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_simpleFDR);

                        cvParamestimated_qvalue.setAccession("MS:1001868");
                        cvParamestimated_qvalue.setName(
                                "distinct peptide-level q-value");
                        cvParamestimated_qvalue.setValue(sii_arr[1]);
                        cvParamestimated_qvalue.setCv(cv);

                        sii.getCvParam().add(cvParamestimated_qvalue);

                        cvParamfdrscore.setAccession("MS:1002360");
                        cvParamfdrscore.setName(
                                "distinct peptide-level FDRScore");
                        cvParamfdrscore.setValue(sii_arr[2]);
                        cvParamfdrscore.setCv(cv);
                        sii.getCvParam().add(cvParamfdrscore);

                        // Added by Fawaz Ghali: update to mzid 1.2 20/01/2015
                        CvParam peptidePassesThreshold = new CvParam();
                        peptidePassesThreshold.setAccession("MS:1002500");
                        peptidePassesThreshold.setName(
                                "peptide passes threshold");
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
                Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup
                        = mzIdentMLUnmarshaller
                        .unmarshalCollectionFromXpath(
                                MzIdentMLElement.ProteinAmbiguityGroup);
                while (iterProteinAmbiguityGroup.hasNext()) {
                    ProteinAmbiguityGroup proteinAmbiguityGroup
                            = iterProteinAmbiguityGroup.next();

                    for (ProteinDetectionHypothesis proteinDetectionHypothesis
                            : proteinAmbiguityGroup
                            .getProteinDetectionHypothesis()) {
                        List<CvParam> cvParamListproteinDetectionHypothesis
                                = proteinDetectionHypothesis.getCvParam();
                        for (CvParam cvParam
                                : cvParamListproteinDetectionHypothesis) {
                            String accession = cvParam.getAccession();

                            if (accession.equals("MS:1001591")) {
                                String sortedValue = sorted_stringMap.get(
                                        proteinAmbiguityGroup.getId());
                                CvParam cvParamestimated_simpleFDR
                                        = new CvParam();
                                CvParam cvParamestimated_qvalue = new CvParam();
                                String[] sii_arr = sortedValue.split(":");
                                Cv cv = MzidLibUtils.getpsiCV(
                                        mzIdentMLUnmarshaller);
                                cvParamestimated_simpleFDR.setAccession(
                                        "MS:1002370");
                                cvParamestimated_simpleFDR.setName(
                                        "protein group-level local FDR");
                                cvParamestimated_simpleFDR.setValue(sii_arr[0]);
                                cvParamestimated_simpleFDR.setCv(cv);

                                proteinAmbiguityGroup.getCvParam().add(
                                        cvParamestimated_simpleFDR);

                                cvParamestimated_qvalue.setAccession(
                                        "MS:1002373");
                                cvParamestimated_qvalue.setName(
                                        "protein group-level q-value");
                                cvParamestimated_qvalue.setValue(sii_arr[1]);
                                cvParamestimated_qvalue.setCv(cv);

                                proteinAmbiguityGroup.getCvParam().add(
                                        cvParamestimated_qvalue);
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

            System.out.println("Output written to " + csvFileName);
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className
                    + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    /**
     * Compute FDR score and q value etc using the method described in Jones et
     * al. Proteomics, 2009,9, 1220-1229
     * Parallel computing using Stream
     */
    public void computeFDRusingJonesMethodPar() {
        getEvalueSortedPeptideListPar();
        computeSimpleFDRPar();
        computeQValuesPar();
        computeFDRScorePar();
    }

    /**
     * Sort the data using evalue
     * Make simple arrayLists from complicated data structures returned after
     * reading mzIdentML file. These simple arrayLists will be used for further
     * calculations of various scores.
     * Parallel computing using Stream
     */
    private void getEvalueSortedPeptideListPar() {

        // Call the sorting routine to find the indices of sorted evalues
        TreeSortForIndices sortClass = new TreeSortForIndices();
        Integer[] sortOrderForEvalues = sortClass.sortTheValueColumn(
                fdrgReader.getEvalues().toArray(), betterScoresAreLower);

        // Arrange the values in the order determined by the sorting
        // operation, such that, each index an arraylist can be map to
        // the entries in other arraylists for the "same" index
        Stream.of(sortOrderForEvalues).parallel().forEachOrdered(ind -> {
            sorted_spectrumResult.add(fdrgReader.getSpectrumResult().get(
                    ind));
            sorted_spectrumItem.add(fdrgReader.getSpectrumItem().get(ind));
            sorted_peptideNames.add(fdrgReader.getPeptideNames().get(ind));
            sorted_evalues.add(fdrgReader.getEvalues().get(ind));
            sorted_decoyOrNot.add(fdrgReader.getDecoyOrNot().get(ind));
        });

        // Clear the memory for the items no more needed
        fdrgReader.getSpectrumResult().clear();
        fdrgReader.getPeptideNames().clear();
        fdrgReader.getSpectrumItem().clear();
        fdrgReader.getEvalues().clear();
        fdrgReader.getDecoyOrNot().clear();
    }

    /**
     * Compute simple FDR
     * Parallel computing using Stream
     */
    private void computeSimpleFDRPar() {

        final MutableInt falsePositiveCount = new MutableInt();
        final MutableInt allTargets = new MutableInt();
        final MutableDouble falsePositiveDivRatio = new MutableDouble();
        final MutableDouble simpleFDR = new MutableDouble();

        IntStream.range(0, sorted_peptideNames.size()).parallel()
                .forEachOrdered(idx -> {
                    if (sorted_decoyOrNot.get(idx).equals("true")) {
                        falsePositiveCount.increment();
                    } else {
                        allTargets.increment();
                    }
                    falsePositiveDivRatio.setValue((double) falsePositiveCount.
                            intValue() / decoyRatio);

                    simpleFDR.setValue(falsePositiveDivRatio.doubleValue()
                            / (double) allTargets.intValue());
                    estimated_simpleFDR.add(simpleFDR.doubleValue());
                    estimated_qvalue.add(0d);

                    fp.add(falsePositiveDivRatio.doubleValue());
                    double tpValue = allTargets.intValue()
                            - falsePositiveDivRatio.
                            doubleValue();
                    tp.add(tpValue);
                });

        if (fp.isEmpty() || (fp.get(fp.size() - 1) == 0)) {
            System.out.println(
                    "No decoys found for search engine Mascot|Omssa|tandem - likely caused by: wrong decoy regex, database doesn't contain decoys or the search reported only identifications for stringent e-values - please allow identifications up to e-value = 10 for omssa and tandem. See omssa and tandem documentation for how to do change this setting");
        }
    }

    /**
     * Compute q-value
     * Parallel computing using Stream
     */
    private void computeQValuesPar() {
        double immediateMinFdr = estimated_simpleFDR.get(estimated_simpleFDR.
                size() - 1);
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
     * Parallel computing using Stream
     */
    private void computeFDRScorePar() {
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
        int i = 0;
        for (; i < sorted_peptideNames.size(); i++) {
            double current_evalue = sorted_evalues.get(i);
            double current_qvalue = estimated_qvalue.get(i);

            if (current_qvalue > prev_qvalue) {

                double slope;
                double intercept;

                // Work out the slope and the intercept
                // Tells us which co-ordinates to use for
                if (current_evalue != prev_evalue) {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue
                            - prev_evalue);
                    intercept = prev_qvalue - slope * prev_evalue;
                } else {
                    slope = (current_qvalue - prev_qvalue) / (current_evalue
                            - prev_prev_evalue);
                    intercept = prev_qvalue - slope * prev_prev_evalue;
                }

                if (counter_backwardStep > 0) { // compute the FDR score for flat q-value region
                    for (int k = 0; k <= counter_backwardStep; k++) {
                        int index = i - counter_backwardStep + k;
                        double fdrScore = slope * sorted_evalues.get(index)
                                + intercept;
                        estimated_fdrscore.set(index, fdrScore);
                        lastNonZeroScoreIdx = fdrScore > 0 && index
                                > lastNonZeroScoreIdx ? index
                                        : lastNonZeroScoreIdx;
                    }
                } else { // In case an immediate increment in q value is found
                    double fdrScore = slope * current_evalue + intercept;
                    estimated_fdrscore.set(i, fdrScore);
                    lastNonZeroScoreIdx = fdrScore > 0 && i
                            > lastNonZeroScoreIdx ? i : lastNonZeroScoreIdx;
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
                System.out.println("Last FDR Value Found at " + i
                        + " and pre-recorded version is: " + lastNonZeroScoreIdx);
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
    public List<String> getSorted_spectrumResult() {
        return sorted_spectrumResult;
    }

    @SuppressWarnings("unused")
    public List<String> getSorted_peptideNames() {
        return sorted_peptideNames;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_evalues() {
        return sorted_evalues;
    }

    @SuppressWarnings("unused")
    public List<String> getSorted_decoyOrNot() {
        return sorted_decoyOrNot;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_simpleFDR() {
        return estimated_simpleFDR;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_qValues() {
        return estimated_qvalue;
    }

    @SuppressWarnings("unused")
    public List<Double> getSorted_estimatedFDR() {
        return estimated_fdrscore;
    }

    @SuppressWarnings("unused")
    public List<Double> getTP() {
        return tp;
    }

    @SuppressWarnings("unused")
    public List<Double> getFP() {
        return fp;
    }

    /**
     * Clear all the data structures in the class to release all the memory
     */
    @SuppressWarnings("unused")
    public void clearAllData() {
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
    // Write the sorted data into a file
    public void writeToTsvFile(String fileName)
            throws Exception {

        try (Writer out = new BufferedWriter(new FileWriter(fileName))) {
            String outStr;
            for (int i = 0; i < sorted_peptideNames.size(); i++) {

                outStr = sorted_spectrumResult.get(i) + "\t"
                        + sorted_peptideNames.
                        get(i) + "\t" + sorted_decoyOrNot.get(i)
                        + "\t"
                        + estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.
                        get(i) + "\t" + estimated_fdrscore.get(i)
                        + "\n";

                out.write(outStr);
            }
        }
    }

    // Write the sorted data into a file
    @SuppressWarnings("unused")
    public void writeToCsvFile(String fileName)
            throws Exception {

        try (Writer out = new BufferedWriter(new FileWriter(fileName))) {
            String outStrHead
                    = "sorted_spectrumResult.get(i),sorted_peptideNames.get(i) , sorted_decoyOrNot.get(i) ,  sorted_evalues.get(i).toString() , + sorted_scores.get(i).toString() , estimated_simpleFDR.get(i) , estimated_qvalue.get(i) , estimated_fdrscore.get(i) \n";

            for (int i = 0; i < sorted_peptideNames.size(); i++) {

                String outStr = sorted_spectrumResult.get(i) + ","
                        + sorted_peptideNames.get(i) + ","
                        + sorted_decoyOrNot.get(i) + ","
                        + sorted_evalues.get(i).
                        toString() + ","
                        + estimated_simpleFDR.get(i) + "," + estimated_qvalue.
                        get(i)
                        + "," + estimated_fdrscore.get(i)
                        + "\n";

                out.write(outStr);
            }
        }
    }

    //  Write the sorted data into a file
    @SuppressWarnings("unused")
    public void writeTheSortedDataToFile(String fileName)
            throws Exception {

        // Basic check to see that each peptide has a evalue
        if (sorted_peptideNames.size() != sorted_evalues.size()) {
            throw (new Exception("Number of entries = " + sorted_peptideNames.
                    size()
                    + "in sorted_peptideNames don't match with the number of entries = "
                    + sorted_evalues.size()
                    + "in sorted_evalues"));
        }

        try (Writer out = new BufferedWriter(new FileWriter(fileName))) {
            for (int i = 0; i < sorted_peptideNames.size(); i++) {
                String outStr = sorted_spectrumResult.get(i) + "\t"
                        + sorted_spectrumItem.get(i) + "\t"
                        + sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.
                        get(
                                i) + "\t"
                        + sorted_evalues.get(i).toString() + "\n";

                out.write(outStr);
            }
        }
    }

    private UserParam makeUserParam(String name, String value) {
        UserParam userParam = new UserParam();
        userParam.setName(name);
        userParam.setValue(value);
        return userParam;
    }

}
