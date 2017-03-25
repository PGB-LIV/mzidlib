package uk.ac.liv.mzidlib;

import uk.ac.liv.mzidlib.util.MzidLibUtils;
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

import javax.xml.bind.JAXBException;

import org.apache.commons.lang.StringEscapeUtils;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author jonesar
 */
public class ThresholdMzid {

    private String inputFile = "resources/55merge_omssa.mzid";
    private String outputFile = "test_out.mzid";
    private boolean psmThreshold = true;
    private String cvAccessionScoreThreshold = "MS:1001328";
    private double thresholdValue = 1.0;
    private boolean scoresGoLowToHigh = true;
    // Removed by FG to avoid out of memory
    // private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private boolean deleteUnderThreshold = false;
    private Map<String, DBSequence> dbSequenceHashMap = new HashMap<>();
    private Map<String, Peptide> peptideHashMap = new HashMap<>();
    private Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
    private Cv psiCV;
    private MzidLibUtils utils = new MzidLibUtils();        //Utils class containing makeCvParam methods
    private String thresholdCvName = null;      //to be set while reading the file

    // Added by Fawaz Ghali 14/5/2013
    private String scoreLevel = "PDH"; //PAG or PDH

    /*
     * For testing only
     */
    public static void main(String[] args) {

        //Boolean requireSIIsToPassThreshold = Boolean.parseBoolean(args[0]);
        //TODO - tidy up command line arguments etc.
        ThresholdMzid thresholdMzid = new ThresholdMzid();

    }

    /*
     * Placeholder for code to add thresholds to SpectrumIdentificationItem or
     * ProteinDetectionHypothesis
     *
     */
    public ThresholdMzid() {

        System.out.println("Running with defaults for testing");
        this.init();

    }

    /*
     * Placeholder for code to add thresholds to SpectrumIdentificationItem or
     * ProteinDetectionHypothesis
     *
     */
    public ThresholdMzid(String inFile, String outFile, boolean psmThresh, String cvAccScoreThreshold, double threshValue, boolean scoresLowHigh, boolean deleteUnderThreshold, String scoreLevel) {

        inputFile = inFile;
        outputFile = outFile;
        this.psmThreshold = psmThresh;
        this.cvAccessionScoreThreshold = cvAccScoreThreshold;
        this.thresholdValue = threshValue;
        this.scoresGoLowToHigh = scoresLowHigh;
        this.deleteUnderThreshold = deleteUnderThreshold;
        this.scoreLevel = scoreLevel;
        this.init();
    }


    @Deprecated
    public ThresholdMzid(String inFile, String outFile, boolean psmThresh, String cvAccScoreThreshold, double threshValue, boolean scoresLowHigh, boolean deleteUnderThreshold) {
        this(inFile, outFile, psmThresh, cvAccScoreThreshold, threshValue, scoresLowHigh, deleteUnderThreshold, "");
    }

    private void init() {

        if (psmThreshold) {
            System.out.print("Setting pass threshold to false for PSMs ");
        } else {

            System.out.print("Setting pass threshold to false for PDHs ");
        }

        if (scoresGoLowToHigh) {
            System.out.println(cvAccessionScoreThreshold + " > " + thresholdValue);
        } else {
            System.out.println(cvAccessionScoreThreshold + " < " + thresholdValue);
        }

        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputFile));

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        // Removed by FG
//        if(psmThreshold){
//            setPSMThreshold();
//        }
//        else{
//            
//            setPDHThreshold();
//        }
//        
        writeToMzIdentMLFile();
    }

//    private void setPSMThreshold(){
//        
//        List<SpectrumIdentificationList> siListList = mzIdentML.getDataCollection().getAnalysisData().getSpectrumIdentificationList();
//        for(SpectrumIdentificationList siList : siListList){
//            for(SpectrumIdentificationResult sir : siList.getSpectrumIdentificationResult()){
//                for(SpectrumIdentificationItem sii : sir.getSpectrumIdentificationItem()){
//                    for(CvParam cvParam : sii.getCvParam()){
//                        if(cvParam.getAccession().equals(cvAccessionScoreThreshold) ){
//                            if((Double.parseDouble(cvParam.getValue()) <= thresholdValue && scoresGoLowToHigh) || (Double.parseDouble(cvParam.getValue()) >= thresholdValue && !scoresGoLowToHigh)){
//                                sii.setPassThreshold(true);                               
//                            }
//                            else{
//                                sii.setPassThreshold(false); 
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    private void setPDHThreshold(){
//        
//        for(ProteinAmbiguityGroup pag : mzIdentML.getDataCollection().getAnalysisData().getProteinDetectionList().getProteinAmbiguityGroup()){
//            for(ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()){
//                for(CvParam cvParam : pdh.getCvParam()){
//                      if(cvParam.getAccession().equals(cvAccessionScoreThreshold) ){
//                          //System.out.println("Found:" + cvParam.getName());
//                            if((Double.parseDouble(cvParam.getValue()) <= thresholdValue && scoresGoLowToHigh) || (Double.parseDouble(cvParam.getValue()) >= thresholdValue && !scoresGoLowToHigh)){
//                                    pdh.setPassThreshold(true);
//                            }
//                            else{
//                                 pdh.setPassThreshold(false);
//                            }
//                      }
//                }
//            }
//        }
//    }
    // Write the new data into a file
    private void writeToMzIdentMLFile() {
        try {
            String outFile = outputFile;
//            if (!outFile.endsWith(".mzid")) {
//                outFile = outFile + ".mzid";
//            }
            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            //cvList = mzIdentML.getCvList();
            CvList cvList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
            Iterator<Cv> iterCv = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                }
            }
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            AnalysisSoftwareList analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            AuditCollection auditCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            Provider provider = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            AnalysisProtocolCollection analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            AnalysisCollection analysisCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);

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
            Cv psiCV;
            psiCV = utils.getpsiCV(mzIdentMLUnmarshaller);
            param.setParam(utils.makeCvParam("MS:1002237", "mzidLib", psiCV));
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

            //SequenceCollection sequenceCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Iterator<SpectrumIdentificationResult> sirIter = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (sirIter.hasNext()) {

                SpectrumIdentificationResult sr = sirIter.next();
                List<SpectrumIdentificationItem> siiList = sr.getSpectrumIdentificationItem();
                Iterator<SpectrumIdentificationItem> siiIter = siiList.iterator();
                while (siiIter.hasNext()) {
                    SpectrumIdentificationItem spectrumIdentificationItem = siiIter.next(); // must be called before you can call i.remove()
                    if (psmThreshold) {
                        boolean exit = true;
                        for (CvParam cvParam : spectrumIdentificationItem.getCvParam()) {
                            if (cvParam.getAccession().equals(cvAccessionScoreThreshold)) {
                                exit = false;

                                if (thresholdCvName == null) {
                                    thresholdCvName = cvParam.getName();
                                }

                                if ((Double.parseDouble(cvParam.getValue()) <= thresholdValue && scoresGoLowToHigh) || (Double.parseDouble(cvParam.getValue()) >= thresholdValue && !scoresGoLowToHigh)) {
                                    spectrumIdentificationItem.setPassThreshold(true);
                                } else {
                                    spectrumIdentificationItem.setPassThreshold(false);

                                }
                            }
                        }
                        if (exit) {
                            System.out.println("Bad CVParam i.e one that is not found in the file");
                        }
                    }
                    // delete under threshold
                    if (deleteUnderThreshold && !spectrumIdentificationItem.isPassThreshold()) {
                        siiIter.remove();
                    } else {
                        List<PeptideEvidenceRef> PeptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                        for (PeptideEvidenceRef peptideEvidenceRef : PeptideEvidenceRefList) {
                            if (peptideEvidenceRef != null) {

                                PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                                if (peptideEvidence != null) {
                                    peptideEvidenceHashMap.put(peptideEvidence.getId(), peptideEvidence);
                                    DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, peptideEvidence.getDBSequenceRef());
                                    if (dBSequence != null) {
                                        dbSequenceHashMap.put(dBSequence.getId(), dBSequence);
                                    }
                                    String ref = StringEscapeUtils.escapeXml(peptideEvidence.getPeptideRef());
                                    Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, ref);
                                    if (peptide != null) {
                                        peptideHashMap.put(peptide.getId(), peptide);
                                    }
                                    String peptideref = StringEscapeUtils.escapeXml(spectrumIdentificationItem.getPeptideRef());
                                    Peptide peptideSII = mzIdentMLUnmarshaller.unmarshal(Peptide.class, peptideref);

                                    if (peptideSII != null) {
                                        peptideHashMap.put(peptideSII.getId(), peptideSII);
                                    }
                                }
                            }
                        }
                    }

                }
                if (!siiList.isEmpty()) {
                    siList.getSpectrumIdentificationResult().add(sr);
                }
            }

            SequenceCollection sequenceCollection = new SequenceCollection();

            Iterator<Entry<String, DBSequence>> itDbSeq = dbSequenceHashMap.entrySet().iterator();
            List<DBSequence> dbSeqList = new ArrayList<>();
            while (itDbSeq.hasNext()) {
                Entry<String, DBSequence> pairs = itDbSeq.next();
                dbSeqList.add(pairs.getValue());
                itDbSeq.remove();
            }

            Iterator<Entry<String, Peptide>> itPeptide = peptideHashMap.entrySet().iterator();
            List<Peptide> peptideList = new ArrayList<>();
            while (itPeptide.hasNext()) {
                Entry<String, Peptide> pairs = itPeptide.next();
                peptideList.add(pairs.getValue());
                itPeptide.remove();
            }

            Iterator<Entry<String, PeptideEvidence>> itPeptideEvidence = peptideEvidenceHashMap.entrySet().iterator();
            List<PeptideEvidence> peptideEvidenceList = new ArrayList<>();
            while (itPeptideEvidence.hasNext()) {
                Entry<String, PeptideEvidence> pairs = itPeptideEvidence.next();
                peptideEvidenceList.add(pairs.getValue());
                itPeptideEvidence.remove();
            }

            sequenceCollection.getPeptideEvidence().addAll(peptideEvidenceList);
            sequenceCollection.getDBSequence().addAll(dbSeqList);
            sequenceCollection.getPeptide().addAll(peptideList);

            String proteinDetectionListRef = "";
            if (analysisCollection.getProteinDetection() != null) {
                proteinDetectionListRef = analysisCollection.getProteinDetection().getProteinDetectionListRef();
            }
            ProteinDetectionList pdList;
                pdList = new ProteinDetectionList();

                pdList.setId(proteinDetectionListRef);
                        
            Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {

                ProteinAmbiguityGroup pag = iterProteinAmbiguityGroup.next();
                boolean setPassThreshold = false;
                if (!psmThreshold && scoreLevel.equals("PDH")) {

                    for (ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()) {
                        boolean exit = true;
                        for (CvParam cvParam : pdh.getCvParam()) {
                            if (cvParam.getAccession().equals(cvAccessionScoreThreshold)) {
                                exit = false;

                                if (thresholdCvName == null) {
                                    thresholdCvName = cvParam.getName();
                                }

                                if ((Double.parseDouble(cvParam.getValue()) <= thresholdValue && scoresGoLowToHigh) || (Double.parseDouble(cvParam.getValue()) >= thresholdValue && !scoresGoLowToHigh)) {
                                    pdh.setPassThreshold(true);
                                    setPassThreshold = true;
                                } else {
                                    pdh.setPassThreshold(false);
                                    setPassThreshold = false;
                                }
                            }
                        }
                        if (exit) {
                            System.out.println("Bad CVParam i.e one that is not found in the file");
                        }
                    }
                } else if (!psmThreshold && scoreLevel.equals("PAG")) {

                    for (CvParam cvParam : pag.getCvParam()) {
                        if (cvParam.getAccession().equals(cvAccessionScoreThreshold)) {

                            if (thresholdCvName == null) {
                                thresholdCvName = cvParam.getName();
                            }

                            if ((Double.parseDouble(cvParam.getValue()) <= thresholdValue && scoresGoLowToHigh) || (Double.parseDouble(cvParam.getValue()) >= thresholdValue && !scoresGoLowToHigh)) {
                                setPassThreshold = true;
                            }
                        }
                    }
                    String anchorProteinDetectionHypothesis = null;
                    boolean anchorProtein = false;
                    for (int j = 0; j < pag.getProteinDetectionHypothesis().size(); j++) {
                        ProteinDetectionHypothesis proteinDetectionHypothesis = pag.getProteinDetectionHypothesis().get(j);
                        List<CvParam> cvParamListproteinDetectionHypothesis = proteinDetectionHypothesis.getCvParam();
                        for (int s = 0; s < cvParamListproteinDetectionHypothesis.size(); s++) {
                            CvParam cvParam = cvParamListproteinDetectionHypothesis.get(s);
                            String accession = cvParam.getAccession();
                            if (accession.equals("MS:1001591")) {
                                anchorProteinDetectionHypothesis = proteinDetectionHypothesis.getId();
                                anchorProtein = true;
                                break;
                            }
                        }

                    }
                    for (int j = 0; j < pag.getProteinDetectionHypothesis().size(); j++) {
                        ProteinDetectionHypothesis proteinDetectionHypothesis = pag.getProteinDetectionHypothesis().get(j);
                        if (anchorProtein && anchorProteinDetectionHypothesis.equals(proteinDetectionHypothesis.getId())) {
                            proteinDetectionHypothesis.setPassThreshold(setPassThreshold);
                        } else if (!anchorProtein && j == 0) {
                            proteinDetectionHypothesis.setPassThreshold(setPassThreshold);
                        } else {
                            proteinDetectionHypothesis.setPassThreshold(false);
                        }

                    }
                }
                if (!deleteUnderThreshold || setPassThreshold) {
                    pdList.getProteinAmbiguityGroup().add(pag);
                }
            }

            if (!deleteUnderThreshold) {
                SequenceCollection sequenceCollection1 = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
                if (sequenceCollection1 != null) {
                    marshaller.marshal(sequenceCollection1, writer);
                }
            } else {

                if (sequenceCollection != null) {
                    marshaller.marshal(sequenceCollection, writer);
                }
            }
            writer.write("\n");

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {

                if (psmThreshold) {
                    SpectrumIdentificationProtocol protocol = analysisProtocolCollection.getSpectrumIdentificationProtocol().get(0);
                    ParamList paramList = protocol.getThreshold();
                    paramList.getCvParam().clear();

                    if (thresholdCvName != null) {
                        paramList.getCvParam().add(utils.makeCvParam(cvAccessionScoreThreshold, thresholdCvName, psiCV, "" + thresholdValue));
                    }
                } else {
                    ProteinDetectionProtocol protocol = analysisProtocolCollection.getProteinDetectionProtocol();
                    ParamList paramList = protocol.getThreshold();
                    paramList.getCvParam().clear();
                    //paramList.getCvParam().add(utils.makeCvParam(cvAccessionScoreThreshold, "name", psiCV, "" + thresholdValue));
                    if (thresholdCvName != null) {
                        paramList.getCvParam().add(utils.makeCvParam(cvAccessionScoreThreshold, thresholdCvName, psiCV, "" + thresholdValue));
                    }
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
            // Added by FG 12/11/2014 avoiding empty SIL, invalid mzid file
            if (siList.getSpectrumIdentificationResult().isEmpty()) {
                throw new RuntimeException("SIL is empty!");
            }
            marshaller.marshal(siList, writer);
            writer.write("\n");

            if (!pdList.getProteinAmbiguityGroup().isEmpty()){
                marshaller.marshal(pdList, writer);
                writer.write("\n");
            }
            
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (JAXBException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 04 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

}
