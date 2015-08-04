package uk.ac.liv.mzidlib;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.compare.CompareDBSequence;
import uk.ac.liv.mzidlib.compare.ComparePeptide;
import uk.ac.liv.mzidlib.compare.ComparePeptideEvidence;

/**
 *
 * @author Fawaz Ghali
 */
public class CombinePSMMzidFiles {

    private Set<SpectrumIdentificationList> combinedSpectrumIdentificationList = new HashSet<>();

    private Map<String, DBSequence> outputDBSeqMap = new HashMap<>();
    private Map<String, Peptide> outputPeptideMap = new HashMap<>();
    //private Map<String, Set<PeptideEvidence>> pepEMap = new HashMap<>();

//    private Map<String, DBSequence> dbSeqOldNewMap = new HashMap<>();
//    private Map<String, Peptide> pepOldNewMap = new HashMap<>();
    //private Map<String, PeptideEvidence> pepEOldNewMap = new HashMap<>();
    private Map<String, PeptideEvidence> inputPEMap = new HashMap<>();
    private Map<String, PeptideEvidence> outputPEMap = new HashMap<>();

    //metadata
    private AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();
    private AuditCollection auditCollection = new AuditCollection();
    private Provider provider = new Provider();
    private AnalysisProtocolCollection analysisProtocolCollection = new AnalysisProtocolCollection();
    private CvList cvList = new CvList();
    private AnalysisCollection combinedAnalysisCollection = new AnalysisCollection();
    private Set<SpectrumIdentification> combinedSpectrumIdentification = new HashSet<>();
    private Set<SpectraData> combinedSpectraData = new HashSet<>();
    private Inputs combinedInputs = new Inputs();
    private AnalysisData combinedAnalysisData = new AnalysisData();
    private final boolean combineFractions;
    long startTime, stopTime,elapsedTime;

    private boolean verbose = true;

    public CombinePSMMzidFiles(String inputs, String output, boolean combineFractions) {

        this.combineFractions = combineFractions;
        File file = new File(inputs);
        startTime = System.currentTimeMillis();
        if (file.isDirectory()) {
            File[] mzidFiles = file.listFiles();
            int j = 0;
            for (int i = 0; i < mzidFiles.length; i++) {
                if (mzidFiles[i].getAbsolutePath().endsWith(".mzid")) {
                    MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzidFiles[i].getAbsolutePath()));
                    System.out.println("Reading " + mzidFiles[i].getAbsolutePath());
                    readMzid(mzIdentMLUnmarshaller, j);
                    j = j + 1;
                    System.out.println("Reading " + mzidFiles[i].getAbsolutePath() + " completed.");
                }
            }
        } else {
            String[] mzidFiles = inputs.split(";");
            int j = 0;
            for (int i = 0; i < mzidFiles.length; i++) {

                MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzidFiles[i]));
                System.out.println("Reading " + mzidFiles[i]);
                readMzid(mzIdentMLUnmarshaller, j);
                j = j + 1;
                System.out.println("Reading " + mzidFiles[i] + " completed.");

            }

        }
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("-----------------------------------------------");
            System.out.println("Reading and combining time " + elapsedTime / 1000 + " Seconds");
            System.out.println("-----------------------------------------------");
        }
        startTime = System.currentTimeMillis();
        System.out.println("Writing " + output);
        writeMzidFile(output);
        System.out.println("Writing " + output + " completed.");
        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            System.out.println("-----------------------------------------------");
            System.out.println("Writing time " + elapsedTime / 1000 + " Seconds");
            System.out.println("-----------------------------------------------");        }
    }

    private void readMzid(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, int i) {
        try {
            combineMetaData(mzIdentMLUnmarshaller, i);
            combineSequenceCollection(mzIdentMLUnmarshaller, i);
            combineSpectrumIdentification(mzIdentMLUnmarshaller, i);

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    public void writeMzidFile(String csvFileName) {
        try {
            String outFile = csvFileName;
            Writer writer = new FileWriter(outFile);
            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();
            writer.write(marshaller.createXmlHeader() + "\n");
            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");

            marshaller.marshal(cvList, writer);
            writer.write("\n");

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

            sequenceCollection.getDBSequence().addAll(outputDBSeqMap.values());
            sequenceCollection.getPeptide().addAll(outputPeptideMap.values());
            sequenceCollection.getPeptideEvidence().addAll(outputPEMap.values());

            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");
            ProteinDetectionList pdl = null;
            if (combinedAnalysisCollection != null) {
                ProteinDetection pd = new ProteinDetection();
                pd.setId("PD_1");
                pd.setProteinDetectionProtocol(analysisProtocolCollection.getProteinDetectionProtocol());
                pdl = new ProteinDetectionList();
                pdl.setId("PDL_1");
                pd.setProteinDetectionList(pdl);
                InputSpectrumIdentifications isi = new InputSpectrumIdentifications();
                isi.setSpectrumIdentificationList(combinedAnalysisData.getSpectrumIdentificationList().get(0));
                pd.getInputSpectrumIdentifications().add(isi);
                combinedAnalysisCollection.setProteinDetection(pd);
                marshaller.marshal(combinedAnalysisCollection, writer);
            }
            writer.write("\n");

            marshaller.marshal(analysisProtocolCollection, writer);
            writer.write("\n");

            writer.write(marshaller.createDataCollectionStartTag() + "\n");

            writer.write("\n");

            if (combinedInputs != null) {
                marshaller.marshal(combinedInputs, writer);
            }
            writer.write("\n");

            if (combinedAnalysisData != null) {
                combinedAnalysisData.setProteinDetectionList(pdl);
                marshaller.marshal(combinedAnalysisData, writer);

            }

            // writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
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

    public static void main(String args[]) {

        CombinePSMMzidFiles CombineMzidFiles = new CombinePSMMzidFiles(args[0], args[1], Boolean.valueOf(args[2]));
    }

    private void combineMetaData(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, int i) {
        //metadata
        //cvList, analysisSoftware List, auditCollection and provider should be identical already, so we can use the first instance found and discard the rest.
        if (i == 0) {
            cvList = (CvList) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
            auditCollection = (AuditCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            analysisSoftwareList = (AnalysisSoftwareList) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            analysisProtocolCollection = (AnalysisProtocolCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            provider = (Provider) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);

        }
    }

    private void combineSequenceCollection(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, int i) {
        /*
         The entire algorithm for this part can be as follows.

         Have one map:

         <String accession, DBSequence object>

         When reading a new file, check if accession is present, 

         if yes: retrieve the DBSequence id (retrieved ID)value from the referenced object
         -- Then check whether DBSequence retrieved id is different from current id
         ---- In almost all cases it will be the same, no more work to do

         If no: simply add to the Map, no more work to do

         If current_id doesn't equal retrieved_id,

         then you need to populate dbSeqOldNewMap. 

         In many cases (probably all current csaes), dbSeqOldNewMap will be empty, thus check if empty before trying to use, if empty no work to do.

         */

        // Only for OMSSA and Tandem
        Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        while (iterDBSequence.hasNext()) {
            DBSequence dBSequence = iterDBSequence.next();
            String accession = dBSequence.getAccession();
            if (!outputDBSeqMap.containsKey(accession)) {
                outputDBSeqMap.put(accession, dBSequence);
            }
        }

        Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
        while (iterPeptide.hasNext()) {
            Peptide peptide = iterPeptide.next();
            String key = peptide.getPeptideSequence() + "_" + peptide.getId();
            if (!outputPeptideMap.containsKey(key)) {
                outputPeptideMap.put(key, peptide);
            }
        }

        inputPEMap.clear();
        Iterator<PeptideEvidence> iterPeptideEvidence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
        while (iterPeptideEvidence.hasNext()) {
            PeptideEvidence pe = iterPeptideEvidence.next();
            inputPEMap.put(pe.getId(), pe);
        }

    }

    private void combineSpectrumIdentification(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, int i) {
        AnalysisCollection ac = (AnalysisCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
        List<SpectrumIdentification> spectrumIdentificationList = ac.getSpectrumIdentification();

        // Handle inputs by combining spectraData after changing the id
        Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
        List<SpectraData> spectraDataList = inputs.getSpectraData();
        for (SpectraData spectraData : spectraDataList) {
            String newKey = spectraData.getId() + "_combined_" + i;
            spectraData.setId(newKey);
            combinedSpectraData.add(spectraData);
            if (i == 0) {
                combinedInputs.getSourceFile().addAll(inputs.getSourceFile());
                combinedInputs.getSearchDatabase().addAll(inputs.getSearchDatabase());
            }
        }
        combinedInputs.getSpectraData().addAll(combinedSpectraData);

        // Handle AnalysisCollection
        for (SpectrumIdentification spectrumIdentification : spectrumIdentificationList) {
            String newId = spectrumIdentification.getId() + "_combined_" + i;
            spectrumIdentification.setId(newId);

            // Handle InputSpectra 
            List<InputSpectra> inputSpectraList = spectrumIdentification.getInputSpectra();
            for (InputSpectra inputSpectra : inputSpectraList) {
                String newKeyRef = inputSpectra.getSpectraDataRef() + "_combined_" + i;
                for (SpectraData spectraData : combinedSpectraData) {
                    if (spectraData.getId().equals(newKeyRef)) {
                        inputSpectra.setSpectraData(spectraData);
                        break;
                    }
                }
            }

            // Handle SpectrumIdentificationList
            String newSilID = spectrumIdentification.getSpectrumIdentificationListRef() + "_combined_" + i;
            AnalysisData ad = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisData);

            List<SpectrumIdentificationList> spectrumIdentificationListList = ad.getSpectrumIdentificationList();
            for (SpectrumIdentificationList spectrumIdentificationList1 : spectrumIdentificationListList) {
                String newID = spectrumIdentificationList1.getId() + "_combined_" + i;
                spectrumIdentificationList1.setId(newID);

                FragmentationTable fragmentationTable = spectrumIdentificationList1.getFragmentationTable();
                if (fragmentationTable != null) {
                    List<Measure> measureList = fragmentationTable.getMeasure();
                    for (Measure measure : measureList) {
                        String measureNewID = measure.getId() + "_combined_" + i;
                        measure.setId(measureNewID);
                    }
                }
                List<SpectrumIdentificationResult> spectrumIdentificationResultList = spectrumIdentificationList1.getSpectrumIdentificationResult();
                for (SpectrumIdentificationResult spectrumIdentificationResult : spectrumIdentificationResultList) {
                    String newSIRID = spectrumIdentificationResult.getId() + "_combined_" + i;
                    spectrumIdentificationResult.setId(newSIRID);
                    String newSpectraDataRef = spectrumIdentificationResult.getSpectraDataRef() + "_combined_" + i;
                    for (SpectraData spectraData : combinedSpectraData) {
                        if (spectraData.getId().equals(newSpectraDataRef)) {
                            spectrumIdentificationResult.setSpectraData(spectraData);
                            break;
                        }

                    }
                    List<SpectrumIdentificationItem> spectrumIdentificationItemList = spectrumIdentificationResult.getSpectrumIdentificationItem();
                    for (SpectrumIdentificationItem spectrumIdentificationItem : spectrumIdentificationItemList) {
                        String newSIIID = spectrumIdentificationItem.getId() + "_combined_" + i;
                        spectrumIdentificationItem.setId(newSIIID);

                        // Update old PE with new PE
                        List<PeptideEvidenceRef> oldPeptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                        List<PeptideEvidenceRef> newPeptideEvidenceRefList = new ArrayList<PeptideEvidenceRef>();
                        for (int j = 0; j < oldPeptideEvidenceRefList.size(); j++) {
                            PeptideEvidenceRef peptideEvidenceRef = oldPeptideEvidenceRefList.get(j);
                            PeptideEvidence pe = inputPEMap.get(peptideEvidenceRef.getPeptideEvidenceRef());

                            String key = pe.getPeptideRef() + "_" + pe.getDBSequenceRef() + "_" + pe.getStart() + "_" + pe.getEnd();
                            pe.setId(key);

                            if (!outputPEMap.containsKey(key)) {
                                outputPEMap.put(key, pe);

                            } else {
                                pe = outputPEMap.get(key);
                            }

                            PeptideEvidenceRef per = new PeptideEvidenceRef();
                            per.setPeptideEvidence(pe);
                            newPeptideEvidenceRefList.add(per);

                        }

                        spectrumIdentificationItem.getPeptideEvidenceRef().clear();
                        spectrumIdentificationItem.getPeptideEvidenceRef().addAll(newPeptideEvidenceRefList);
                    }
                }
                if (newID.equals(newSilID)) {
                    spectrumIdentification.setSpectrumIdentificationList(spectrumIdentificationList1);
                }
                combinedAnalysisData.getSpectrumIdentificationList().add(spectrumIdentificationList1);
            }
            combinedSpectrumIdentification.add(spectrumIdentification);
        }
        if (combineFractions) {
            combinedAnalysisCollection.getSpectrumIdentification().addAll(combinedSpectrumIdentification);
        } else {
            if (combinedAnalysisCollection.getSpectrumIdentification() == null || combinedAnalysisCollection.getSpectrumIdentification().size() == 0) {
                combinedAnalysisCollection.getSpectrumIdentification().add(combinedSpectrumIdentification.iterator().next());
            } else {
                for (SpectrumIdentification si : combinedSpectrumIdentification) {
                    combinedAnalysisCollection.getSpectrumIdentification().get(0).getInputSpectra().addAll(si.getInputSpectra());

                }
            }
        }

        combinedSpectraData.clear();
        combinedSpectrumIdentification.clear();

    }
}
