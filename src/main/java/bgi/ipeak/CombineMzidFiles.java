package bgi.ipeak;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author duchaoqin
 */
public class CombineMzidFiles {

    private HashMap<String, DBSequence> dBSequenceHashMap = new HashMap<String, DBSequence>();
    private HashMap<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<String, PeptideEvidence>();
    private HashMap<String, Peptide> peptideHashMap = new HashMap<String, Peptide>();
    private HashMap<String, SpectraData> id_spectraData = new HashMap<String, SpectraData>();
    //metadata
    private List<CvList> cvList = new ArrayList<CvList>();
    private List<AnalysisSoftwareList> analysisSoftwareList = new ArrayList<AnalysisSoftwareList>();
    private List<Provider> provider = new ArrayList<Provider>();
    private List<AuditCollection> auditCollection = new ArrayList<AuditCollection>();
    private List<SequenceCollection> sequenceCollections = new ArrayList<SequenceCollection>();
    private List<AnalysisCollection> analysisCollection = new ArrayList<AnalysisCollection>();
    private List<AnalysisProtocolCollection> analysisProtocolCollection = new ArrayList<AnalysisProtocolCollection>();
    private List<DataCollection> dataCollections = new ArrayList<DataCollection>();
    private List<SpectrumIdentificationList> combinedSpectrumIdentificationList = new ArrayList<SpectrumIdentificationList>();

    //id change infor;
    HashMap<String, DBSequence> db_id_replace = new HashMap<String, DBSequence>();
    HashMap<String, Peptide> pep_id_replace = new HashMap<String, Peptide>();
    HashMap<String, PeptideEvidence> peptideEvidence_id_replace = new HashMap<String, PeptideEvidence>();

    private Vector<String> list_file;
    private HashMap<String, Integer> name_id;

    @Option(name = "-dir", required = true, usage = "(required),input dir or input list file(text file:one line one mzid file)")
    private String input_dir;
    @Option(name = "-out", required = true, usage = "(required),out file path ")
    private String output;
    private String threshold_score_acc = "MS:1001491";
    private Double threshold_score = 1.0;

    public static void main(String args[]) throws IOException {

        CombineMzidFiles com = new CombineMzidFiles();
        CmdLineParser parser = new CmdLineParser(com);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-11)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.CombineMzidFiles [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        com.get_file_list();
        com.combined_mzid();
    }

    public CombineMzidFiles(Vector<String> list_file, String output, HashMap<String, Integer> name_id,
            Double threshold_score, String threshold_score_cv) {
        this.list_file = list_file;
        this.output = output;
        this.name_id = name_id;
        this.threshold_score_acc = threshold_score_cv;
        this.threshold_score = threshold_score;
    }

    public CombineMzidFiles(Vector<String> list_file, String output, HashMap<String, Integer> name_id) {
        this.list_file = list_file;
        this.output = output;
        this.name_id = name_id;
    }

    public CombineMzidFiles() {

    }

    private void get_file_list() throws IOException {
        File in_file = new File(input_dir);
        if (in_file.isDirectory()) {
            File[] mzid_lis = in_file.listFiles();
            list_file = new Vector<String>();
            name_id = new HashMap<String, Integer>();
            Integer id = 0;
            for (File mzid : mzid_lis) {
                list_file.add(mzid.getAbsolutePath());
                name_id.put(mzid.getName(), id);
            }
        } else if (in_file.isFile()) {
            BufferedReader reader = new BufferedReader(new FileReader(in_file));
            String line;
            list_file = new Vector<String>();
            name_id = new HashMap<String, Integer>();
            while ((line = reader.readLine()) != null) {
                String[] infor = line.split("\t");
                list_file.add(infor[0]);
                File list = new File(infor[0]);
                name_id.put(list.getName(), Integer.valueOf(infor[1]));
            }
            reader.close();
        }
    }

    public void combined_mzid() {
        Vector<File> mzidFiles = new Vector<File>();
        for (String file : list_file) {
            File the_file = new File(file);
            mzidFiles.add(the_file);
        }
        for (File mzfile : mzidFiles) {
            Integer i = name_id.get(mzfile.getName());
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzfile.getAbsolutePath()));
            System.out.println("Reading " + mzfile.getAbsolutePath() + " File id=" + i);
            readMzid(mzIdentMLUnmarshaller, i);
            System.out.println("Reading " + mzfile.getAbsolutePath() + " completed.");
        }
        System.out.println("Writing " + output);
        writeMzidFile(output);
        System.out.println("Writing " + output + " completed.");
    }

    private void readMzid(MzIdentMLUnmarshaller mzIdentMLUnmarshaller, int i) {
        try {

            cvList.add((CvList) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList));

            analysisSoftwareList.add((AnalysisSoftwareList) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList));

            provider.add((Provider) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider));

            auditCollection.add((AuditCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection));

            SequenceCollection sc = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
            List<DBSequence> iterDBSequence = sc.getDBSequence();
            for (DBSequence dBSequence : iterDBSequence) {
                String old_id = dBSequence.getId();
                String db_newID = "dbseq_" + dBSequence.getAccession();
                dBSequence.setId(db_newID);
                if (!dBSequenceHashMap.containsKey(dBSequence.getId())) {
                    dBSequenceHashMap.put(dBSequence.getId(), dBSequence);
                }
                db_id_replace.put("com_" + i + "_" + old_id, dBSequenceHashMap.get(dBSequence.getId()));
            }

            List<Peptide> iterPeptide = sc.getPeptide();
            for (Peptide peptide : iterPeptide) {
                String pep_old_ID = peptide.getId();
                String pep_new_ID = peptide.getPeptideSequence();
                peptide.getSubstitutionModification();
                List<Modification> pep_modList = peptide.getModification();
                for (Modification modification : pep_modList) {
                    if (modification.getMonoisotopicMassDelta() != null) {
                        pep_new_ID += "_" + String.valueOf(modification.getMonoisotopicMassDelta()) + "@" + String.valueOf(modification.getLocation());
                    } else if (modification.getAvgMassDelta() != null) {
                        pep_new_ID += "_" + String.valueOf(modification.getAvgMassDelta()) + "@" + String.valueOf(modification.getLocation());
                    }
                }
                peptide.setId(pep_new_ID);
                if (!peptideHashMap.containsKey(pep_new_ID)) {
                    peptideHashMap.put(peptide.getId(), peptide);
                }
                pep_id_replace.put("com_" + i + "_" + pep_old_ID, peptideHashMap.get(pep_new_ID));

            }
            List<PeptideEvidence> iterPeptideEvidence = sc.getPeptideEvidence();
            for (PeptideEvidence peptideEvidence : iterPeptideEvidence) {
                peptideEvidence.setDBSequence(db_id_replace.get("com_" + i + "_" + peptideEvidence.getDBSequenceRef()));

                peptideEvidence.setPeptide(pep_id_replace.get("com_" + i + "_" + peptideEvidence.getPeptideRef()));
                String pe_old_id = peptideEvidence.getId();
                String pe_new_id = peptideEvidence.getDBSequenceRef() + "_" + peptideEvidence.getPeptideRef();
                peptideEvidence.setId(pe_new_id);
                if (!peptideEvidenceHashMap.containsKey(pe_new_id)) {
                    peptideEvidenceHashMap.put(peptideEvidence.getId(), peptideEvidence);
                }
                peptideEvidence_id_replace.put("com_" + i + "_" + pe_old_id, peptideEvidenceHashMap.get(pe_new_id));
            }
            sequenceCollections.add(sc);

            AnalysisCollection ac = (AnalysisCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            List<SpectrumIdentification> spectrum_identification = ac.getSpectrumIdentification();
            for (int k = 0; k < spectrum_identification.size(); k++) {
                SpectrumIdentification si = spectrum_identification.get(k);
                si.setId("com_" + i + "_" + si.getId());
            }
            analysisCollection.add(ac);

            analysisProtocolCollection.add((AnalysisProtocolCollection) mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection));

            DataCollection dc = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.DataCollection);
            List<SpectraData> spectraData_id = dc.getInputs().getSpectraData();
            for (SpectraData spectraData : spectraData_id) {
                String sd_old_id = spectraData.getId();
                String sd_id = "SID_" + i;
                spectraData.setId(sd_id);
                id_spectraData.put("com_" + i + "_" + sd_old_id, spectraData);
            }
            List<SpectrumIdentificationList> spectrumIdentificationList22 = dc.getAnalysisData().getSpectrumIdentificationList();

//			}
            for (int j = 0; j < spectrumIdentificationList22.size(); j++) {
                SpectrumIdentificationList sil = spectrumIdentificationList22.get(j);
                List<SpectrumIdentificationResult> sirs = sil.getSpectrumIdentificationResult();

                for (SpectrumIdentificationResult spectrumIdentificationResult : sirs) {
                    spectrumIdentificationResult.setSpectraData(id_spectraData.get("com_" + i + "_" + spectrumIdentificationResult.getSpectraDataRef()));
                    spectrumIdentificationResult.setId("com_" + i + "_" + spectrumIdentificationResult.getId());
                    spectrumIdentificationResult.setSpectrumID("com_" + i + "_" + spectrumIdentificationResult.getSpectrumID());
                    List<SpectrumIdentificationItem> sii = spectrumIdentificationResult.getSpectrumIdentificationItem();
                    List<SpectrumIdentificationItem> sii_new = new ArrayList<SpectrumIdentificationItem>();
                    for (SpectrumIdentificationItem spectrumIdentificationItem : sii) {
                        Boolean pass_item = false;
                        if (spectrumIdentificationItem.getRank() > 1) {
                            pass_item = true;
                        }
                        List<CvParam> cvParams = spectrumIdentificationItem.getCvParam();
                        for (CvParam cvParam : cvParams) {
                            if (cvParam.getAccession().equals(threshold_score_acc)) {
                                if (Double.valueOf(cvParam.getValue()) > threshold_score) {
                                    pass_item = true;
                                }
                            }
                        }
                        spectrumIdentificationItem.setId("com_" + i + "_" + spectrumIdentificationItem.getId());
                        spectrumIdentificationItem.setPeptide(pep_id_replace.get("com_" + i + "_" + spectrumIdentificationItem.getPeptideRef()));
                        List<PeptideEvidenceRef> peptiideER = spectrumIdentificationItem.getPeptideEvidenceRef();
                        for (PeptideEvidenceRef peptideEvidenceRef : peptiideER) {
                            String com_refID = peptideEvidenceRef.getPeptideEvidenceRef();
                            peptideEvidenceRef.setPeptideEvidence(peptideEvidence_id_replace.get("com_" + i + "_" + com_refID));
                        }
                        if (!pass_item) {
                            sii_new.add(spectrumIdentificationItem);
                        }
                    }
                    spectrumIdentificationResult.getSpectrumIdentificationItem().clear();
                    spectrumIdentificationResult.getSpectrumIdentificationItem().addAll(sii_new);
                }

                combinedSpectrumIdentificationList.add(sil);
            }
            dataCollections.add(dc);
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
            System.exit(1);
        }
    }

    private void writeMzidFile(String combinemzidFileName) {
        try {
            String outFile = combinemzidFileName;
            FileWriter writer = new FileWriter(outFile);
            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            writer.write(marshaller.createXmlHeader() + "\n");
            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");

            CvList combineCvList = new CvList();
            HashMap<String, Cv> cvHashMap = new HashMap<String, Cv>();
            for (int i = 0; i < cvList.size(); i++) {
                CvList cvList1 = cvList.get(i);
                for (int j = 0; j < cvList1.getCv().size(); j++) {
                    cvHashMap.put(cvList1.getCv().get(j).getId(), cvList1.getCv().get(j));
                }
            }
            combineCvList.getCv().addAll(cvHashMap.values());
            if (combineCvList != null) {
                marshaller.marshal(combineCvList, writer);
            }
            writer.write("\n");

            AnalysisSoftwareList combineanalysisSoftwareList = new AnalysisSoftwareList();
            HashMap<String, AnalysisSoftware> analysisSoftwareHashMap = new HashMap<String, AnalysisSoftware>();
            for (int i = 0; i < analysisSoftwareList.size(); i++) {
                AnalysisSoftwareList analysisSoftwareList1 = analysisSoftwareList.get(i);
                for (int j = 0; j < analysisSoftwareList1.getAnalysisSoftware().size(); j++) {
                    analysisSoftwareHashMap.put(analysisSoftwareList1.getAnalysisSoftware().get(j).getId(), analysisSoftwareList1.getAnalysisSoftware().get(j));
                }
            }
            combineanalysisSoftwareList.getAnalysisSoftware().addAll(analysisSoftwareHashMap.values());
            if (combineanalysisSoftwareList != null) {
                marshaller.marshal(combineanalysisSoftwareList, writer);
            }
            writer.write("\n");

            if (provider.get(0) != null) {
                marshaller.marshal(provider.get(0), writer);
            }
            writer.write("\n");

            if (auditCollection.get(0) != null) {
                marshaller.marshal(auditCollection.get(0), writer);
            }
            writer.write("\n");

            SequenceCollection sequenceCollection = new SequenceCollection();
            sequenceCollection.getDBSequence().addAll(dBSequenceHashMap.values());
            sequenceCollection.getPeptide().addAll(peptideHashMap.values());
            sequenceCollection.getPeptideEvidence().addAll(peptideEvidenceHashMap.values());

            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");

            AnalysisCollection combinedanalysisCollection = new AnalysisCollection();
            HashMap<String, SpectrumIdentification> spectrumIdentificationHashMap = new HashMap<String, SpectrumIdentification>();
            for (int i = 0; i < analysisCollection.size(); i++) {
                AnalysisCollection analysisCollection1 = analysisCollection.get(i);
                for (int j = 0; j < analysisCollection1.getSpectrumIdentification().size(); j++) {
                    spectrumIdentificationHashMap.put(analysisCollection1.getSpectrumIdentification().get(j).getId(), analysisCollection1.getSpectrumIdentification().get(j));
                }
            }
            combinedanalysisCollection.getSpectrumIdentification().addAll(spectrumIdentificationHashMap.values());
            if (combinedanalysisCollection != null) {
                marshaller.marshal(combinedanalysisCollection, writer);
            }
            writer.write("\n");

            AnalysisProtocolCollection combinedanalysisProtocolCollection = new AnalysisProtocolCollection();
            HashMap<String, SpectrumIdentificationProtocol> spectrumIdentificationProtocolHashMap = new HashMap<String, SpectrumIdentificationProtocol>();
            for (int i = 0; i < analysisProtocolCollection.size(); i++) {
                AnalysisProtocolCollection analysisProtocolCollection1 = analysisProtocolCollection.get(i);
                for (int j = 0; j < analysisProtocolCollection1.getSpectrumIdentificationProtocol().size(); j++) {
                    spectrumIdentificationProtocolHashMap.put(analysisProtocolCollection1.getSpectrumIdentificationProtocol().get(j).getId(), analysisProtocolCollection1.getSpectrumIdentificationProtocol().get(j));
                }
            }
            combinedanalysisProtocolCollection.getSpectrumIdentificationProtocol().addAll(spectrumIdentificationProtocolHashMap.values());
            if (combinedanalysisProtocolCollection != null) {
                marshaller.marshal(combinedanalysisProtocolCollection, writer);
            }
            writer.write("\n");

            DataCollection combined_dc = dataCollections.get(0);
            for (SpectraData spectraData : id_spectraData.values()) {
                if (combined_dc.getInputs().getSpectraData().contains(spectraData)) {
                    continue;
                }
                combined_dc.getInputs().getSpectraData().add(spectraData);
            }

            if (!name_id.values().contains(1)) {
                SpectraData defaultSpectraData = new SpectraData();
                defaultSpectraData.setId("SID_1");
                defaultSpectraData.setLocation("dir");
                combined_dc.getInputs().getSpectraData().add(defaultSpectraData);
            }

            SpectrumIdentificationList list = new SpectrumIdentificationList();
            list.setId(combinedSpectrumIdentificationList.get(0).getId());
            for (int i = 0; i < combinedSpectrumIdentificationList.size(); i++) {
                List<SpectrumIdentificationResult> spectrumIdentificationResult = combinedSpectrumIdentificationList.get(i).getSpectrumIdentificationResult();
                for (SpectrumIdentificationResult spectrumIdentificationResult2 : spectrumIdentificationResult) {
                    if (spectrumIdentificationResult2.getSpectrumIdentificationItem().isEmpty()) {
                        continue;
                    }
                    list.getSpectrumIdentificationResult().add(spectrumIdentificationResult2);
                }
            }
            combined_dc.getAnalysisData().getSpectrumIdentificationList().clear();
            combined_dc.getAnalysisData().getSpectrumIdentificationList().add(list);
            marshaller.marshal(combined_dc, writer);
            writer.write("\n");

            writer.write(marshaller.createMzIdentMLClosingTag());
            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
            System.exit(1);
        }
    }
}
