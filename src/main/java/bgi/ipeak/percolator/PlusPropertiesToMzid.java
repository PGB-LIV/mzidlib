package bgi.ipeak.percolator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Cv;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class PlusPropertiesToMzid {

    @Option(name = "-mzid", required = true, usage = "(required) Path and file name of input mzIdentML file.")
    private String mzid;
    @Option(name = "-plus", required = true, usage = "(required) Path and file name of input plus file.")
    private String featureFile;
    @Option(name = "-out", required = false, usage = "(optional) Path and file name of output mzIdentML file.")
    private String outFile;
    @Option(name = "-user", required = false, usage = "User-defined plus.If true, only add UserParam to mzid. Conflict with -replace.")
    private boolean user = false;
//	private Double threshold_score=0.2;

    /**
     * Constructors for AddFeaturesToMzid.
     *
     * @param mzid mzIdentML file.
     * @param per_file Percolator file.
     * @param outfile Output file.
     * @throws IOException Exception thrown by the constructor.
     */
    public PlusPropertiesToMzid(String mzid, String per_file, String outfile) throws IOException {
        this.mzid = mzid;
        this.featureFile = per_file;
        this.outFile = outfile;
    }

    public PlusPropertiesToMzid(String mzid, String per_file) throws IOException {
        this.mzid = mzid;
        this.featureFile = per_file;
        String addFeatures = mzid.substring(0, mzid.lastIndexOf(".")) + "AddP.mzid";
        this.outFile = addFeatures;
    }

    public PlusPropertiesToMzid() {

    }

    /**
     * Add features to mzid.
     * @param args Arguments to the main method.
     * @throws IOException IOException thrown by main method.
     */
    public static void main(String[] args) throws IOException {
        PlusPropertiesToMzid ad = new PlusPropertiesToMzid();
        CmdLineParser parser = new CmdLineParser(ad);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp iPeak.jar bgi.ipeak.PlusPropertiesToMzid [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        ad.export();
    }

    public void export() throws IOException {
        System.out.println("Add percolator scores " + featureFile + " to " + this.mzid);
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(this.mzid));

        //These are the main structures to be output by the main writing method
        SequenceCollection sequenceCollection = unmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
        CvList cvList = unmarshaller.unmarshal(MzIdentMLElement.CvList);
        AnalysisSoftwareList analysisSoftwareList = unmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
        Provider provider = unmarshaller.unmarshal(MzIdentMLElement.Provider);
        AuditCollection auditCollection = unmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
        //AnalysisSampleCollection analysisSampleCollection=unmarshaller.unmarshal(MzIdentMLElement.AnalysisSampleCollection);
        AnalysisCollection analysisCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
        AnalysisProtocolCollection analysisProtocolCollection = unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
        DataCollection dataCollection = unmarshaller.unmarshal(MzIdentMLElement.DataCollection);
        if (dataCollection == null) {
            System.err.println("Failed! DataCollection is null!");
            System.exit(0);
        }

        BufferedReader ff = new BufferedReader(new FileReader(new File(this.featureFile)));
        if (this.outFile == null) {
            this.outFile = this.mzid.substring(0, this.mzid.lastIndexOf(".")) + "_plus.mzid";
        }
        BufferedWriter wf = new BufferedWriter(new FileWriter(new File(this.outFile)));

        String fline = ff.readLine();
        if (fline == null) {
            wf.close();
            ff.close();
            System.err.println("Error with empty feature file: " + this.featureFile);
            return;
        }
        String[] head = fline.split("\t");
        HashMap<String, String[]> featuresMap = new HashMap<String, String[]>();
        while ((fline = ff.readLine()) != null) {
            String[] ltmp = fline.split("\t");
            //String index=ltmp[0].split(":")[1]+"."+ltmp[1];
            String index = ltmp[0].split(":")[1];
            featuresMap.put(index, ltmp);
        }
        HashMap<String, PeptideEvidence> pepEviMap = new HashMap<String, PeptideEvidence>();
        for (PeptideEvidence pepEvi : sequenceCollection.getPeptideEvidence()) {
            pepEviMap.put(pepEvi.getId(), pepEvi);
        }
        List<SpectrumIdentificationList> sil = dataCollection.getAnalysisData().getSpectrumIdentificationList();
//		List <SpectrumIdentificationList> threshold_sil= new ArrayList<SpectrumIdentificationList>();
        for (SpectrumIdentificationList sIdentList : sil) {
            for (SpectrumIdentificationResult spectrumIdentResult
                    : sIdentList.getSpectrumIdentificationResult()) {

                // Get the name of SpectrumIdentificationResult
                String spectrumID = spectrumIdentResult.getSpectrumID();
                String index = spectrumID.split("=")[1];
                for (SpectrumIdentificationItem spectrumIdentItem : spectrumIdentResult.getSpectrumIdentificationItem()) {
                    index += "." + spectrumIdentItem.getRank();
                    if (featuresMap.containsKey(index)) {
                        if (this.user) {
                            for (int i = 2; i < head.length - 2; i++) {
                                UserParam up = new UserParam();
                                up.setName(head[i]);
                                up.setValue(featuresMap.get(index)[i]);
                                spectrumIdentItem.getUserParam().add(up);
                            }
                        } else {
                            Cv psiCV;
                            psiCV = new Cv();
                            psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
                            psiCV.setId("PSI-MS");
                            psiCV.setVersion("2.25.0");
                            psiCV.setFullName("PSI-MS");
                            //replace other score with percolator score.

                            Boolean has_score = false;
                            Boolean has_pep = false;
                            Boolean has_qvalue = false;
                            for (CvParam cvp : spectrumIdentItem.getCvParam()) {
                                if (cvp.getAccession().equals("MS:1001492")) {
                                    cvp.setValue(featuresMap.get(index)[2]);
                                    has_score = true;
                                } else if (cvp.getAccession().equals("MS:1001491")) {
                                    cvp.setValue(featuresMap.get(index)[3]);
                                    has_pep = true;
                                } else if (cvp.getAccession().equals("MS:1001493")) {
                                    cvp.setValue(featuresMap.get(index)[4]);
                                    has_qvalue = true;
                                }
                            }
                            if (!has_score) {

                                CvParam cp = new CvParam();
                                cp.setAccession("MS:1001492");
                                cp.setName("percolator:score");
                                cp.setValue(featuresMap.get(index)[2]);
                                cp.setCv(psiCV);
                                spectrumIdentItem.getCvParam().add(cp);
                            }
                            if (!has_qvalue) {
                                //add percolator qvalue.
                                CvParam cpq = new CvParam();
                                cpq.setAccession("MS:1001491");
                                cpq.setName("percolator:Q value");
                                cpq.setValue(featuresMap.get(index)[3]);
                                cpq.setCv(psiCV);
                                spectrumIdentItem.getCvParam().add(cpq);
                            }
                            if (!has_pep) {
                                //add percolator PEP.
                                CvParam cpp = new CvParam();
                                cpp.setAccession("MS:1001493");
                                cpp.setName("percolator:PEP");
                                cpp.setValue(featuresMap.get(index)[4]);
                                cpp.setCv(psiCV);
                                spectrumIdentItem.getCvParam().add(cpp);
                            }
                            if (this.user) {
                                for (int i = 2; i < head.length - 2; i++) {
                                    UserParam up = new UserParam();
                                    up.setName(head[i]);
                                    up.setValue("NaN");
                                    spectrumIdentItem.getUserParam().add(up);
                                }
                            }
                        }
                    } else {
                        Cv psiCV;
                        psiCV = new Cv();
                        psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
                        psiCV.setId("PSI-MS");
                        psiCV.setVersion("2.25.0");
                        psiCV.setFullName("PSI-MS");
                        //replace other score with percolator score.
                        Boolean has_score = false;
                        Boolean has_pep = false;
                        Boolean has_qvalue = false;
                        for (CvParam cvp : spectrumIdentItem.getCvParam()) {
                            if (cvp.getAccession().equals("MS:1001492")) {
                                cvp.setValue("-100");
                                has_score = true;
                            } else if (cvp.getAccession().equals("MS:1001491")) {
                                cvp.setValue("1");
                                has_pep = true;
                            } else if (cvp.getAccession().equals("MS:1001493")) {
                                cvp.setValue("1");
                                has_qvalue = true;
                            }
                        }
                        if (!has_score) {

                            CvParam cp = new CvParam();
                            cp.setAccession("MS:1001492");
                            cp.setName("percolator:score");
                            cp.setValue("-100");
                            cp.setCv(psiCV);
                            spectrumIdentItem.getCvParam().add(cp);
                        }
                        if (!has_qvalue) {
                            //add percolator qvalue.
                            CvParam cpq = new CvParam();
                            cpq.setAccession("MS:1001491");
                            cpq.setName("percolator:Q value");
                            cpq.setValue("1");
                            cpq.setCv(psiCV);
                            spectrumIdentItem.getCvParam().add(cpq);
                        }
                        if (!has_pep) {
                            //add percolator PEP.
                            CvParam cpp = new CvParam();
                            cpp.setAccession("MS:1001493");
                            cpp.setName("percolator:PEP");
                            cpp.setValue("1");
                            cpp.setCv(psiCV);
                            spectrumIdentItem.getCvParam().add(cpp);
                        }
                        if (this.user) {
                            for (int i = 2; i < head.length - 2; i++) {
                                UserParam up = new UserParam();
                                up.setName(head[i]);
                                up.setValue("NaN");
                                spectrumIdentItem.getUserParam().add(up);
                            }
                        }
                    }
                }
            }
            ff.close();
            MzIdentMLMarshaller m = new MzIdentMLMarshaller();
            wf.write(m.createXmlHeader() + "\n");

            // mzIdentML start tag
            wf.write(m.createMzIdentMLStartTag("12345") + "\n");
            if (cvList != null) {
                wf.write(m.marshal(cvList) + "\n");
            }
            if (analysisSoftwareList != null) {
                wf.write(m.marshal(analysisSoftwareList) + "\n");
            }
            if (provider != null) {
                wf.write(m.marshal(provider) + "\n");
            }
            if (auditCollection != null) {
                wf.write(m.marshal(auditCollection) + "\n");
            }
            if (sequenceCollection != null) {
                wf.write(m.marshal(sequenceCollection) + "\n");
            }
            if (analysisCollection != null) {
                wf.write(m.marshal(analysisCollection) + "\n");
            }
            if (analysisProtocolCollection != null) {
                wf.write(m.marshal(analysisProtocolCollection) + "\n");
            }
            if (dataCollection != null) {
                wf.write(m.marshal(dataCollection) + "\n");
            }
            wf.write(m.createMzIdentMLClosingTag());
            wf.close();
        }
    }
}
