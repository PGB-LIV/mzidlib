
package uk.ac.liv.mzidlib.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
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

import javax.xml.bind.JAXBException;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
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
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;

import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.gff.CdsInformation;
import uk.ac.liv.mzidlib.gff.ProteinResults;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 *
 * @author Fawaz Ghali
 */
public class AddGenomeCoordinatesForPeptides {

    private String inputMzid;
    private String outputMzid;
    private String inputGff;
    private String outputGff;
    private Cv psiCV;
    private MzIdentMLVersion version;

    private final List unmappedAccessions = new ArrayList();

    private final Map<String, List<CdsInformation>> cdsRecords = new HashMap<>();

    public AddGenomeCoordinatesForPeptides(String inputMzid, String outputMzid,
                                           String inputGff, String outputGff,
                                           String ver) {
        this.inputMzid = inputMzid;
        this.outputMzid = outputMzid;
        this.outputGff = outputGff;
        this.inputGff = inputGff;
        this.version = MzIdentMLVersion.getVersion(ver);

    }

    /**
     *
     * @param inputGffFile
     * @param outputGff
     * @param proteinHits
     */
    private void writingTheGFFfile(String inputGffFile, String outputGffFile,
                                   Map<String, ProteinResults> proteinHits) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(inputGffFile));
            String line;
            boolean fastaRegionFlag = false;

            Writer fstream = new FileWriter(outputGffFile);
            Writer out = new BufferedWriter(fstream);

            while ((line = in.readLine()) != null) {

                if (line.startsWith("##")) {
                    if (line.equals("##FASTA")) {
                        fastaRegionFlag = true;
                    }
                }

                if (fastaRegionFlag == false) {

                    String[] lineArray = line.split("\t");
                    // Added by Fawaz Ghali 22/08/2016 to handle ab inito gff3
                    if (inputGffFile.contains("abinitio") && lineArray.length
                            == 9 && lineArray[2].equals("exon")) {
                        out.write(line + "\n");

                        String accession = "";
                        String[] gffArrayLine = line.split("\t");
                        String seqId = gffArrayLine[0];
                        String source = gffArrayLine[1];
                        String type = gffArrayLine[2];
                        String start = gffArrayLine[3];
                        String end = gffArrayLine[4];
                        long startPos = Long.parseLong(start);
                        long endPos = Long.parseLong(end);
                        String score = gffArrayLine[5];
                        String strand = gffArrayLine[6];
                        String phase = gffArrayLine[7];
                        String attribute = gffArrayLine[8];

                        String[] splits = attribute.split(";");
                        for (int i = 0; i < splits.length; i++) {
                            String string = splits[i];
                            if (string.contains("Parent=transcript:")) {
                                accession = string.replace("Parent=transcript:",
                                                           "");
                            }
                        }
                        accession = accession.toLowerCase();
                        CdsInformation cdsObj = new CdsInformation(seqId,
                                                                   source,
                                                                   startPos,
                                                                   endPos,
                                                                   strand,
                                                                   phase,
                                                                   attribute);
                        List<CdsInformation> cdsColl = cdsRecords.get(
                                accession);

                        if (cdsColl == null) {
                            cdsColl = new ArrayList<>();
                        }
                        cdsColl.add(cdsObj);
                        cdsRecords.put(accession, cdsColl);
                    }

                    if (lineArray.length == 9 && lineArray[2].equals("CDS")) {
                        out.write(line + "\n");

                        String[] gffArrayLine = line.split("\t");

                        String attribute = gffArrayLine[8];

                        String accession = null;
                        accession = attribute.split(";")[0];
                        if (accession.startsWith("ID=apidb|")) {
                            accession = accession.substring(9);
                        }
                        if (accession.startsWith("transcript_id")) {
                            accession = accession.substring(13);
                            accession = accession.trim();
                            accession = accession.replaceAll("\"", "");
                        }
                        if (accession.startsWith("ID=cds_")) {
                            String[] splits = attribute.split(";");
                            accession = splits[splits.length - 1];
                            accession = accession.replaceAll("Parent=", "");
                        }
                        if (accession.startsWith("ID=gb|")) {
                            accession = accession.substring(6);

                        }
                        if (accession.startsWith("genscan_id=")) {
                            accession = accession.substring(11);
                        }
                        if (accession.startsWith("ID=CDS:")) {
                            String[] splits = attribute.split(";");
                            for (int i = 0; i < splits.length; i++) {
                                String string = splits[i];
                                if (string.startsWith("protein_id=")) {
                                    String[] proteinId = string.split("=");
                                    accession = proteinId[1];
                                }
                            }
                        }

                        if (accession.startsWith("ID=")) {
                            accession = accession.substring(3);
                        }
                        accession = accession.toLowerCase();
                        List<CdsInformation> cdsColl = cdsRecords.get(
                                accession);
                        String seqId = gffArrayLine[0];
                        String source = gffArrayLine[1];
                        String type = gffArrayLine[2];
                        String start = gffArrayLine[3];
                        String end = gffArrayLine[4];
                        String score = gffArrayLine[5];
                        String strand = gffArrayLine[6];
                        String phase = gffArrayLine[7];
                        long startPos = Long.parseLong(start);
                        long endPos = Long.parseLong(end);
                        CdsInformation cdsObj = new CdsInformation(seqId,
                                                                   source,
                                                                   startPos,
                                                                   endPos,
                                                                   strand,
                                                                   phase,
                                                                   attribute);
                        if (cdsColl == null) {
                            cdsColl = new ArrayList<>();
                        }
                        cdsColl.add(cdsObj);
                        cdsRecords.put(accession, cdsColl);
                    }

                }

                // Come out when ##FASTA encountered
                if (fastaRegionFlag) {
                    break;
                }

            } // end of while

            // Exit if no CDS found
            if (cdsRecords.isEmpty()) {
                String exMessage = "No CDS field found.....exiting";
                System.out.println(exMessage);
                in.close();
                out.close();
                throw new Exception(exMessage);
            }

            in.close();
            // Perform Chromosome mapping
            // for (int i = 0; i < proteinHits.size(); i++) {

            for (Map.Entry<String, ProteinResults> pair : proteinHits.entrySet()) {
                ProteinResults pr = pair.getValue();
                String[] coords = new String[2];
                long startMap;
                long endMap;

                coords = mapToGff(pr);
                if (coords != null) {
                    // Take care of the negative strand reporting in GFF
                    startMap = Long.parseLong(coords[0]);
                    endMap = Long.parseLong(coords[1]);
                    if (endMap < startMap) {
                        long tmp = startMap;
                        startMap = endMap;
                        endMap = tmp;
                    }

                    String accession = pr.getAccession();
                    accession = accession.toLowerCase();
                    List<CdsInformation> cdsColl = cdsRecords.get(accession);

                    String feature = "peptide";
                    String score = ".";
                    String mappedStartLocation = Long.toString(startMap);
                    String mappedEndLocation = Long.toString(endMap);

                    String attr = "ID=pep_" + accession + "_"
                            + ";description=peptide_seq:" + pr.getPeptideSeq()
                            + "_fdr_score:" + pr.getfdrScore()
                            + "_peptide_start: " + pr.getstart()
                            + "_peptide_end: " + pr.getEnd() + ";Derives_from="
                            + accession;
                    // Just use the information from the first CDS entry 
                    //with mapped start and end locations
                    String gffEntry = cdsColl.get(0).getSeqId() + "\t"
                            + cdsColl.get(0).getSource() + "\t" + feature + "\t"
                            + mappedStartLocation
                            + "\t" + mappedEndLocation + "\t" + score + "\t"
                            + cdsColl.get(0).getStrand() + "\t"
                            + cdsColl.get(0).getSePhase() + "\t" + attr + "\n";
                    out.write(gffEntry);
                }
            }
            out.close();

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void writeMappingResults() {
        Map<String, ProteinResults> prList = getProteinResults(inputMzid);
        writingTheGFFfile(inputGff, outputGff, prList);

        try {
            String outFile = outputMzid;

            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller(version);

            //cvList = mzIdentML.getCvList();
            MzIdentMLUnmarshaller mzIdentMlUnmarshaller
                    = new MzIdentMLUnmarshaller(new File(inputMzid));

            CvList cvList = mzIdentMlUnmarshaller.unmarshal(
                    MzIdentMLElement.CvList);
            Iterator<Cv> iterCv = mzIdentMlUnmarshaller
                    .unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                }
            }
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            AnalysisSoftwareList analysisSoftwareList = mzIdentMlUnmarshaller.
                    unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            AuditCollection auditCollection = mzIdentMlUnmarshaller.unmarshal(
                    MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            Provider provider = mzIdentMlUnmarshaller.unmarshal(
                    MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            AnalysisProtocolCollection analysisProtocolCollection
                    = mzIdentMlUnmarshaller.unmarshal(
                            MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            AnalysisCollection analysisCollection = mzIdentMlUnmarshaller.
                    unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            Inputs inputs = mzIdentMlUnmarshaller.unmarshal(
                    MzIdentMLElement.Inputs);

            writer.write(marshaller.createXmlHeader() + "\n");

            String mzId = mzIdentMlUnmarshaller.getMzIdentMLId();
            if (mzId != null) {
                writer.write(marshaller.createMzIdentMLStartTag(mzId) + "\n");
            } else {
                writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            }

            if (cvList != null) {
                marshaller.marshal(cvList, writer);
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

            param.setParam(CvConstants.MZIDLIB);
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

            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection
                        .getSpectrumIdentification().get(0)
                        .getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable
                    = mzIdentMlUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Iterator<SpectrumIdentificationResult> sirIter
                    = mzIdentMlUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (sirIter.hasNext()) {

                SpectrumIdentificationResult sr = sirIter.next();
                siList.getSpectrumIdentificationResult().add(sr);

            }

            String proteinDetectionListRef = "";
            if (analysisCollection.getProteinDetection() != null) {
                proteinDetectionListRef = analysisCollection
                        .getProteinDetection().getProteinDetectionListRef();
            }
            ProteinDetectionList pdList;
            pdList = new ProteinDetectionList();

            pdList.setId(proteinDetectionListRef);

            Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup
                    = mzIdentMlUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {
                ProteinAmbiguityGroup pag = iterProteinAmbiguityGroup.next();
                pdList.getProteinAmbiguityGroup().add(pag);
            }

            // here 
            SequenceCollection sequenceCollection1 = mzIdentMlUnmarshaller
                    .unmarshal(MzIdentMLElement.SequenceCollection);

            SequenceCollection sq = new SequenceCollection();
            List<DBSequence> dbSeqList = sequenceCollection1.getDBSequence();
            HashMap dbSeqHm = new HashMap();
            for (int i = 0; i < dbSeqList.size(); i++) {
                DBSequence get = dbSeqList.get(i);
                dbSeqHm.put(get.getId(), get);

            }
            sq.getPeptide().addAll(sequenceCollection1.getPeptide());

            List<PeptideEvidence> peList = sequenceCollection1
                    .getPeptideEvidence();
            List<PeptideEvidence> newPeList = new ArrayList<>();
            for (int i = 0; i < peList.size(); i++) {
                PeptideEvidence peptideEvidence = peList.get(i);
                DBSequence dbSeq = mzIdentMlUnmarshaller.unmarshal(
                        DBSequence.class, peptideEvidence.getDBSequenceRef());
                String accession = dbSeq.getAccession();
                if (accession.startsWith("generic")) {
                    accession = accession.split("\\|")[1];
                    accession = accession.substring(2);
                    // Update by Fawaz Ghali 04/08/2016 fix Ensembl mapping version 85
                    if (accession.contains(".")) {
                        accession = accession.substring(0, accession
                                                        .indexOf("."));
                    }
                } else {
                    throw new RuntimeException("Not generic accession");
                }
                accession = accession.toLowerCase();
                List<CdsInformation> gffData = cdsRecords.get(accession);
                if (gffData != null && gffData.size() > 0 && !peptideEvidence
                        .isIsDecoy()) {
                    ProteinResults pr = prList.get(peptideEvidence.getId());
                    String[] coords = new String[2];
                    long startMap;
                    long endMap;

                    coords = mapToGff(pr);

                    // unmapped peptide
                    if (coords == null || unmappedAccessions.contains(pr
                            .getAccession())) {
                        if (!dbSeq.getCvParam().contains(
                                CvConstants.UNMAPPED_PROTEIN)) {
                            dbSeq.getCvParam().add(MzidLibUtils.makeCvParam(
                                    "MS:1002741",
                                    "unmapped protein",
                                    gffData.get(0).getSeqId(), psiCV));
                        }
                        if (!peptideEvidence.getCvParam().contains(
                                CvConstants.UNMAPPED_PEPTIDE)) {
                            peptideEvidence.getCvParam().add(MzidLibUtils
                                    .makeCvParam(
                                            "MS:1002740", "unmapped peptide",
                                            gffData.get(0).getSeqId(), psiCV));
                        }
                    } else {
                        startMap = Long.parseLong(coords[0]);
                        endMap = Long.parseLong(coords[1]);
                        if (endMap < startMap) {
                            long tmp = startMap;
                            startMap = endMap;
                            endMap = tmp;
                        }
                        // Added by Fawaz Ghali 13/8/2015
                        // Update the code to cover the case where a peptide covers multiple exons
                        List<CdsInformation> sortedCds
                                = sortCdsAccordingToStartPosition(gffData);
                        List<CdsInformation> outputCds = new ArrayList();
                        int countCds = 0;
                        boolean sorted = false;
                        if (gffData.get(0).getStrand().equals("-")) {
                            sortedCds = sortCdsAccordingToStartPositionTemp(
                                    sortedCds);
                        }
                        for (int j = 0; j < sortedCds.size(); j++) {
                            CdsInformation object = sortedCds.get(j);
                            long s = object.getStart();
                            long e = object.getEnd();
                            long newStart, newEnd;
                            countCds = j;
                            if (startMap > e) {
                                continue;
                            } else {
                                newStart = startMap;
                                if (endMap <= e) {
                                    newEnd = endMap;
                                    sorted = true;
                                } else {
                                    newEnd = e;
                                }

                                CdsInformation temp = new CdsInformation("",
                                                                         "",
                                                                         newStart,
                                                                         newEnd,
                                                                         "",
                                                                         "",
                                                                         "");
                                outputCds.add(temp);
                                break;
                            }
                        }
                        if (!sorted) {
                            for (int j = countCds + 1; j < sortedCds.size(); j++) {
                                CdsInformation object = sortedCds.get(j);
                                long s = object.getStart();
                                long e = object.getEnd();
                                long newStart, newEnd;

                                if (endMap > e) {
                                    CdsInformation temp = new CdsInformation(
                                            "", "", s, e, "", "", "");
                                    outputCds.add(temp);
                                } else {
                                    newStart = s;
                                    newEnd = endMap;

                                    CdsInformation temp = new CdsInformation(
                                            "", "", newStart, newEnd, "", "", "");
                                    outputCds.add(temp);
                                    break;
                                }
                            }
                        }
                        int pepLen = 3 * (peptideEvidence.getEnd()
                                - peptideEvidence.getStart() + 1);
                        int coordLen = 0;

                        long s, e = 0;
                        String startsList = "";
                        String positionsList = "";

                        for (int j = 0; j < outputCds.size(); j++) {
                            CdsInformation cdsInfo = outputCds.get(j);
                            s = cdsInfo.getStart() - 1;

                            startsList = startsList + s + ",";
                            if (j == 0) {
                                e = cdsInfo.getEnd() + 2;
                            } else {
                                e = cdsInfo.getEnd();
                            }

                            coordLen = coordLen + MzidLibUtils.safeLongToInt(e
                                    - s);
                            positionsList = positionsList + String
                                    .valueOf(e - s) + ",";

                        }

                        dbSeq.getCvParam().add(MzidLibUtils.makeCvParam(
                                "MS:1002637",
                                "chromosome name",
                                gffData.get(0).getSeqId(), psiCV));

                        dbSeq.getCvParam().add(MzidLibUtils.makeCvParam(
                                "MS:1002638",
                                "chromosome strand",
                                gffData.get(0).getStrand(), psiCV));
                        dbSeq.getCvParam().add(MzidLibUtils.makeCvParam(
                                "MS:1002644",
                                "genome reference version",
                                new File(inputGff).getName(), psiCV));
                        dbSeqHm.put(dbSeq.getId(), dbSeq);

                        peptideEvidence.getCvParam().add(MzidLibUtils
                                .makeCvParam(
                                        "MS:1002640",
                                        "peptide end on chromosome",
                                        String.valueOf(endMap), psiCV));

                        peptideEvidence.getCvParam().add(MzidLibUtils
                                .makeCvParam(
                                        "MS:1002641", "peptide exon count",
                                        String.valueOf(outputCds.size()), psiCV));
                        if (startsList.endsWith(",")) {
                            startsList = startsList.substring(0, startsList
                                                              .length() - 1);
                        }
                        if (positionsList.endsWith(",")) {
                            positionsList = positionsList.substring(0,
                                                                    positionsList
                                                                    .length()
                                                                    - 1);
                        }
                        peptideEvidence.getCvParam().add(MzidLibUtils
                                .makeCvParam(
                                        "MS:1002642",
                                        "peptide exon nucleotide sizes",
                                        positionsList, psiCV));
                        peptideEvidence.getCvParam().add(MzidLibUtils
                                .makeCvParam(
                                        "MS:1002643",
                                        "peptide start positions on chromosome",
                                        startsList, psiCV));

                    }
                } else if (gffData == null && !peptideEvidence.isIsDecoy()) {
                    dbSeq.getCvParam().add(MzidLibUtils
                            .makeCvParam("MS:1002741",
                                         "unmapped protein", psiCV));
                    peptideEvidence.getCvParam().add(MzidLibUtils.makeCvParam(
                            "MS:1002740",
                            "unmapped peptide",
                            psiCV));
                    dbSeqHm.put(dbSeq.getId(), dbSeq);
                }
                newPeList.add(peptideEvidence);
            }
            sq.getDBSequence().addAll(dbSeqHm.values());
            sq.getPeptideEvidence().addAll(newPeList);

            marshaller.marshal(sq, writer);

            writer.write("\n");

            marshaller.marshal(analysisCollection, writer);

            writer.write("\n");

            if (analysisProtocolCollection != null) {

                List<SpectrumIdentificationProtocol> sipList
                        = analysisProtocolCollection
                        .getSpectrumIdentificationProtocol();
                for (SpectrumIdentificationProtocol sip : sipList) {
                    List<CvParam> cvs = sip.getAdditionalSearchParams()
                            .getCvParam();
                    cvs.remove(CvConstants.NO_SPECIAL_PROCESSING);
                    if (!cvs.contains(CvConstants.PROTEOGENOMICS_SEARCH)) {
                        cvs.add(CvConstants.PROTEOGENOMICS_SEARCH);
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

            marshaller.marshal(siList, writer);
            writer.write("\n");

            if (!pdList.getProteinAmbiguityGroup().isEmpty()) {
                marshaller.marshal(pdList, writer);
                writer.write("\n");
            }

            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");
            writer.write(marshaller.createMzIdentMLClosingTag());
            writer.close();
            System.out.println("Output written to " + outFile);
        } catch (JAXBException | IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1]
                    .getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + e
                    .getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information "
                    + "on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    public String[] mapToGff(ProteinResults pr) {
        String[] gffMapping = new String[3];

        // Get accession from the protein object
        String accession = pr.getAccession();
        accession = accession.toLowerCase();
        // return if no key found
        if (!cdsRecords.containsKey(accession)) {
            return null;
        }

        //...find the cds this range falls into...
        List<CdsInformation> cdsColl = cdsRecords.get(accession);

        long[] mappedCords = determineTheLocationOfSeqOnCds(pr, cdsColl);

        gffMapping[0] = Long.toString(mappedCords[0]);
        gffMapping[1] = Long.toString(mappedCords[1]);
        gffMapping[2] = cdsRecords.get(accession).get(0).getSeqId();

        return gffMapping;
    }

    /**
     * Determine the location of sequence on CDS.
     *
     * @param pr      protein results
     * @param cdsColl CDS column list
     *
     * @return an array of locations
     */
    long[] determineTheLocationOfSeqOnCds(ProteinResults pr,
                                          List<CdsInformation> cdsColl) {
        long[] gffEntry = new long[2];

        // Get locations from the protein object
        long start = pr.getstart();
        long end = pr.getEnd();

        // sort the CDS collection according to the strand
        List<CdsInformation> sortedCds
                = sortCdsAccordingToStartPosition(cdsColl);

        long mappedStart = getMappedCordinates(start, sortedCds, pr);
        long mappedEnd = getMappedCordinates(end, sortedCds, pr);
        if (mappedStart == -1) {
            System.out.println("For accession: " + pr.getAccession()
                    + " and a peptide eveidence: " + pr.getPeptideEvidenceId()
                    + " start coordinates couldn't be mapped.");
            unmappedAccessions.add(pr.getAccession());
        }
        if (mappedEnd == -1) {
            System.out.println("For accession: " + pr.getAccession()
                    + " and a peptide eveidence: " + pr.getPeptideEvidenceId()
                    + " the end coordinates couldn't be mapped.");
            mappedEnd = mappedStart + ((end - start) * 3);
            System.out.println("mapped_end modified value= " + mappedEnd);
            unmappedAccessions.add(pr.getAccession());

        }
        gffEntry[0] = mappedStart;
        gffEntry[1] = mappedEnd;

        return gffEntry;
    }

    /**
     * Get mapped coordinate.
     *
     * @param number  number
     * @param cdsColl CDs column list
     * @param pr      protein results
     *
     * @return a coordinate
     */
    long getMappedCordinates(long number, List<CdsInformation> cdsColl,
                             ProteinResults pr) {
        long mappedCord = -1;

        // compute the cumm array
        long[] cummArray = cumulativeStartPositions(cdsColl);
        //long number_toMap = (number - 1) * 3;
        long numberToMap = number * 3;

        int idx = determineTheIndexInCummulativeArray(numberToMap, cummArray);

        if (idx != -1) {

            long shift = cummArray[idx] - numberToMap;

            if (cdsColl.get(0).getStrand().contains("+")) {
                long endCds = cdsColl.get(idx).getEnd();
                mappedCord = endCds - shift;
            } else {
                long startCds = cdsColl.get(idx).getStart();
                mappedCord = startCds + shift;
            }
        }

        return mappedCord;
    }

    /**
     * Sort the CDS record, with the ascending or descending order of the start
     * position.
     * If the CDS is on +ve strand then sort in ascending order,
     * otherwise, sort in descending order.
     *
     * @param cdsColl CDS column list
     *
     * @return list of CDS information
     */
    List<CdsInformation> sortCdsAccordingToStartPosition(
            List<CdsInformation> cdsCollection) {

        List<CdsInformation> cdsColl = new ArrayList<>(cdsCollection);

        String strand = cdsCollection.get(0).getStrand();
        boolean sortAscending = false;

        if (strand.contains("+")) {
            sortAscending = true;
        }

        // Ascending sort...
        for (int i = 0; i < cdsColl.size() - 1; i++) {
            for (int j = i + 1; j < cdsColl.size(); j++) {
                CdsInformation cds1 = cdsColl.get(i);
                CdsInformation cds2 = cdsColl.get(j);
                if (cds1.getStart() > cds2.getStart()) {
                    CdsInformation temp = cds1;
                    cdsColl.set(i, cds2);
                    cdsColl.set(j, temp);
                }
            }
        }
        // If we need it in descending order, then reverse the sorted array
        if (!sortAscending) {
            for (int end = cdsColl.size() - 1; end >= 0; end--) {
                int start = cdsColl.size() - 1 - end;
                if (start < end) {
                    CdsInformation temp = cdsColl.get(end);
                    cdsColl.set(end, cdsColl.get(start));
                    cdsColl.set(start, temp);
                }
            }
        }
        return cdsColl;

    }

    List<CdsInformation> sortCdsAccordingToStartPositionTemp(
            List<CdsInformation> cdsCollection) {

        List<CdsInformation> cdsColl = new ArrayList<>(cdsCollection);

        for (int i = cdsColl.size() - 1; i >= 0; i--) {
            int j = cdsColl.size() - 1 - i;
            if (j < i) {
                CdsInformation temp = cdsColl.get(i);
                cdsColl.set(i, cdsColl.get(j));
                cdsColl.set(j, temp);
            }
        }
        return cdsColl;

    }

    long[] cumulativeStartPositions(List<CdsInformation> cdsCollection) {
        long[] cummStartPosition = new long[cdsCollection.size()];

        // compute Diff = end - start
        for (int i = 0; i < cdsCollection.size(); i++) {
            
            cummStartPosition[i] = cdsCollection.get(i).getEnd()
                    - cdsCollection.get(i).getStart() + 1;
        }

        // Form cumulative array 
        for (int i = 1; i < cummStartPosition.length; i++) {
            cummStartPosition[i] = cummStartPosition[i] + cummStartPosition[i
                    - 1];
        }
        return cummStartPosition;
    }

    int determineTheIndexInCummulativeArray(long number, long[] cummArray) {
        int idx = -1;

        try {
            for (int i = 0; i < cummArray.length; i++) {
                if (number <= cummArray[i]) {
                    idx = i;
                    break;
                }
            }
        } catch (Exception e) {
            String methodName = Thread.currentThread().getStackTrace()[1]
                    .getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + e
                    .getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information "
                    + "on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

        return idx;
    }

    /**
     * Get protein results from input mzid file.
     *
     * @param summaryFile input mzid file name.
     *
     * @return map of PeptideEvidence id to ProteinResults object.
     */
    public Map<String, ProteinResults> getProteinResults(String summaryFile) {
        Map<String, ProteinResults> prList = new HashMap<>();
        ProteinResults pr = null;

        Map<String, DBSequence> dbSequenceHashMap = new HashMap<>();
        Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
        Map<String, Peptide> peptideHashMap = new HashMap<>();
        try {
            MzIdentMLUnmarshaller mzIdentMlUnmarshaller
                    = new MzIdentMLUnmarshaller(new File(summaryFile));
            Iterator<DBSequence> iterDdSequence = mzIdentMlUnmarshaller
                    .unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDdSequence.hasNext()) {
                DBSequence dbSequence = iterDdSequence.next();
                dbSequenceHashMap.put(dbSequence.getId(), dbSequence);
            }

            Iterator<PeptideEvidence> iterPeptideEvidence
                    = mzIdentMlUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence pe = iterPeptideEvidence.next();
                peptideEvidenceHashMap.put(pe.getId(), pe);
            }
            Iterator<Peptide> iterPeptide = mzIdentMlUnmarshaller
                    .unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
                peptideHashMap.put(peptide.getId(), peptide);
            }

            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult
                    = mzIdentMlUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult
                        = iterSpectrumIdentificationResult.next();
                for (int j = 0; j < spectrumIdentificationResult
                        .getSpectrumIdentificationItem().size(); j++) {
                    SpectrumIdentificationItem spectrumIdentItem
                            = spectrumIdentificationResult
                            .getSpectrumIdentificationItem().get(j);

                    if (spectrumIdentItem.isPassThreshold()) {
                        List<CvParam> cvParamList = spectrumIdentItem
                                .getCvParam();
                        String fdrScore = "";
                        for (int i = 0; i < cvParamList.size(); i++) {
                            CvParam cvParam = cvParamList.get(i);
                            if (cvParam.getAccession().equals("MS:1002356")) {
                                fdrScore = cvParam.getValue();
                                break;
                            }
                        }
                        List<PeptideEvidenceRef> peptideEvidenceRefList
                                = spectrumIdentItem.getPeptideEvidenceRef();

                        for (PeptideEvidenceRef peRef : peptideEvidenceRefList) {
                            PeptideEvidence peptideEvidence
                                    = peptideEvidenceHashMap.get(peRef
                                            .getPeptideEvidenceRef());
                            DBSequence dbSequence = dbSequenceHashMap.get(
                                    peptideEvidence.getDBSequenceRef());
                            String accession = dbSequence.getAccession();
                            if (accession.startsWith("generic")) {
                                accession = accession.split("\\|")[1];
                                accession = accession.substring(2);
                                // Update by Fawaz Ghali 04/08/2016 fix Ensembl mapping version 85
                                if (accession.contains(".")) {
                                    accession = accession.substring(0, accession
                                                                    .indexOf("."));
                                }
                            } else {
                                throw new RuntimeException(
                                        "Not a generic accession");
                            }
                            Peptide peptide = peptideHashMap.get(
                                    peptideEvidence.getPeptideRef());
                            pr = new ProteinResults(accession, peptideEvidence
                                                    .isIsDecoy(),
                                                    peptideEvidence.getStart(),
                                                    peptideEvidence.getEnd(),
                                                    peptideEvidence.getId(),
                                                    fdrScore, peptide
                                                    .getPeptideSequence());
                            prList.put(peptideEvidence.getId(), pr);
                        }
                    }
                }
            }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1]
                    .getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + error
                    .getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information "
                    + "on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        return prList;
    }

}
