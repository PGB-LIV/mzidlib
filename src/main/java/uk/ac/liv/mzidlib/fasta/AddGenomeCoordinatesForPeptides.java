
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
import uk.ac.ebi.jmzidml.model.MzIdentMLObject;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.gff.CDS_Information;
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

    private List unmappedAccessions = new ArrayList();

    private Map<String, List<CDS_Information>> cdsRecords = new HashMap<>();

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

                        String[] gffArrayLine = line.split("\t");
                        //ArrayList<CDS_Information> gffDetails = new ArrayList();
                        String seqId = gffArrayLine[0];
                        String source = gffArrayLine[1];
                        String type = gffArrayLine[2];
                        String start = gffArrayLine[3];
                        String end = gffArrayLine[4];
                        String score = gffArrayLine[5];
                        String strand = gffArrayLine[6];
                        String phase = gffArrayLine[7];
                        String attribute = gffArrayLine[8];
                        long startPos;
                        long endPos;
                        startPos = Long.parseLong(start);
                        endPos = Long.parseLong(end);
                        CDS_Information cdsObj = new CDS_Information(seqId,
                                                                     source,
                                                                     startPos,
                                                                     endPos,
                                                                     strand,
                                                                     phase,
                                                                     attribute);
                        String accession = "";

                        String[] splits = attribute.split(";");
                        for (int i = 0; i < splits.length; i++) {
                            String string = splits[i];
                            if (string.contains("Parent=transcript:")) {
                                accession = string.replace("Parent=transcript:",
                                                           "");
                            }
                        }

                        accession = accession.toLowerCase();

                        List<CDS_Information> cdsColl = cdsRecords.
                                get(accession);

                        if (cdsColl == null) {
                            cdsColl = new ArrayList<>();
                        }
                        cdsColl.add(cdsObj);
                        cdsRecords.put(accession, cdsColl);
                    }
                    if (lineArray.length == 9 && lineArray[2].equals("CDS")) {
                        out.write(line + "\n");

                        String[] gffArrayLine = line.split("\t");
                        //ArrayList<CDS_Information> gffDetails = new ArrayList();
                        String seqId = gffArrayLine[0];
                        String source = gffArrayLine[1];
                        String type = gffArrayLine[2];
                        String start = gffArrayLine[3];
                        String end = gffArrayLine[4];
                        String score = gffArrayLine[5];
                        String strand = gffArrayLine[6];
                        String phase = gffArrayLine[7];
                        String attribute = gffArrayLine[8];
                        long startPos;
                        long endPos;
                        startPos = Long.parseLong(start);
                        endPos = Long.parseLong(end);
                        CDS_Information cdsObj = new CDS_Information(seqId,
                                                                     source,
                                                                     startPos,
                                                                     endPos,
                                                                     strand,
                                                                     phase,
                                                                     attribute);
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
                                    String[] protein_id = string.split("=");
                                    accession = protein_id[1];
                                }
                            }
                        }
                        accession = accession.toLowerCase();
                        List<CDS_Information> cdsColl = cdsRecords.
                                get(accession);

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
            Iterator it = proteinHits.entrySet().iterator();
            while (it.hasNext()) {
                Map.Entry pair = (Map.Entry) it.next();
                ProteinResults pr = (ProteinResults) pair.getValue();
                String[] co_ords = new String[2];
                long start_map, end_map;

                co_ords = mapToGff(pr);
                if (co_ords != null) {
                    // Take care of the negative strand reporting in GFF
                    start_map = Long.parseLong(co_ords[0]);
                    end_map = Long.parseLong(co_ords[1]);
                    if (end_map < start_map) {
                        long tmp = start_map;
                        start_map = end_map;
                        end_map = tmp;
                    }

                    String accession = pr.getAccession();
                    accession = accession.toLowerCase();
                    List<CDS_Information> cdsColl = cdsRecords.get(accession);

                    String feature = "peptide";
                    String score = ".";
                    String mappedStartLocation = Long.toString(start_map);
                    String mappedEndLocation = Long.toString(end_map);

                    String attr = "ID=pep_" + accession + "_"
                            + ";description=peptide_seq:" + pr.getPeptideSeq()
                            + "_fdr_score:" + pr.getfdrScore()
                            + "_peptide_start: " + pr.getstart()
                            + "_peptide_end: " + pr.getEnd() + ";Derives_from="
                            + accession;
                    // Just use the information from the first CDS entry with mapped start and end locations
                    String gffEntry = cdsColl.get(0).getSeqID() + "\t"
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
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller
                    = new MzIdentMLUnmarshaller(new File(inputMzid));

            CvList cvList = mzIdentMLUnmarshaller.unmarshal(
                    MzIdentMLElement.CvList);
            Iterator<Cv> iterCv = mzIdentMLUnmarshaller.
                    unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                }
            }
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            AnalysisSoftwareList analysisSoftwareList = mzIdentMLUnmarshaller.
                    unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            AuditCollection auditCollection = mzIdentMLUnmarshaller.unmarshal(
                    MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            Provider provider = mzIdentMLUnmarshaller.unmarshal(
                    MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            AnalysisProtocolCollection analysisProtocolCollection
                    = mzIdentMLUnmarshaller.unmarshal(
                            MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            AnalysisCollection analysisCollection = mzIdentMLUnmarshaller.
                    unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            Inputs inputs = mzIdentMLUnmarshaller.unmarshal(
                    MzIdentMLElement.Inputs);

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
                spectrumIdentificationListRef = analysisCollection.
                        getSpectrumIdentification().get(0).
                        getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Iterator<SpectrumIdentificationResult> sirIter
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (sirIter.hasNext()) {

                SpectrumIdentificationResult sr = sirIter.next();
                siList.getSpectrumIdentificationResult().add(sr);

            }

            String proteinDetectionListRef = "";
            if (analysisCollection.getProteinDetection() != null) {
                proteinDetectionListRef = analysisCollection.
                        getProteinDetection().getProteinDetectionListRef();
            }
            ProteinDetectionList pdList;
            pdList = new ProteinDetectionList();

            pdList.setId(proteinDetectionListRef);

            Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {
                ProteinAmbiguityGroup pag = iterProteinAmbiguityGroup.next();
                pdList.getProteinAmbiguityGroup().add(pag);
            }

            // here 
            SequenceCollection sequenceCollection1 = mzIdentMLUnmarshaller.
                    unmarshal(MzIdentMLElement.SequenceCollection);

            SequenceCollection sq = new SequenceCollection();
            List<DBSequence> dbSeqList = sequenceCollection1.getDBSequence();
            HashMap dbSeqHM = new HashMap();
            for (int i = 0; i < dbSeqList.size(); i++) {
                DBSequence get = dbSeqList.get(i);
                dbSeqHM.put(get.getId(), get);

            }
            sq.getPeptide().addAll(sequenceCollection1.getPeptide());

            List<PeptideEvidence> peList = sequenceCollection1.
                    getPeptideEvidence();
            List<PeptideEvidence> newPeList = new ArrayList<>();
            for (int i = 0; i < peList.size(); i++) {
                PeptideEvidence peptideEvidence = peList.get(i);
                DBSequence dbSeq = mzIdentMLUnmarshaller.unmarshal(
                        DBSequence.class, peptideEvidence.getDBSequenceRef());
                String accession = dbSeq.getAccession();
                if (accession.startsWith("generic")) {
                    accession = accession.split("\\|")[1];
                    accession = accession.substring(2);
                    // Update by Fawaz Ghali 04/08/2016 fix Ensembl mapping version 85
                    if (accession.contains(".")) {

                        accession = accession.substring(0, accession.
                                                        indexOf("."));
                    }
                } else {
                    throw new RuntimeException("Not generic accession");
                }
                accession = accession.toLowerCase();
                List<CDS_Information> gffData = cdsRecords.get(accession);
                if (gffData != null && gffData.size() > 0 && !peptideEvidence.
                        isIsDecoy()) {
                    ProteinResults pr = prList.get(peptideEvidence.getId());
                    String[] co_ords = new String[2];
                    long start_map, end_map;

                    co_ords = mapToGff(pr);

                    // unmapped peptide
                    if (co_ords == null || unmappedAccessions.contains(pr.
                            getAccession())) {
                        dbSeq.getCvParam().add(makeCvParam("MS:1002741",
                                                           "unmapped protein",
                                                           gffData.get(0).
                                                           getSeqID(), psiCV));
                        peptideEvidence.getCvParam().add(makeCvParam(
                                "MS:1002740", "unmapped peptide",
                                gffData.get(0).getSeqID(), psiCV));
                    } // mapped peptide
                    else {
                        start_map = Long.parseLong(co_ords[0]);
                        end_map = Long.parseLong(co_ords[1]);
                        if (end_map < start_map) {
                            long tmp = start_map;
                            start_map = end_map;
                            end_map = tmp;
                        }
                        // Added by Fawaz Ghali 13/8/2015
                        // Update the code to cover the case where a peptide covers multiple exons
                        List<CDS_Information> sortedCDS
                                = sortCDSAccordingToStartPosition(gffData);
                        List<CDS_Information> outputCDS = new ArrayList();
                        int countCDS = 0;
                        boolean sorted = false;
                        if (gffData.get(0).getStrand().equals("-")) {
                            sortedCDS = sortCDSAccordingToStartPositionTemp(
                                    sortedCDS);
                        }
                        for (int j = 0; j < sortedCDS.size(); j++) {
                            CDS_Information object = sortedCDS.get(j);
                            long s = object.getStart();
                            long e = object.getEnd();
                            long newStart, newEnd;
                            countCDS = j;
                            if (start_map > e) {
                                continue;
                            } else {
                                newStart = start_map;
                                if (end_map <= e) {
                                    newEnd = end_map;
                                    sorted = true;
                                } else {
                                    newEnd = e;
                                }

                                CDS_Information temp = new CDS_Information("",
                                                                           "",
                                                                           newStart,
                                                                           newEnd,
                                                                           "",
                                                                           "",
                                                                           "");
                                outputCDS.add(temp);
                                break;
                            }
                        }
                        if (!sorted) {
                            for (int j = countCDS + 1; j < sortedCDS.size(); j++) {
                                CDS_Information object = sortedCDS.get(j);
                                long s = object.getStart();
                                long e = object.getEnd();
                                long newStart, newEnd;

                                if (end_map > e) {
                                    CDS_Information temp = new CDS_Information(
                                            "", "", s, e, "", "", "");
                                    outputCDS.add(temp);
                                } else {
                                    newStart = s;
                                    newEnd = end_map;

                                    CDS_Information temp = new CDS_Information(
                                            "", "", newStart, newEnd, "", "", "");
                                    outputCDS.add(temp);
                                    break;
                                }
                            }
                        }
                        int pepLen = 3 * (peptideEvidence.getEnd()
                                - peptideEvidence.getStart() + 1);
                        int coordLen = 0;
//                        for (int j = 0; j < outputCDS.size(); j++) {
//                            CDS_Information cDS_Information = outputCDS.get(j);
//                            cDS_Information.getStart();
//                            cDS_Information.getEnd();
//                            coordLen = coordLen + toIntExact(cDS_Information.getEnd() - cDS_Information.getStart());
//                        }
                        long s, e = 0;
                        String startsList = "";
                        String positionsList = "";

                        //System.out.println("END");
                        // System.out.println("pepLen "+ pepLen);
                        //System.out.println("outputCDS.size() "+ outputCDS.size());
                        for (int j = 0; j < outputCDS.size(); j++) {
                            CDS_Information cDS_Information = outputCDS.get(j);
                            s = cDS_Information.getStart() - 1;
                            // peptideEvidence.getUserParam().add(makeUserParam("start_map", String.valueOf(s)));
                            startsList = startsList + s + ",";
                            if (j == 0) {
                                e = cDS_Information.getEnd() + 2;
                            } else {
                                e = cDS_Information.getEnd();
                            }
                            // peptideEvidence.getUserParam().add(makeUserParam("end_map", String.valueOf(e)));
                            coordLen = coordLen + MzidLibUtils.safeLongToInt(e
                                    - s);
                            positionsList = positionsList + String.
                                    valueOf(e - s) + ",";
//                            System.out.println("toIntExact(e -s) "+ toIntExact(e -s));
                        }
//                            System.out.println("coordLen "+ coordLen);
//                        if (pepLen==coordLen){
//                            System.out.println("PASS "+ peptideEvidence.getPeptideRef() +" pepLen "+pepLen +" coordLen "+ coordLen);
//                        }else{
//                            System.out.println("FAILED "+ peptideEvidence.getPeptideRef()+" pepLen "+pepLen +" coordLen "+ coordLen);
//                        }
//                        System.out.println("");

//                        peptideEvidence.getUserParam().add(makeUserParam("chr", gffData.get(0).getSeqID()));
                        dbSeq.getCvParam().add(makeCvParam("MS:1002637",
                                                           "chromosome name",
                                                           gffData.get(0).
                                                           getSeqID(), psiCV));
//                        peptideEvidence.getUserParam().add(makeUserParam("strand", gffData.get(0).getStrand()));
                        dbSeq.getCvParam().add(makeCvParam("MS:1002638",
                                                           "chromosome strand",
                                                           gffData.get(0).
                                                           getStrand(), psiCV));
                        dbSeq.getCvParam().add(makeCvParam("MS:1002644",
                                                           "genome reference version",
                                                           new File(inputGff).
                                                           getName(), psiCV));
                        dbSeqHM.put(dbSeq.getId(), dbSeq);
                        //peptideEvidence.getCvParam().add(makeCvParam("MS:1002639", "peptide start on chromosome", String.valueOf(start_map), psiCV));
                        peptideEvidence.getCvParam().add(makeCvParam(
                                "MS:1002640", "peptide end on chromosome",
                                String.valueOf(end_map), psiCV));

                        peptideEvidence.getCvParam().add(makeCvParam(
                                "MS:1002641", "peptide exon count", String.
                                valueOf(outputCDS.size()), psiCV));
                        if (startsList.endsWith(",")) {
                            startsList = startsList.substring(0, startsList.
                                                              length() - 1);
                        }
                        if (positionsList.endsWith(",")) {
                            positionsList = positionsList.substring(0,
                                                                    positionsList.
                                                                    length() - 1);
                        }
                        peptideEvidence.getCvParam().add(makeCvParam(
                                "MS:1002642", "peptide exon nucleotide sizes",
                                positionsList, psiCV));
                        peptideEvidence.getCvParam().add(makeCvParam(
                                "MS:1002643",
                                "peptide start positions on chromosome",
                                startsList, psiCV));

                    }
                } else if (gffData == null && !peptideEvidence.isIsDecoy()) {
                    dbSeq.getCvParam().add(makeCvParam("MS:1002741",
                                                       "unmapped protein", psiCV));
                    peptideEvidence.getCvParam().add(makeCvParam("MS:1002740",
                                                                 "unmapped peptide",
                                                                 psiCV));
                    dbSeqHM.put(dbSeq.getId(), dbSeq);
                }
                newPeList.add(peptideEvidence);
            }
            sq.getDBSequence().addAll(dbSeqHM.values());
            sq.getPeptideEvidence().addAll(newPeList);

            marshaller.marshal(sq, writer);

            writer.write("\n");

            marshaller.marshal(analysisCollection, writer);

            writer.write("\n");

            if (analysisProtocolCollection != null) {

                List<SpectrumIdentificationProtocol> sipList
                        = analysisProtocolCollection.
                        getSpectrumIdentificationProtocol();
                for (SpectrumIdentificationProtocol sip : sipList) {
                    List<CvParam> cvs = sip.getAdditionalSearchParams().
                            getCvParam();
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
        } catch (JAXBException ex) {
            ex.printStackTrace();
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + e.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    public UserParam makeUserParam(String name, String value) {
        UserParam userParam = new UserParam();
        userParam.setName(name);
        userParam.setValue(value);
        return userParam;
    }

    public CvParam makeCvParam(String accession, String name, String value,
                               Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setValue(value);
        cvParam.setCv(cv);
        return cvParam;
    }

    public CvParam makeCvParam(String accession, String name, Cv cv) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        cvParam.setCv(cv);
        return cvParam;
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
        List<CDS_Information> cdsColl = cdsRecords.get(accession);

        long[] mapped_cords = determineTheLocationOfSeqOnCds(pr, cdsColl);

        gffMapping[0] = Long.toString(mapped_cords[0]);
        gffMapping[1] = Long.toString(mapped_cords[1]);
        gffMapping[2] = cdsRecords.get(accession).get(0).getSeqID();

        return gffMapping;
    }

    /**
     * Determine the location of sequence on CDs.
     *
     * @param pr      protein results
     * @param cdsColl CDS column list
     *
     * @return an array of locations
     */
    long[] determineTheLocationOfSeqOnCds(ProteinResults pr,
                                          List<CDS_Information> cdsColl) {
        long[] gffEntry = new long[2];

        // Get locations from the protein object
        long start = pr.getstart();
        long end = pr.getEnd();

        // sort the CDS collection according to the strand
        List<CDS_Information> sortedCDS = sortCDSAccordingToStartPosition(
                cdsColl);

        long mapped_start = getMappedCordinates(start, sortedCDS, pr);
        long mapped_end = getMappedCordinates(end, sortedCDS, pr);
        if (mapped_start == -1) {
            System.out.println("For accession: " + pr.getAccession()
                    + " and a peptide eveidence: " + pr.getPeptideEvidenceID()
                    + " start coordinates couldn't be mapped.");
            unmappedAccessions.add(pr.getAccession());
        }
        if (mapped_end == -1) {
            System.out.println("For accession: " + pr.getAccession()
                    + " and a peptide eveidence: " + pr.getPeptideEvidenceID()
                    + " the end coordinates couldn't be mapped.");
            mapped_end = mapped_start + ((end - start) * 3);
            System.out.println("mapped_end modified value= " + mapped_end);
            unmappedAccessions.add(pr.getAccession());

        }

        gffEntry[0] = mapped_start;
        gffEntry[1] = mapped_end;

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
    long getMappedCordinates(long number, List<CDS_Information> cdsColl,
                             ProteinResults pr) {
        long mappedCord = -1;

        // compute the cumm array
        long[] cummArray = cumulativeStartPositions(cdsColl);
        //long number_toMap = (number - 1) * 3;
        long number_toMap = number * 3;

        int idx = determineTheIndexInCummulativeArray(number_toMap, cummArray);

        if (idx != -1) {

            long shift = cummArray[idx] - number_toMap;

            if (cdsColl.get(0).getStrand().contains("+")) {
                long end_cds = cdsColl.get(idx).getEnd();
                mappedCord = end_cds - shift;
            } else {
                long start_cds = cdsColl.get(idx).getStart();
                mappedCord = start_cds + shift;
            }
        }

        return mappedCord;
    }

    /**
     * Sort the CDS record, with the ascending or descending order of the start
     * position. If the CDS is on +ve strand then sort in ascending order,
     * otherwise, sort in descending order.
     *
     * @param cdsColl CDS column list
     *
     * @return list of CDS information
     */
    List<CDS_Information> sortCDSAccordingToStartPosition(
            List<CDS_Information> cdsCollection) {

        List<CDS_Information> cdsColl = new ArrayList<>(cdsCollection);

        String strand = cdsCollection.get(0).getStrand();
        boolean sortAscending = false;

        if (strand.contains("+")) {
            sortAscending = true;
        }

        // Ascending sort...
        for (int i = 0; i < cdsColl.size() - 1; i++) {
            for (int j = i + 1; j < cdsColl.size(); j++) {
                CDS_Information cds_i = cdsColl.get(i);
                CDS_Information cds_j = cdsColl.get(j);
                if (cds_i.getStart() > cds_j.getStart()) {
                    CDS_Information temp = cds_i;
                    cdsColl.set(i, cds_j);
                    cdsColl.set(j, temp);
                }
            }
        }

        // If we need it in descending order, then reverse the sorted array
        if (!sortAscending) {
            for (int i = cdsColl.size() - 1; i >= 0; i--) {
                int j = cdsColl.size() - 1 - i;
                if (j < i) {
                    CDS_Information temp = cdsColl.get(i);
                    cdsColl.set(i, cdsColl.get(j));
                    cdsColl.set(j, temp);
                }
            }
        }

        return cdsColl;

    }

    List<CDS_Information> sortCDSAccordingToStartPositionTemp(
            List<CDS_Information> cdsCollection) {

        List<CDS_Information> cdsColl = new ArrayList<>(cdsCollection);

        for (int i = cdsColl.size() - 1; i >= 0; i--) {
            int j = cdsColl.size() - 1 - i;
            if (j < i) {
                CDS_Information temp = cdsColl.get(i);
                cdsColl.set(i, cdsColl.get(j));
                cdsColl.set(j, temp);
            }
        }

        return cdsColl;

    }

    long[] cumulativeStartPositions(List<CDS_Information> cdsCollection) {
        long[] cummStartPosition = new long[cdsCollection.size()];

        // compute Diff = end - start
        for (int i = 0; i < cdsCollection.size(); i++) {
            //cummStartPosition[i] = cdsCollection.get(i).getEnd() - cdsCollection.get(i).getStart();
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
            System.out.println(e.getMessage());
            e.printStackTrace();
        }

        return idx;
    }

    public Map<String, ProteinResults> getProteinResults(String summaryFile) {
        Map<String, ProteinResults> prList = new HashMap();
        ProteinResults pr = null;

        Map<String, DBSequence> dBSequenceHashMap = new HashMap<>();
        Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
        Map<String, Peptide> peptideHashMap = new HashMap<>();
        Map<String, String> peptideIdAndSequenceHash = new HashMap<>();
        try {
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller
                    = new MzIdentMLUnmarshaller(new File(summaryFile));
            Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.
                    unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {
                DBSequence dBSequence = iterDBSequence.next();
                dBSequenceHashMap.put(dBSequence.getId(), dBSequence);
            }

            Iterator<PeptideEvidence> iterPeptideEvidence
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence pe = iterPeptideEvidence.next();
                peptideEvidenceHashMap.put(pe.getId(), pe);
            }
            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.
                    unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();

                String pepId = peptide.getId();
                String pepSeq = peptide.getPeptideSequence();
                peptideIdAndSequenceHash.put(pepId, pepSeq);
                peptideHashMap.put(peptide.getId(), peptide);
            }

            Iterator<MzIdentMLObject> iterSpectrumIdentificationResult
                    = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult
                        = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.
                        next();
                for (int j = 0; j < spectrumIdentificationResult.
                        getSpectrumIdentificationItem().size(); j++) {
                    SpectrumIdentificationItem spectrumIdentItem
                            = spectrumIdentificationResult.
                            getSpectrumIdentificationItem().get(j);

                    if (spectrumIdentItem.isPassThreshold()) {
                        List<CvParam> cvParamList = spectrumIdentItem.
                                getCvParam();
                        String fdrScore = "";
                        for (int i = 0; i < cvParamList.size(); i++) {
                            CvParam cvParam = cvParamList.get(i);
                            if (cvParam.getAccession().equals("MS:1002356")) {
                                fdrScore = cvParam.getValue();
                                break;
                            }

                        }
                        List<PeptideEvidenceRef> peptideEvidenceRef
                                = spectrumIdentItem.getPeptideEvidenceRef();

                        for (int i = 0; i < peptideEvidenceRef.size(); i++) {
                            PeptideEvidenceRef peptideEvidenceRef1
                                    = peptideEvidenceRef.get(i);
                            PeptideEvidence peptideEvidence
                                    = peptideEvidenceHashMap.get(
                                            peptideEvidenceRef1.
                                            getPeptideEvidenceRef());
                            DBSequence dBSequence = dBSequenceHashMap.get(
                                    peptideEvidence.getDBSequenceRef());
                            String accession = dBSequence.getAccession();
                            if (accession.startsWith("generic")) {
                                accession = accession.split("\\|")[1];
                                accession = accession.substring(2);
                                // Update by Fawaz Ghali 04/08/2016 fix Ensembl mapping version 85
                                if (accession.contains(".")) {

                                    accession = accession.substring(0,
                                                                    accession.
                                                                    indexOf("."));
                                }
                            } else {
                                throw new RuntimeException(
                                        "Not a generic accession");
                            }
                            Peptide peptide = peptideHashMap.get(
                                    peptideEvidence.getPeptideRef());
                            pr = new ProteinResults(accession, peptideEvidence.
                                                    isIsDecoy(),
                                                    peptideEvidence.getStart(),
                                                    peptideEvidence.getEnd(),
                                                    peptideEvidence.getId(),
                                                    fdrScore, peptide.
                                                    getPeptideSequence());
                            prList.put(peptideEvidence.getId(), pr);

                        }
                    }

                }
            }

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of " + error.
                    getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        return prList;
    }

}
