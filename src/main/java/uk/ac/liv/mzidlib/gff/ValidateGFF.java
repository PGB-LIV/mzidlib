/**
 * The code to process a GFF3 file for creating FASTA and for creating a new
 * GFF3 with mapped peptides added. The code works in two different modes which
 * is listed in the main().
 *
 * The code is pretty strict in parsing the GFF3 format, as it demands the GFF3
 * to have certain fields - - ## FASTA - the third column in GFF3 denoting type
 * must be equal to string "CDS" denoting coding region - The ninth column in
 * GFF3 denoting attribute must have an "ID" - The FASTA sequences in the GFF3
 * must have the same identifier as the corresponding ID of the CDS
 *
 * The code does the following ignoring and error reporting - - if no ##FASTA
 * section present in the file - report error and exit - if any CDS is missing
 * the "ID" in attribute - report error and exit - if extra fasta sequences
 * present than the ones with GFF3 entry - ignore extra sequences, continue - if
 * extra CDS present in GFF3 entries, but no corresponding fasta seq - ignore
 * those CDS, continue
 */
package uk.ac.liv.mzidlib.gff;

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
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

public class ValidateGFF {

    Map<String, List<CDS_Information>> cdsRecords;
    String nameOfTempFastaFile;
    String extOfTempFastaFile;
    File temp;
    private MzidLibUtils mzidLibUtils;

    public ValidateGFF() {
        cdsRecords = new HashMap<String, List<CDS_Information>>();
        nameOfTempFastaFile = "TempFastaFromGFF";
        extOfTempFastaFile = "fasta";
        mzidLibUtils = new MzidLibUtils();
    }

    /**
     * Processes the GFF file and creates a temp FASTA file containing the
     * sequences from the ##FASTA section of the GFF
     *
     * @param inputGff
     * @return Name of the temp Fasta file created
     * @throws Exception
     */
    public String processGffFileAndCreateFasta(String gffFile) throws Exception {

        try {
            BufferedReader in = new BufferedReader(new FileReader(gffFile));
            String line;
            boolean fastaRegionFlag = false;

            // Create temp file for storing the Fasta Sequences in the GFF file.
            // REMOVED by Fawaz
            //temp = File.createTempFile(nameOfTempFastaFile, extOfTempFastaFile);
            // Added by Fawaz: store the fasta file in same directory of the gff
            temp = new File(gffFile.substring(0, gffFile.lastIndexOf('.')) + ".fasta");
            // 
            //temp.deleteOnExit();		   

            Writer temp_out = new BufferedWriter(new FileWriter(temp));

            // Variables needed during processing of ##FASTA file
            String accn = new String();
            String seq = new String();
            int fastaCount = 0;

            while ((line = in.readLine()) != null) {

                //System.out.println(line);
                if (line.startsWith("##")) {
                    if (line.equals("##FASTA")) {
                        fastaRegionFlag = true;
                    }
                }

                if (fastaRegionFlag == false) {

                    String[] records = line.split("\t", -1);
                    if (records.length != 9) {
                        System.out.println("Skipping... - " + line);
                        continue;
                    }

                    String seqId = records[0];
                    String source = records[1];
                    String type = records[2];
                    String start = records[3];
                    String end = records[4];
                    String score = records[5];
                    String strand = records[6];
                    String phase = records[7];
                    String attribute = records[8];

                    if (type.trim().matches("CDS")) {
                        long startPos;
                        long endPos;
                        try {
                            startPos = Long.parseLong(start);
                            endPos = Long.parseLong(end);
                            //A simple check...
                            //if(endPos < startPos){
                            //	System.out.println("End position = " + endPos + " greater than start position " + startPos+ " in Gff file");
                            //	throw new Exception();
                            //}
                        } catch (NumberFormatException e) {
                            String exMessage = "ProteoAnnotator : Problem in processing the start and end fields of Line - \n" + line + "\n...exiting";
                            System.out.println(exMessage);
                            in.close();
                            temp_out.close();
                            throw new Exception(exMessage);
                        }

                        CDS_Information cdsObj = new CDS_Information(seqId, source, startPos, endPos, strand, phase, attribute);
                        String cdsId = parseAttributeFieldToExtractId(attribute);

                        List<CDS_Information> cdsColl = cdsRecords.get(cdsId);

                        if (cdsColl == null) {
                            cdsColl = new ArrayList<CDS_Information>();
                        }

                        cdsColl.add(cdsObj);

                        if (cdsId != null) {
                            this.cdsRecords.put(cdsId, cdsColl);
                        }
                    }

                } else {

                    // Now we are in the FASTA section, so we should have all the cds information
                    if (cdsRecords.size() == 0) {
                        String exMessage = "No CDS field found.....exiting";
                        System.out.println(exMessage);
                        in.close();
                        temp_out.close();
                        throw new Exception(exMessage);
                    }

                    // If all well, then process fasta - skip the line having ##FASTA 
                    if (line.equals("##FASTA")) {
                        continue;
                    }

                    // Write all other lines
                    if (line.contains(">")) {
                        boolean cdsAccnFound = verifyNeededAccnOrNot(accn);

                        // Write to file only if it is a cds Accession
                        if (cdsAccnFound) {
                            temp_out.write("\n" + accn + "\n" + seq);
                            fastaCount++;
                        }
                        accn = line;
                        seq = "";
                    } else {
                        seq = seq.concat(line);
                    }

                }

            } // end of while

            //System.out.println("Out of while loop");
            // Signal error if no ##FASTA is encountered...
            if (fastaRegionFlag == false) {
                String exMessage = "No ##FASTA directive encountered.....Quitting";
                System.out.println(exMessage);
                in.close();
                temp_out.close();
                throw new Exception(exMessage);
            }

            // Write the last fasta sequence in the file
            if (line == null) {
                boolean cdsAccnFound = verifyNeededAccnOrNot(accn);
                if (cdsAccnFound) {
                    temp_out.write("\n" + accn + "\n" + seq);
                    fastaCount++;
                }
            }

            System.out.println("Total CDS entries found - " + cdsRecords.size());

            in.close();
            temp_out.close();

            // return the name of the file with fasta sequences 
            return temp.getAbsolutePath();

        } catch (Exception e) {
            // All the exceptions caught during making of a FASTA file are reported here and the program exists
            System.out.println("ProteoAnnotator : Exception encountered in the GFF - " + e.getMessage());
            throw e;
        }

    }

    /**
     *
     * @param attribute
     * @return
     */
    String parseAttributeFieldToExtractId(String attribute) throws Exception {
        String id = new String();

        String[] splitOnColon = attribute.split(";", -1);
        if (splitOnColon.length == 0) {
            String[] keyVal = attribute.split("=", -1);
            if (keyVal.length == 0) {
                System.out.println("No ID found for the CDS entry");
            } else if (keyVal[0].compareToIgnoreCase("ID") < 0) {
                System.out.println("No ID found for the CDS entry");
            } else {
                id = keyVal[1];
            }
        } else {
            boolean found = false;
            for (int i = 0; i < splitOnColon.length; i++) {

                String subAttr = splitOnColon[i];
                String[] keyVal = subAttr.split("=", -1);
                if (keyVal.length == 0) {
                    System.out.println("No ID found for the CDS entry");
                } else if (keyVal[0].compareToIgnoreCase("ID") < 0) {
                    System.out.println("No ID found for the CDS entry");
                } else {
                    id = keyVal[1];
                    found = true;
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                System.out.println("No ID found for the CDS entry");
                throw new Exception();
            }

        }

        return id;
    }

    /**
     *
     * @param accnToCheck
     * @return
     */
    boolean verifyNeededAccnOrNot(String accnToCheck) {
        boolean found = false;

        // Remove the ">" from accnToCheck
        if (accnToCheck.contains(">")) {
            accnToCheck = accnToCheck.replace(">", "").trim();
        }

        found = cdsRecords.containsKey(accnToCheck);

        return found;
    }

    /**
     *
     * @param inputGffFile
     * @param outputGff
     * @param proteinHits
     */
    public void writingTheGFFfile(String inputGffFile, String outputGffFile, List<ProteinResults> proteinHits)
            throws Exception {
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

                    String[] records = line.split("\t", -1);

                    out.write(line + "\n");

                    if (records.length != 9) {
                        System.out.println("Skipping... - " + line);
                        continue;
                    }

                    String seqId = records[0];
                    String source = records[1];
                    String type = records[2];
                    String start = records[3];
                    String end = records[4];
                    String score = records[5];
                    String strand = records[6];
                    String phase = records[7];
                    String attribute = records[8];

                    if (type.trim().matches("CDS")) {
                        long startPos;
                        long endPos;
                        try {
                            startPos = Long.parseLong(start);
                            endPos = Long.parseLong(end);

                            //A simple check...
                            //if(endPos < startPos){
                            //	System.out.println("End position = " + endPos + " greater than start position " + startPos+ " in Gff file");
                            //	throw new Exception();
                            //}
                        } catch (NumberFormatException e) {
                            String exMessage = "ProteoAnnotator : Problem in processing the start and end fields of Line - \n" + line + "\n...exiting";
                            System.out.println(exMessage);
                            in.close();
                            out.close();
                            throw new Exception(exMessage);
                        }

                        CDS_Information cdsObj = new CDS_Information(seqId, source, startPos, endPos, strand, phase, attribute);
                        String cdsId = parseAttributeFieldToExtractId(attribute);

                        List<CDS_Information> cdsColl = cdsRecords.get(cdsId);

                        if (cdsColl == null) {
                            cdsColl = new ArrayList<CDS_Information>();
                        }

                        cdsColl.add(cdsObj);
                        this.cdsRecords.put(cdsId, cdsColl);
                    }

                }

                // Come out when ##FASTA encountered
                if (fastaRegionFlag) {
                    break;
                }

            } // end of while

            // Exit if no CDS found
            if (cdsRecords.size() == 0) {
                String exMessage = "No CDS field found.....exiting";
                System.out.println(exMessage);
                in.close();
                out.close();
                throw new Exception(exMessage);
            }

            in.close();

            // Perform Chromosome mapping
            for (int i = 0; i < proteinHits.size(); i++) {
                ProteinResults pr = proteinHits.get(i);
                String[] co_ords = new String[2];
                long start_map, end_map;
                try {
                    co_ords = mapToGff(pr);

                    // Take care of the negative strand reporting in GFF
                    start_map = Long.parseLong(co_ords[0]);
                    end_map = Long.parseLong(co_ords[1]);
                    if (end_map < start_map) {
                        long tmp = start_map;
                        start_map = end_map;
                        end_map = tmp;
                    }
                } catch (Exception e) {
                    //System.out.println("Exception - " + pr.getAccession() +" \t"+ pr.getstart() +" \t"+ pr.getEnd());
                    //e.printStackTrace();
                    continue;
                }

                String accession = pr.getAccession();
                List<CDS_Information> cdsColl = cdsRecords.get(accession);
                String feature = "peptide";
                String score = ".";
                String mappedStartLocation = Long.toString(start_map);
                String mappedEndLocation = Long.toString(end_map);

                String attr = "ID=pep_" + accession + "_" + i + ";description=peptide;Derives_from=" + accession;
                // Just use the information from the first CDS entry with mapped start and end locations
                String gffEntry = cdsColl.get(0).getSeqID() + "\t" + cdsColl.get(0).getSource() + "\t" + feature + "\t" + mappedStartLocation
                        + "\t" + mappedEndLocation + "\t" + score + "\t" + cdsColl.get(0).getStrand() + "\t" + cdsColl.get(0).getSePhase() + "\t" + attr + "\n";
                out.write(gffEntry);
            }
            out.close();

        } catch (IOException e) {
            System.out.println("ProteoAnnotator : IO Excpetion : Unable to create/read file");
            throw e;
        } catch (Exception e) {
            System.out.println("ProteoAnnotator : Exception encountered in the GFF - " + e.getMessage());
            throw e;
        }
    }

    /**
     * Map the co-ordinates Returns a String array - idx=0=start, idx=1=end and
     * idx=2=chr
     */
    public String[] mapToGff(ProteinResults pr) {
        String[] gffMapping = new String[3];

        // Get accession from the protein object
        String accession = pr.getAccession();

        // return if no key found
        if (!cdsRecords.containsKey(accession)) {
            return null;
        }

        //...find the cds this range falls into...
        List<CDS_Information> cdsColl = cdsRecords.get(accession);

        long[] mapped_cords = determineTheLocationOfSeqOnCds(pr, cdsColl);

        gffMapping[0] = Long.toString(mapped_cords[0]);
        gffMapping[1] = Long.toString(mapped_cords[1]);
        gffMapping[2] = cdsRecords.get(accession).get(0).seqid;

        return gffMapping;
    }

    /**
     *
     * @param pr
     * @param cdsColl
     * @return
     */
    long[] determineTheLocationOfSeqOnCds(ProteinResults pr, List<CDS_Information> cdsColl) {
        long[] gffEntry = new long[2];

        // Get locations from the protein object
        long start = pr.getstart();
        long end = pr.getEnd();

        // sort the CDS collection according to the strand
        List<CDS_Information> sortedCDS = sortCDSAccordingToStartPosition(cdsColl);

        long mapped_start = getMappedCordinates(start, sortedCDS);
        long mapped_end = getMappedCordinates(end, sortedCDS);

        gffEntry[0] = mapped_start;
        gffEntry[1] = mapped_end;

        return gffEntry;
    }

    /**
     *
     * @param number
     * @param cdsColl
     * @return
     */
    long getMappedCordinates(long number, List<CDS_Information> cdsColl) {
        long mappedCord = 0;

        // compute the cumm array
        long[] cummArray = cumulativeStartPositions(cdsColl);
        //long number_toMap = (number - 1) * 3;
        long number_toMap = number * 3;

        int idx = determineTheIndexInCummulativeArray(number_toMap, cummArray);

        long shift = cummArray[idx] - number_toMap;

        if (cdsColl.get(0).getStrand().contains("+")) {
            long end_cds = cdsColl.get(idx).getEnd();
            mappedCord = end_cds - shift;
        } else {
            long start_cds = cdsColl.get(idx).getStart();
            mappedCord = start_cds + shift;
        }

        return mappedCord;
    }

    /**
     * Sort the CDS record, with the ascending or descending order of the start
     * position. If the CDS is on +ve strand then sort in ascending order,
     * otherwise, sort in descending order.
     *
     * @param cdsColl
     * @return
     */
    List<CDS_Information> sortCDSAccordingToStartPosition(List<CDS_Information> cdsCollection) {

        List<CDS_Information> cdsColl = new ArrayList<CDS_Information>(cdsCollection);

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

    /**
     *
     * @param cdsCollection
     * @return
     */
    long[] cumulativeStartPositions(List<CDS_Information> cdsCollection) {
        long[] cummStartPosition = new long[cdsCollection.size()];

        // compute Diff = end - start
        for (int i = 0; i < cdsCollection.size(); i++) {
            //cummStartPosition[i] = cdsCollection.get(i).getEnd() - cdsCollection.get(i).getStart();
            cummStartPosition[i] = cdsCollection.get(i).getEnd() - cdsCollection.get(i).getStart() + 1;
        }

        // Form cumulative array 
        for (int i = 1; i < cummStartPosition.length; i++) {
            cummStartPosition[i] = cummStartPosition[i] + cummStartPosition[i - 1];
        }

        return cummStartPosition;
    }

    /**
     *
     * @param number
     * @param cummArray
     * @return
     */
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

    /**
     * Test function 1
     *
     * @param args
     * @throws Exception
     */
    public void createFastaFromGFF(String gffFile) throws Exception {

        System.out.println("Fasta File produced -" + processGffFileAndCreateFasta(gffFile));
        System.out.println("Press enter to exit...");
        System.in.read();

    }

    public List<ProteinResults> getProteinResults(String summaryFile, boolean passThreshold) {
        List<ProteinResults> prList = null;
        ProteinResults pr = null;

        Map<String, DBSequence> dBSequenceHashMap = new HashMap<>();
        Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
        Map<String, Peptide> peptideHashMap = new HashMap<>();
        Map<String, String> peptideIdAndSequenceHash = new HashMap<String, String>();
        try {
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(summaryFile));
            Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while (iterDBSequence.hasNext()) {
                DBSequence dBSequence = iterDBSequence.next();
                dBSequenceHashMap.put(dBSequence.getId(), dBSequence);
            }

            Iterator<PeptideEvidence> iterPeptideEvidence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
            while (iterPeptideEvidence.hasNext()) {
                PeptideEvidence pe = iterPeptideEvidence.next();
                peptideEvidenceHashMap.put(pe.getId(), pe);
            }
            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();

                String pepId = peptide.getId();
                String pepSeq = peptide.getPeptideSequence();
                peptideIdAndSequenceHash.put(pepId, pepSeq);
                peptideHashMap.put(peptide.getId(), peptide);
            }

            Iterator<MzIdentMLObject> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.next();
                for (int j = 0; j < spectrumIdentificationResult.getSpectrumIdentificationItem().size(); j++) {
                    SpectrumIdentificationItem spectrumIdentItem = spectrumIdentificationResult.getSpectrumIdentificationItem().get(j);
                    if (passThreshold) {
                        if (spectrumIdentItem.isPassThreshold()) {

                            List<PeptideEvidenceRef> peptideEvidenceRef = spectrumIdentItem.getPeptideEvidenceRef();

                            for (int i = 0; i < peptideEvidenceRef.size(); i++) {
                                PeptideEvidenceRef peptideEvidenceRef1 = peptideEvidenceRef.get(i);
                                PeptideEvidence peptideEvidence = peptideEvidenceHashMap.get(peptideEvidenceRef1);
                                DBSequence dBSequence = dBSequenceHashMap.get(peptideEvidence.getDBSequenceRef());
                                pr = new ProteinResults(dBSequence.getAccession(), peptideEvidence.isIsDecoy(), peptideEvidence.getStart(), peptideEvidence.getEnd(), peptideEvidence.getId(), "", "");
                                // TODO: Null Pointer here! prList never initialised
                                prList.add(pr);

                            }
                        }

                    } else {
                        List<PeptideEvidenceRef> peptideEvidenceRef = spectrumIdentItem.getPeptideEvidenceRef();

                        for (int i = 0; i < peptideEvidenceRef.size(); i++) {
                            PeptideEvidenceRef peptideEvidenceRef1 = peptideEvidenceRef.get(i);
                            PeptideEvidence peptideEvidence = peptideEvidenceHashMap.get(peptideEvidenceRef1);
                            DBSequence dBSequence = dBSequenceHashMap.get(peptideEvidence.getDBSequenceRef());
                            pr = new ProteinResults(dBSequence.getAccession(), peptideEvidence.isIsDecoy(), peptideEvidence.getStart(), peptideEvidence.getEnd(), peptideEvidence.getId(), "", "");
                            // TODO: Null Pointer here! prList never initialised
                            prList.add(pr);
                        }
                    }

                }
            }

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        return prList;
    }

    /**
     *
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {

        ValidateGFF vgf = new ValidateGFF();

        String inputGff = args[0];

        // GFF Reading part here...
        if (args.length == 1) {
            System.out.println("Creating Fasta file from the GFF.....");
            vgf.createFastaFromGFF(inputGff);
            return;
        }

        try {
            // GFF writing part here...
            if (args.length == 5) {
                String inputMzid = args[1];
                boolean passThreshold = Boolean.valueOf(args[2]);
                String outputGff = args[3];
                String outputMzid = args[4];

                //ReadProteinRecordFromSummaryFile rip = new ReadProteinRecordFromSummaryFile(inputMzid, decoyIdentifier);
                //ArrayList<ProteinResults> prList = rip.readProteinResultsFromFile(fdrThreshold);
                List<ProteinResults> prList = vgf.getProteinResults(inputMzid, passThreshold);

                vgf.writingTheGFFfile(inputGff, outputGff, prList);
                vgf.writingTheMzidfile(inputMzid, outputMzid, prList);
            }
        } catch (NumberFormatException e) {
            String errMsg = "Exception in the arguments : Check the arguments";
            System.out.println(errMsg);
            throw e;
        } catch (Exception e) {
            String errMsg = "Exception caused in Mapping of peptides on GFF file";
            System.out.println(errMsg);
            throw e;
        }

    }

    private void writingTheMzidfile(String inputMzid, String outputMzid, List<ProteinResults> prList) {

        try {
            String outFile = outputMzid;

            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            //cvList = mzIdentML.getCvList();
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));
            Cv psiCV;
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
            if (analysisSoftwareList != null) {
                AnalysisSoftware analysisSoftware = new AnalysisSoftware();
                analysisSoftware.setName(this.getClass().getSimpleName());
                Date date = new Date();
                SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
                analysisSoftware.setName(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
                analysisSoftware.setId(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
                analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);
                marshaller.marshal(analysisSoftwareList, writer);
            }
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
                siList.getSpectrumIdentificationResult().add(sr);

            }

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
                pdList.getProteinAmbiguityGroup().add(pag);
            }

            // here 
            SequenceCollection sequenceCollection1 = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);

            SequenceCollection sq = new SequenceCollection();
            sq.getDBSequence().addAll(sequenceCollection1.getDBSequence());
            sq.getPeptide().addAll(sq.getPeptide());

            List<PeptideEvidence> peList = sq.getPeptideEvidence();
            // Perform Chromosome mapping
            for (int i = 0; i < prList.size(); i++) {
                ProteinResults pr = prList.get(i);
                String[] co_ords = new String[2];
                long start_map, end_map;
                co_ords = mapToGff(pr);
                // Take care of the negative strand reporting in GFF
                start_map = Long.parseLong(co_ords[0]);
                end_map = Long.parseLong(co_ords[1]);
                if (end_map < start_map) {
                    long tmp = start_map;
                    start_map = end_map;
                    end_map = tmp;
                }
                PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, pr.getPeptideEvidenceID());
                if (peList.contains(peptideEvidence)) {
                    peList.remove(peptideEvidence);
                    peptideEvidence.getUserParam().add(mzidLibUtils.makeUserParam("start_map", String.valueOf(start_map)));
                    peptideEvidence.getUserParam().add(mzidLibUtils.makeUserParam("end_map", String.valueOf(end_map)));
                    peList.add(peptideEvidence);
                }

            }
            sq.getPeptideEvidence().addAll(peList);

            if (sequenceCollection1 != null) {
                marshaller.marshal(sq, writer);
            }
            writer.write("\n");

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
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

            marshaller.marshal(pdList, writer);
            writer.write("\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");
            writer.write(marshaller.createMzIdentMLClosingTag());
            writer.close();
            System.out.println("Output written to " + outFile);
        } catch (JAXBException ex) {
            ex.printStackTrace();
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    
}
