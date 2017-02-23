package uk.ac.liv.mzidlib.fasta;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

import java.util.regex.Pattern;

import uk.ac.liv.mzidlib.util.FileHandler;

/**
 *
 * @author Fawaz Ghali
 */
public class GenericFasta {

    private String inputFasta_A;
    private String output;
    private String accession_regex;
    private String inputGff;

    public GenericFasta(String input, String output, String accession_regex, String inputGff) {
        this.inputFasta_A = input;
        System.out.println("Input Fasta A: " + inputFasta_A);
        this.output = output;
        System.out.println("output: " + output);
        this.accession_regex = accession_regex;
        System.out.println("accession_regex: " + accession_regex);
        this.inputGff = inputGff;
        System.out.println("Input GFF: " + inputGff);
        //this.description_regex = description_regex;
        if (inputGff != null && !inputGff.equals("")) {
            createGenericeFastaFromGff();
        } else {
            createGenericFasta();
        }

    }

    private void createGenericFasta() {
        Writer temp_out = null;
        try {

            System.out.println("Converting " + inputFasta_A + " to a generic fasta format");
            BufferedReader in = new BufferedReader(new FileReader(inputFasta_A));
            String line;
            File temp = new File(output);
            temp_out = new BufferedWriter(new FileWriter(temp));
            while ((line = in.readLine()) != null) {
                if (line.startsWith(">")) {
                    temp_out.write(genericHeader(line) + "\n");
                } else {
                    temp_out.write(line + "\n");
                }
            }
            in.close();
            temp_out.close();
            System.out.println("Done. Check: " + output);
            // return the name of the file with fasta sequences
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                temp_out.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    public String genericHeader(String header) {
        String newHeader = "";
        String accession = "";
        String description = "";
        Pattern accession_pattern = Pattern.compile(accession_regex);
        //Pattern description_pattern = Pattern.compile(description_regex);

        //Matcher accession_matcher = accession_pattern.matcher(header);
        //Matcher description_matcher = description_pattern.matcher(header);
        // \S+
        //E:\proteogalaxy\pipeline\Data\fasta\ToxoDB-9.0_TgondiiME49_AnnotatedProteins_generice.fasta E:\proteogalaxy\pipeline\Data\fasta\ToxoDB-9.0_TgondiiME49_AnnotatedProteins_generice1.fasta "\|(.*?)\|" ">..|[^|]*|\(.*\)" 
//        if (accession_matcher.find()) {
//            accession = accession_matcher.group();
//            //If there is a description
//            if (accession.length() + 3 < header.length()) {
//                description = header.substring(accession.length() + 3);
//            }
//
//            if (accession.startsWith(">sp|")) {
//                accession = accession.replaceFirst(">sp|", ">generic|");
//            } else {
//                accession = accession.replaceFirst(">", "");
//                if (accession.contains("|")) {
//                    accession = accession.substring(accession.lastIndexOf("|") + 1);
//                }
//                newHeader = ">generic|" + accession + "|" + description;
//            }
//        }
        if (header.contains(" ")) {
            accession = header.split(" ")[0];
            //If there is a description
//            if (accession.length() + 3 < header.length()) {
//                description = header.substring(accession.length() + 1);
//            }

            if (accession.startsWith(">sp|")) {
                accession = accession.replaceFirst(">sp|", ">generic");
                //newHeader = accession  + "|" + description;
                newHeader = accession;
            } else {
                accession = accession.replaceFirst(">", "");
                if (accession.contains("|")) {
                    accession = accession.substring(accession.lastIndexOf("|") + 1);
                }
                newHeader = ">generic|" + accession + "|" + description;
            }
        } else {
            accession = header;
            accession = accession.replaceFirst(">", "");
            newHeader = ">generic|" + accession + "|";
        }

        return newHeader;
    }

    private void createGenericeFastaFromGff() {
        BufferedReader in = null;
        // Variables needed during processing of ##FASTA file
        File tempFasta;
        String tempFastaPath = "";
        Writer temp_out = null;
        boolean fastaRegionFlag = false;
        try {
            File handledGffFile = FileHandler.handleFile(inputGff, true, true); 
            String handledGffFilePath = handledGffFile.getAbsolutePath();
            System.out.println("Parsing the GFF: " + handledGffFilePath);                       
            in = new BufferedReader(new FileReader(handledGffFile));
            String line;

            boolean createFastaFileFlagFromGff = false;

            while ((line = in.readLine()) != null) {
                if (line.equals("##FASTA")) {
                    fastaRegionFlag = true;
                    System.out.println("The GFF: " + handledGffFilePath + " contains FASTA");
                }
                if (fastaRegionFlag) {
                    if (!createFastaFileFlagFromGff) {
                        tempFastaPath = handledGffFilePath.substring(0, handledGffFilePath.lastIndexOf('.')) + ".fasta";
                        tempFasta = new File(tempFastaPath);
                        temp_out = new BufferedWriter(new FileWriter(tempFasta));
                        createFastaFileFlagFromGff = true;
                    }

                    // If all well, then process fasta - skip the line having ##FASTA 
                    if (line.equals("##FASTA")) {
                        continue;
                    }
                    if ((line.startsWith(">") && line.contains("cds_")) || (line.startsWith(">") && line.contains("rna_")) || (line.startsWith(">gb|"))) {
                        temp_out.write(line + "\n");
                        String seq = "";
                        while ((seq = in.readLine()) != null) {
                            if (!seq.startsWith(">")) {
                                temp_out.write(seq + "\n");
                            } else {
                                break;
                            }
                        }

                    }
                }
//                }else{
//                    System.out.println("The GFF file does not contain FASTA header.");
//                }
            }
            if (createFastaFileFlagFromGff) {
                System.out.println("Creating a temp fasta file from the GFF.");
                inputFasta_A = tempFastaPath;
                createGenericFasta();
            } else if ((inputFasta_A == null) || inputFasta_A.equals("")) {
                System.out.println("The GFF does not contain FASTA, and no FASTA file is provided.");
                throw new RuntimeException("The GFF does not contain FASTA, and no FASTA file is provided.");
            } else {
                System.out.println("Create Generic Fasta");
                createGenericFasta();
            }

            System.out.println("Parsing " + handledGffFilePath + " done.");
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                in.close();
                if (temp_out != null) {
                    temp_out.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

//    public void attachGffTag() {
//        System.out.println("Parsing " + inputFasta + " to get the accessions list...");
//        try {
//            BufferedReader in = new BufferedReader(new FileReader(inputFasta));
//            Writer temp_out = null;
//            File temp = new File(outputFasta);
//            temp_out = new BufferedWriter(new FileWriter(temp));
//            String line;
//            while ((line = in.readLine()) != null) {
//                if (line.startsWith(">")) {
//                    String[] lineArray = line.split("\\|");
//                    if (lineArray.length > 1) {
//                        String accession = lineArray[1];
//                        accession = accession.trim();
//                        String exonString = "|gff=##";
//
//                        if (line.endsWith("|")) {
//                            line = line.substring(0, line.length() - 1);
//                        }
//                        line = line + exonString;
//                        temp_out.write(line + "\n");
//
//                    } else {
//                        in.close();
//                        temp_out.close();
//                        throw new RuntimeException("No match");
//                    }
//                } else {
//                    temp_out.write(line + "\n");
//                }
//            }
//            in.close();
//            temp_out.close();
//
//        } catch (IOException ex) {
//            ex.printStackTrace();
//        }
//        System.out.println("Parsing " + inputFasta + " done.");
//
//    }
    public static void main(String[] args) {
        GenericFasta gf = new GenericFasta(args[0], args[1], args[2], args[3]);

    }
}
