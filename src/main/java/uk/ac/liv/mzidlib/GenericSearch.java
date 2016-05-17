package uk.ac.liv.mzidlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

/**
 *
 * @author Fawaz Ghali 08-Sep-2014
 */
public class GenericSearch {

    // Calling the MzidLib Main class
    private MzIdentMLLib mzidLib = null;

    private String outputFolder = null;
    private String spectrum_files = null;
    private String inputFasta_A = null;
    private String searchParameters = null;

    private String prefix = null;
    private String debugFile = "GenericSearch.txt";
    private PrintWriter out = null;

    private String peptideThreshValue = "0.01";
    private String proteinThreshValue = "0.01";

    private boolean enablePercolator = false;

    // Constructor 
    public GenericSearch() {
        mzidLib = new MzIdentMLLib();

    }

    // Constructor 
    public GenericSearch(String inputFasta, String spectrum_files, String outputFolder, String searchParameters, String prefix, String peptideThreshValue, String proteinThreshValue) {

        mzidLib = new MzIdentMLLib();
        this.outputFolder = outputFolder;
        this.spectrum_files = spectrum_files;
        this.inputFasta_A = inputFasta;
        this.searchParameters = searchParameters;
        this.prefix = prefix;
        this.peptideThreshValue = peptideThreshValue;
        this.proteinThreshValue = proteinThreshValue;
    }

    // Run the main function
    public static void main(String[] args) {

    }

    // Run ProteoAnnotator
    public void runGenericSearch() {
        Date date = new Date();
        try {
            File oFolder = new File(outputFolder);

            if (!oFolder.exists()) {
                System.out.println("");
                System.out.println("Creating directory: " + outputFolder);
                System.out.println("");
                boolean result = oFolder.mkdirs();

                if (!result) {
                    throw new RuntimeException("Creating the output folder has failed");
                }
            }
            debugFile = outputFolder + File.separator + debugFile;
            new File(debugFile).delete();

            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));
            out.println("GenericSearch");
            out.println("");

            out.println(date.toString());

            out.println("");
            out.println("Output Folder: " + outputFolder);
            out.println("");

            String outputGenericFastaFile = null;

            if (spectrum_files == null || spectrum_files.equals("")) {
                throw new RuntimeException("The MGF file is missing.");
            }

            if (prefix == null) {
                prefix = "";
            }

            // Check if the GFF contains FASTA and create a generic fasta file from the GFF.
            out.println("Check if the canonical GFF contains FASTA and create a generic fasta file from the GFF.");
            out.println("The input FASTA file will be ignored if the GFF contains FASTA.");
            if (inputFasta_A != null && !inputFasta_A.equals("")) {
                Path p = Paths.get(inputFasta_A);
                String file = p.getFileName().toString();
                outputGenericFastaFile = prefix + file.substring(0, file.lastIndexOf(".")) + "_generic.fasta";
                outputGenericFastaFile = outputFolder + File.separator + outputGenericFastaFile;
                String accession_regex = "\\S+";

                String[] genericFastaInput = {"GenericFasta", inputFasta_A, outputGenericFastaFile, "-accession_regex", accession_regex, "-inputGff", "", "-compress", "false"};
                print(out, genericFastaInput);
                mzidLib.init(genericFastaInput);

            }

            out.close();
            File dir = new File(outputFolder);

            String[] listFiles = dir.list();
            // Create a decoy database
            //out.println("Create a decoy database");
            //out.println("");
            String[] createDecoyDBInput = {"-in", outputGenericFastaFile, "-decoy"};
            SearchGUICLI SearchGUICLI = new SearchGUICLI(outputFolder, debugFile);
            String decoyFasta = outputGenericFastaFile.substring(0, outputGenericFastaFile.lastIndexOf(".")) + "_concatenated_target_decoy.fasta";
            SearchGUICLI.runDeocyCLI(outputGenericFastaFile);

            // Prepare the search identification parameter file
            out.println("Prepare the search identification parameter file");
            out.println("");
            String outputParameterFile = outputFolder + File.separator + prefix + "combined.parameters";

            // Run SearchGUI search
            SearchGUICLI.runParameterFileCLI(decoyFasta, outputParameterFile, searchParameters);

            // Test if spectrum_files is a single MGF file or a folder containing multiple MGF files
            boolean testDir = new File(spectrum_files).isDirectory();
            String omssaoutputFile = "";
            String tandemoutputFile = "";
            if (testDir) {
                File dirMGF = new File(spectrum_files);

                String[] listMGFFiles = dirMGF.list();
                // Loop on all MGF files and for each MGF file call SearchGUI
                for (String string : listMGFFiles) {
                    if (string.endsWith(".mgf")) {
                        out.println("Run SearchGUI search on: " + spectrum_files + File.separator + string);
                        out.println("");
                        // Check MGF file size and spectra count
                        MgfFile dataFile = new MgfFile(new File(spectrum_files + File.separator + string));
                        long fileSize = (new File(spectrum_files + File.separator + string)).length();
                        if (dataFile.getSpectraCount() > 25000 || fileSize > Math.pow(1024, 3)) {
                            throw new RuntimeException("The MGF file is bigger than 1GB or the Spectra Count > 25000.");

                        }
                        String[] searchParamters = {"-spectrum_files", spectrum_files + File.separator + string, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0", "-tide", "0"};
                        SearchGUICLI.runSearchGUICLI(searchParamters);
                        out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));
                    }

                }

                List<String> omssaFiles = new ArrayList();
                List<String> tandemFiles = new ArrayList();

                List<String> tandemFDRThresholdoutputFiles = new ArrayList();
                List<String> omssaFDRThresholdoutputFiles = new ArrayList();

                listFiles = dir.list();
                // loop on output folder and group omssa files together and tandem files together
                for (String string : listFiles) {
                    if (string.endsWith(".omx")) {
                        omssaFiles.add(outputFolder + File.separator + string);
                    }
                    if (string.endsWith(".t.xml")) {
                        tandemFiles.add(outputFolder + File.separator + string);
                    }

                }
                // for every tandem file, convert it to mzid, call fdr and threshold, then combine all mzid files
                for (int i = 0; i < tandemFiles.size(); i++) {
                    String tandemFile = tandemFiles.get(i);
                    tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                    if (enablePercolator) {
                        String[] tandemInput = {"XtandemPercolator", tandemFile, outputFolder, "-decoyRegex", "REVERSED", "-proteinCodeRegex", "\\S+"};
                        print(out, tandemInput);
                        mzidLib.init(tandemInput);
                        File tandem_file = new File(tandemFile);
                        String name = tandem_file.getName().substring(0, tandem_file.getName().lastIndexOf("."));
                        String mzidfile = outputFolder + name + ".mzid";
                        File oldfile = new File(mzidfile);
                        File newfile = new File(tandemoutputFile);

                        if (oldfile.renameTo(newfile)) {
                            System.out.println("Rename succesful");
                        } else {
                            System.out.println("Rename failed");
                        }
                        tandemFDRThresholdoutputFiles.add(tandemoutputFile);
                    } else {
                        String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
                        print(out, tandemInput);
                        mzidLib.init(tandemInput);
                        tandemFDRThresholdoutputFiles.add(tandemoutputFile);
                    }
                }

                // for every omssa file, convert it to mzid, call fdr and threshold, then combine all mzid files
                for (int i = 0; i < omssaFiles.size(); i++) {
                    String omssaFile = omssaFiles.get(i);
                    omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "_omssa.mzid";
                    if (enablePercolator) {
                        String[] omssaInput = {"OmssaPercolator", omssaFile, outputFolder, "-database", decoyFasta, "-decoyRegex", "REVERSED"};
                        print(out, omssaInput);
                        mzidLib.init(omssaInput);
                        File omssa_file = new File(omssaFile);
                        String name = omssa_file.getName().substring(0, omssa_file.getName().lastIndexOf("."));
                        String mzidfile = outputFolder + name + ".mzid";
                        File oldfile = new File(mzidfile);
                        File newfile = new File(omssaoutputFile);

                        if (oldfile.renameTo(newfile)) {
                            System.out.println("Rename succesful");
                        } else {
                            System.out.println("Rename failed");
                        }
                        omssaFDRThresholdoutputFiles.add(omssaoutputFile);
                    } else {
                        String[] omssaInput = {"Omssa2mzid", omssaFile, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
                        print(out, omssaInput);
                        mzidLib.init(omssaInput);
                        omssaFDRThresholdoutputFiles.add(omssaoutputFile);
                    }

                }

                // combine results for omssa and tandem
                omssaoutputFile = outputFolder + File.separator + prefix + "combined_omssa.mzid";
                String omssaInput = "";
                for (int i = 0; i < omssaFDRThresholdoutputFiles.size(); i++) {
                    String string = omssaFDRThresholdoutputFiles.get(i);
                    omssaInput = omssaInput + string + ";";
                }
                if (!omssaInput.equals("")) {
                    omssaInput = omssaInput.substring(0, omssaInput.length() - 1);
                }
                if (omssaFDRThresholdoutputFiles.size() == 1) {
                    omssaoutputFile = omssaInput;
                }
                String[] omssaCompinedInput = {"CombinePSMMzidFiles", omssaInput, omssaoutputFile, "-combineFractions", "true", "-compress", "false"};
                print(out, omssaCompinedInput);
                mzidLib.init(omssaCompinedInput);

                tandemoutputFile = outputFolder + File.separator + prefix + "combined_tandem.mzid";
                String tandemInput = "";
                for (int i = 0; i < tandemFDRThresholdoutputFiles.size(); i++) {
                    String string = tandemFDRThresholdoutputFiles.get(i);
                    tandemInput = tandemInput + string + ";";
                }
                if (!tandemInput.equals("")) {
                    tandemInput = tandemInput.substring(0, tandemInput.length() - 1);
                }

                if (tandemFDRThresholdoutputFiles.size() == 1) {
                    tandemoutputFile = tandemInput;
                }
                String[] tandemCompinedInput = {"CombinePSMMzidFiles", tandemInput, tandemoutputFile, "-combineFractions", "true", "-compress", "false"};
                print(out, tandemCompinedInput);
                mzidLib.init(tandemCompinedInput);

            } else {
                out.println("Run SearchGUI search.");
                out.println("");
                // Check MGF file size and spectra count
                MgfFile dataFile = new MgfFile(new File(spectrum_files));
                long fileSize = (new File(spectrum_files)).length();
                if (dataFile.getSpectraCount() > 25000 || fileSize > Math.pow(1024, 3)) {
                    throw new RuntimeException("The MGF file is bigger than 1GB or the Spectra Count > 25000.");

                }
                if (dataFile.getSpectraCount() < 1000) {
                    throw new RuntimeException("The MGF file is too small to run.");
                }

                String[] searchParamters = {"-spectrum_files", spectrum_files, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0", "-tide", "0"};
                SearchGUICLI.runSearchGUICLI(searchParamters);
                out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));

                // Convert Omssa and Xtadem to Mzid files
                out.println("Convert Omssa and Xtadem to Mzid files.");

                String omssaFile = null;
                String tandemFile = null;
                listFiles = dir.list();
                for (String string : listFiles) {
                    if (string.endsWith(".omx")) {
                        omssaFile = outputFolder + File.separator + string;
                    }
                    if (string.endsWith(".t.xml")) {
                        tandemFile = outputFolder + File.separator + string;
                    }

                }

                if (omssaFile == null || tandemFile == null) {
                    throw new RuntimeException("Omssa or XTandem files do not exist");
                }
                omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "_omssa.mzid";
                if (enablePercolator) {

                    String[] omssaInput = {"OmssaPercolator", omssaFile, outputFolder, "-database", decoyFasta, "-decoyRegex", "REVERSED"};
                    print(out, omssaInput);
                    mzidLib.init(omssaInput);
                    File omssa_file = new File(omssaFile);
                    String name = omssa_file.getName().substring(0, omssa_file.getName().lastIndexOf("."));
                    String mzidfile = outputFolder + name + ".mzid";
                    File oldfile = new File(mzidfile);
                    File newfile = new File(omssaoutputFile);

                    if (oldfile.renameTo(newfile)) {
                        System.out.println("Rename succesful: " + newfile);
                    } else {
                        System.out.println("Rename failed: " + newfile);
                    }

                } else {
                    String[] omssaInput = {"Omssa2mzid", omssaFile, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
                    print(out, omssaInput);
                    mzidLib.init(omssaInput);
                }

                if (enablePercolator) {
                    String[] tandemInput = {"XtandemPercolator", tandemFile, outputFolder, "-decoyRegex", "REVERSED", "-proteinCodeRegex", "\\S+"};
                    print(out, tandemInput);
                    mzidLib.init(tandemInput);
                    File tandem_file = new File(tandemFile);
                    String name = tandem_file.getName().substring(0, tandem_file.getName().lastIndexOf("."));
                    String mzidfile = outputFolder + name + ".mzid";
                    File oldfile = new File(mzidfile);
                    File newfile = new File(tandemoutputFile);

                    if (oldfile.renameTo(newfile)) {
                        System.out.println("Rename succesful");
                    } else {
                        System.out.println("Rename failed");
                    }

                } else {
                    tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                    String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
                    print(out, tandemInput);
                    mzidLib.init(tandemInput);
                }

            }

            // Combine search engines
            out.println("Combine search engines.");

            String combineSearchEnginesOutputFile = outputFolder + File.separator + prefix + "combined_search.mzid";

            String combineSearchEnginesOutputdebugFile = outputFolder + File.separator + prefix + "combined_search_debug.txt";

            String[] combineSearchEnginesInput = {"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
            print(out, combineSearchEnginesInput);
            mzidLib.init(combineSearchEnginesInput);
            // Global FDR Peptide FDR < 0.01
            out.println("Global FDR Peptide FDR < 0.01");

            String falseDiscoveryRateGlobalOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide.mzid";
            String[] falseDiscoveryRateGlobalInput = {"FalseDiscoveryRateGlobal", combineSearchEnginesOutputFile, falseDiscoveryRateGlobalOutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002356", "-betterScoresAreLower", "true", "-fdrLevel", "Peptide", "-proteinLevel", "PAG", "-compress", "false"};
            print(out, falseDiscoveryRateGlobalInput);
            mzidLib.init(falseDiscoveryRateGlobalInput);

            // Threshold 1% Peptide FDR - delete under Threshold
            out.println("Threshold 1% Peptide FDR - delete under Threshold");
            String thresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold.mzid";
            String[] thresholdInput = {"Threshold", falseDiscoveryRateGlobalOutputFile, thresholdOutputFile, "-isPSMThreshold", "true", "-cvAccessionForScoreThreshold", "MS:1002360", "-threshValue", peptideThreshValue, "-betterScoresAreLower", "true", "-deleteUnderThreshold", "true", "-compress", "false"};
            print(out, thresholdInput);
            mzidLib.init(thresholdInput);

            out.println("Run ProteoGrouper");

            String proteoGrouperOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_proteoGrouper.mzid";
            String[] proteoGrouperInput = {"ProteoGrouper", thresholdOutputFile, proteoGrouperOutputFile, "-requireSIIsToPassThreshold", "true", "-verboseOutput", "false", "-cvAccForSIIScore", "MS:1002360", "-logTransScore", "true", "-version1_1", "true", "-useProteoAnnotator", "true", "-compress", "false"};
            print(out, proteoGrouperInput);
            mzidLib.init(proteoGrouperInput);

            //  FDR Global protein level
            out.println("FDR Global protein level");

            String proteoGrouperFDROutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_proteoGrouper_fdr.mzid";
            String[] proteoGrouperFDRInput = {"FalseDiscoveryRateGlobal", proteoGrouperOutputFile, proteoGrouperFDROutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002236", "-betterScoresAreLower", "false", "-fdrLevel", "ProteinGroup", "-proteinLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRInput);
            mzidLib.init(proteoGrouperFDRInput);

            //  Threshold 1% Protein FDR - delete under Threshold
            out.println("Threshold 1% Protein FDR - delete under Threshold");

            String proteoGrouperFDRThresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_proteoGrouper_fdr_Threshold.mzid";
            String[] proteoGrouperFDRThresholdInput = {"Threshold", proteoGrouperFDROutputFile, proteoGrouperFDRThresholdOutputFile, "-isPSMThreshold", "false", "-cvAccessionForScoreThreshold", "MS:1002373", "-threshValue", proteinThreshValue, "-betterScoresAreLower", "true", "-deleteUnderThreshold", "false", "-scoreLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRThresholdInput);
            mzidLib.init(proteoGrouperFDRThresholdInput);

            //  AddEmpaiToMzid 
//            out.println("AddEmpaiToMzid");
//
//            String proteoGrouperFDRThresholdAddEmpaiToMzidOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_fdr_Threshold_AddEmpaiToMzid.mzid";
//            String[] proteoGrouperFDRThresholdInputAddEmpaiToMzid = {"AddEmpaiToMzid", proteoGrouperFDRThresholdOutputFile, proteoGrouperFDRThresholdAddEmpaiToMzidOutputFile, "-fastaFile", decoyFasta, "-accessionSplitRegex", "/ /", "-compress", "false"};
//            print(out, proteoGrouperFDRThresholdInputAddEmpaiToMzid);
//            mzidLib.init(proteoGrouperFDRThresholdInputAddEmpaiToMzid);
            // Simon: Commented out the following as they are not currently relevant to the GenericSearch.
            // exportPSMs to csv using Mzid2Csv
            out.println("exportPSMs to csv using Mzid2Csv");

            String mzid2CsvOutputFile1 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportPSMs.csv";
            String[] mzid2CsvInput1 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile1, "-exportType", "exportPSMs", "-compress", "false"};
            print(out, mzid2CsvInput1);
            mzidLib.init(mzid2CsvInput1);

            // exportProteinGroups to csv using Mzid2Csv
            out.println("exportProteinGroups to csv using Mzid2Csv");

            String mzid2CsvOutputFile2 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportProteinGroups.csv";
            String[] mzid2CsvInput2 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile2, "-exportType", "exportProteinGroups", "-compress", "false"};
            print(out, mzid2CsvInput2);
            mzidLib.init(mzid2CsvInput2);

            // exportRepProteinPerPAGOnly to csv using Mzid2Csv
            out.println("exportRepProteinPerPAGOnly to csv using Mzid2Csv");

            String mzid2CsvOutputFile3 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportRepProteinPerPAGOnly.csv";
            String[] mzid2CsvInput3 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile3, "-exportType", "exportRepProteinPerPAGOnly", "-compress", "false"};
            print(out, mzid2CsvInput3);
            mzidLib.init(mzid2CsvInput3);

            // exportProteinsOnly to csv using Mzid2Csv
            out.println("exportProteinsOnly to csv using Mzid2Csv");

            String mzid2CsvOutputFile4 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportProteinsOnly.csv";
            String[] mzid2CsvInput4 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile4, "-exportType", "exportProteinsOnly", "-compress", "false"};
            print(out, mzid2CsvInput4);
            mzidLib.init(mzid2CsvInput4);
            out.println("");
            out.println("");
            out.println("GenericSearch finished on " + date.toString());
            out.close();

        } catch (Exception ex) {
            ex.printStackTrace();
            out.println("");
            out.println(ex.getMessage());
            out.println("");
            out.println("");
            out.println("ProteoAnnotator finished on " + date.toString());
            out.close();
        }

    }

    private void print(PrintWriter out, String[] fastaGffMapperInput) {
        out.println("");
        out.println("\tAbout to run " + fastaGffMapperInput[0] + " using the following parameters:");
        String s = "";
        for (int i = 1; i < fastaGffMapperInput.length; i++) {
            String string = fastaGffMapperInput[i];
            s = s + string + " ";
        }
        out.println("\t" + s);
        out.println("");
        out.println("");
    }

}
