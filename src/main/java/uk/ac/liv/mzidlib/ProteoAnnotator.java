package uk.ac.liv.mzidlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.liv.mzidlib.Callable.MzidCallable;

/**
 *
 * @author Fawaz Ghali 04-Apr-2014
 */
public class ProteoAnnotator {

    // Calling the MzidLib Main class
    private MzIdentMLLib mzidLib = null;

    private String outputFolder = null;
    private String spectrum_files = null;
    private String inputGFF_A = null;
    private String inputFasta_A = null;
    private String inputPredicted = null;

    private String searchParameters = null;

    private String prefix = null;
    private List<String> fastaFilesList = new ArrayList();
    private String debugFile = "ProteoAnnotator.txt";
    private String performanceFile = "Performance.txt";
    private PrintWriter out = null;

    private boolean useProteoGrouper = true;
    private String peptideThreshValue = "0.01";
    private String proteinThreshValue = "0.01";

    private boolean enablePercolator = false;
    private boolean enableMsgf = false;

    //private boolean deleteFiles = false;
    long startTime, stopTime, elapsedTime;

    private boolean verbose = true;

    //private ExecutorService executorPool;
    //List<MzidCallable> collection;
    //int cores = 1;
    // Constructor 
    public ProteoAnnotator() {
        mzidLib = new MzIdentMLLib();

    }

    // Constructor 
    public ProteoAnnotator(String inputGFF, String inputFasta, String spectrum_files, String outputFolder, String inpuPredicted, String searchParameters, String prefix, String peptideThreshValue, String proteinThreshValue) {

        mzidLib = new MzIdentMLLib();
        this.outputFolder = outputFolder;
        this.spectrum_files = spectrum_files;
        this.inputGFF_A = inputGFF;
        this.inputFasta_A = inputFasta;
        this.inputPredicted = inpuPredicted;
        this.searchParameters = searchParameters;
        this.prefix = prefix;
        this.peptideThreshValue = peptideThreshValue;
        this.proteinThreshValue = proteinThreshValue;
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("========================================");
        System.out.println("No. of availableProcessors: " + cores);
        System.out.println("========================================");
        //collection = new ArrayList<MzidCallable>();

    }

    // Run the main function
    public static void main(String[] args) {

    }

    // Run ProteoAnnotator
    public void runProteoAnnotator() {
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
            debugFile = prefix + debugFile;
            performanceFile = prefix + performanceFile;
            debugFile = outputFolder + File.separator + debugFile;
            new File(debugFile).delete();

            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));

            performanceFile = outputFolder + File.separator + performanceFile;
            new File(performanceFile).delete();

            PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile, true)));
            StringBuffer bf = new StringBuffer();

            out.println("ProteoAnnotator");
            out.println("");

            out.println(date.toString());

            out.println("");
            out.println("Output Folder: " + outputFolder);
            out.println("");

            String outputGenericFastaFile = null;

//            if (inputGFF_A == null || inputGFF_A.equals("")) {
//                throw new RuntimeException("The GFF file is missing.");
//            } else {
//                out.println("The canonical GFF file: " + inputGFF_A);
//                out.println("");
//
//            }
            if (spectrum_files == null || spectrum_files.equals("")) {
                throw new RuntimeException("The MGF file is missing.");
            }

            if (prefix == null) {
                prefix = "";
            }

            // Check if the GFF contains FASTA and create a generic fasta file from the GFF.
            out.println("Check if the canonical GFF contains FASTA and create a generic fasta file from the GFF.");
            out.println("The input FASTA file will be ignored if the GFF contains FASTA.");
//            if (inputGFF_A != null && !inputGFF_A.equals("")) {
//                Path p = Paths.get(inputGFF_A);
//                String file = p.getFileName().toString();
            outputGenericFastaFile = prefix + "_A_generic.fasta";
            outputGenericFastaFile = outputFolder + File.separator + outputGenericFastaFile;
            String accession_regex = "\\S+";
            if (inputFasta_A == null) {
                inputFasta_A = "";
            }
            if (inputGFF_A == null) {
                inputGFF_A = "";
            }
            String[] genericFastaInput = {"GenericFasta", inputFasta_A, outputGenericFastaFile, "-accession_regex", accession_regex, "-inputGff", inputGFF_A, "-compress", "false"};
            print(out, genericFastaInput);
            mzidLib.init(genericFastaInput);
            fastaFilesList.add(outputGenericFastaFile);
//            }

            if (inputPredicted != null && !inputPredicted.equals("")) {
                out.println("Handling non-canonical gene files");
                String[] predictedSets = inputPredicted.split("##");
                String[] altModels = {"B", "C", "D", "E", "F"};
                int countAlt = 0;
                for (String predictedSet : predictedSets) {
                    String[] pairs = predictedSet.split(";");
                    String inpuPredictedGff = null;
                    String inpuPredictedFasta = null;
                    if (pairs != null) {

                        if (pairs.length == 2) {
                            inpuPredictedFasta = pairs[1];
                            inpuPredictedGff = pairs[0];
                        } else if (pairs.length == 1) {
                            if (pairs[0].endsWith("fa") || pairs[0].endsWith("fasta")) {
                                inpuPredictedFasta = pairs[0];
                            }
                            if (pairs[0].endsWith("gff") || pairs[0].endsWith("gff3")) {
                                inpuPredictedGff = pairs[0];
                            }
                        }

                        // Check if the GFF contains FASTA and create a generic fasta file from the GFF.
//                        if (inpuPredictedGff != null && !inpuPredictedGff.equals("")) {
//                            Path p = Paths.get(inpuPredictedGff);
//                            String file = p.getFileName().toString();
                        outputGenericFastaFile = prefix + altModels[countAlt] + "_generic.fasta";
                        countAlt = countAlt + 1;
                        outputGenericFastaFile = outputFolder + File.separator + outputGenericFastaFile;
                        String accession_regex1 = "\\S+";
                        if (inpuPredictedFasta == null) {
                            inpuPredictedFasta = "";
                        }
                         if (inpuPredictedGff == null) {
                            inpuPredictedGff = "";
                        }
                        String[] genericFastaInput1 = {"GenericFasta", inpuPredictedFasta, outputGenericFastaFile, "-accession_regex", accession_regex1, "-inputGff", inpuPredictedGff, "-compress", "false"};
                        print(out, genericFastaInput1);
                        mzidLib.init(genericFastaInput1);
                        fastaFilesList.add(outputGenericFastaFile);
//                        }

                    }
                }
            }

            // Combine previous fasta files
            out.println("Combine fasta files");

            File dir = new File(outputFolder);
            String fastaFiles = "";
            String[] listFiles = dir.list();

            for (int i = 0; i < fastaFilesList.size(); i++) {
                String object = fastaFilesList.get(i);
                fastaFiles = fastaFiles + object + ";";
            }

            String outputCombinedFastaFile = outputFolder + File.separator + prefix + "combined.fasta";
            String[] combineFastaFilesInput = {"CombineFastaFiles", fastaFiles, outputCombinedFastaFile, "-compress", "false"};
            print(out, combineFastaFilesInput);
            mzidLib.init(combineFastaFilesInput);
            out.close();
            // Create a decoy database
            //out.println("Create a decoy database");
            //out.println("");
            String[] createDecoyDBInput = {"-in", outputCombinedFastaFile, "-decoy"};
            SearchGUICLI SearchGUICLI = new SearchGUICLI(outputFolder, debugFile);
            String decoyFasta = outputCombinedFastaFile.substring(0, outputCombinedFastaFile.lastIndexOf(".")) + "_concatenated_target_decoy.fasta";
            SearchGUICLI.runDeocyCLI(outputCombinedFastaFile);

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
            String msgfoutputFile = "";
            if (testDir) {
                File dirMGF = new File(spectrum_files);

                String[] listMGFFiles = dirMGF.list();
                // Loop on all MGF files and for each MGF file call SearchGUI
                startTime = System.currentTimeMillis();
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

                        if (dataFile.getSpectraCount() < 1000) {
                            throw new RuntimeException("The MGF file is too small to run.");
                        }
                        String[] searchParamters;
                        if (enableMsgf) {
                            searchParamters = new String[]{"-spectrum_files", spectrum_files + File.separator + string, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "1", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                        } else {

                            searchParamters = new String[]{"-spectrum_files", spectrum_files + File.separator + string, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                        }
                        SearchGUICLI.runSearchGUICLI(searchParamters);
                        out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));
                    }

                }
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nSearch time " + elapsedTime / 1000 + " Seconds");
                }

                List<String> omssaFiles = new ArrayList();
                List<String> tandemFiles = new ArrayList();
                List<String> msgfFiles = new ArrayList();

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

                    if (string.endsWith(".msgf.mzid")) {
                        msgfFiles.add(outputFolder + File.separator + string);
                    }

                }
                // for every tandem file, convert it to mzid, call fdr and threshold, then combine all mzid files
                startTime = System.currentTimeMillis();
                for (int i = 0; i < tandemFiles.size(); i++) {
                    String tandemFile = tandemFiles.get(i);

                    if (enablePercolator) {
                        tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "AddP.mzid";
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
                        tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                        String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
                        print(out, tandemInput);
                        tandemFDRThresholdoutputFiles.add(tandemoutputFile);
                        mzidLib.init(tandemInput);
                    }
                    //MzidCallable task = new MzidCallable(mzidLib, tandemInput);
                    //collection.add(task);

                }
//                try {
//                    executorPool = Executors.newFixedThreadPool(cores);
//                    executorPool.invokeAll(collection);
//                } catch (Exception e) {
//                    e.printStackTrace();
//                } finally {
//                    executorPool.shutdown();
//                }

                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nTandem2mzid time " + elapsedTime / 1000 + " Seconds");
                    System.out.println("\n\nTandem2mzid time " + elapsedTime / 1000 + " Seconds");
                }

                // for every omssa file, convert it to mzid, call fdr and threshold, then combine all mzid files
                startTime = System.currentTimeMillis();
                // collection.clear();
                for (int i = 0; i < omssaFiles.size(); i++) {
                    String omssaFile = omssaFiles.get(i);

                    if (enablePercolator) {
                        omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "AddP.mzid";
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
                        omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "_omssa.mzid";
                        String[] omssaInput = {"Omssa2mzid", omssaFile, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
                        print(out, omssaInput);
                        omssaFDRThresholdoutputFiles.add(omssaoutputFile);
                        mzidLib.init(omssaInput);
                    }
                    //MzidCallable task = new MzidCallable(mzidLib, omssaInput);
                    //collection.add(task);

                }
//                try {
//                    executorPool = Executors.newFixedThreadPool(cores);
//                    executorPool.invokeAll(collection);
//                } catch (Exception e) {
//                    e.printStackTrace();
//                } finally {
//                    executorPool.shutdown();
//                }
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nOmssa2mzid time " + elapsedTime / 1000 + " Seconds");
                    System.out.println("\n\nOmssa2mzid time " + elapsedTime / 1000 + " Seconds");
                }

                // combine results for oms sa and tandem
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
                startTime = System.currentTimeMillis();
                startTime = System.currentTimeMillis();
                String[] omssaCompinedInput = {"CombinePSMMzidFiles", omssaInput, omssaoutputFile, "-combineFractions", "true", "-compress", "false"};
                print(out, omssaCompinedInput);
                mzidLib.init(omssaCompinedInput);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                    System.out.println("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                }
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
                startTime = System.currentTimeMillis();
                String[] tandemCompinedInput = {"CombinePSMMzidFiles", tandemInput, tandemoutputFile, "-combineFractions", "true", "-compress", "false"};
                print(out, tandemCompinedInput);
                mzidLib.init(tandemCompinedInput);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                    System.out.println("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                }
                if (enableMsgf) {
                    msgfoutputFile = outputFolder + File.separator + prefix + "combined_msgf.mzid";
                    String msgfInput = "";
//                    if (enablePercolator) {
//                        for (int i = 0; i < msgfFiles.size(); i++) {
//                            String stringInputMSGF = msgfFiles.get(i);
//                            String stringOutputMSGF = stringInputMSGF.substring(0, stringInputMSGF.lastIndexOf(".")) + "AddP.mzid";
//                            msgfInput = msgfInput + stringOutputMSGF + ";";
//                            String[] msgfInputs = {"MsgfPercolator", stringInputMSGF, outputFolder, "-decoyRegex", "REVERSED"};
//                            print(out, msgfInputs);
//                            mzidLib.init(msgfInputs);
//                            
//                            
//                            
//                        }
//
//                    } else {

                    for (int i = 0; i < msgfFiles.size(); i++) {
                        String string = msgfFiles.get(i);
                        msgfInput = msgfInput + string + ";";
                    }
                    if (!msgfInput.equals("")) {
                        msgfInput = msgfInput.substring(0, msgfInput.length() - 1);
                    }

                    if (msgfFiles.size() == 1) {
                        msgfoutputFile = msgfInput;
//                        }
                    }

                    startTime = System.currentTimeMillis();
                    String[] msgfCompinedInput = {"CombinePSMMzidFiles", msgfInput, msgfoutputFile, "-combineFractions", "true", "-compress", "false"};
                    print(out, msgfCompinedInput);
                    mzidLib.init(msgfCompinedInput);
                    stopTime = System.currentTimeMillis();
                    elapsedTime = stopTime - startTime;
                    if (verbose) {
                        bf.append("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                        System.out.println("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                    }
                }

            } else {
                out.println("Run SearchGUI search.");
                out.println("");
                // Check MGF file size and spectra count
                MgfFile dataFile = new MgfFile(new File(spectrum_files));
                long fileSize = (new File(spectrum_files)).length();
                if (dataFile.getSpectraCount() > 25000 || fileSize > Math.pow(1024, 3)) {
                    throw new RuntimeException("The MGF file is bigger than 1GB or the Spectra Count > 25000.");

                }
                startTime = System.currentTimeMillis();
                String[] searchParamters = {"-spectrum_files", spectrum_files, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                SearchGUICLI.runSearchGUICLI(searchParamters);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nSearch time " + elapsedTime / 1000 + " Seconds");
                }
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
                startTime = System.currentTimeMillis();
                omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "_omssa.mzid";
                String[] omssaInput = {"Omssa2mzid", omssaFile, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
                print(out, omssaInput);
                mzidLib.init(omssaInput);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nOmssa2mzid time " + elapsedTime / 1000 + " Seconds");
                }
                tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
                print(out, tandemInput);
                startTime = System.currentTimeMillis();
                mzidLib.init(tandemInput);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nTandem2mzid time " + elapsedTime / 1000 + " Seconds");
                }

            }
            // Combine search engines
            out.println("Combine search engines.");

            String combineSearchEnginesOutputFile = outputFolder + File.separator + prefix + "combined_search.mzid";

            String combineSearchEnginesOutputdebugFile = outputFolder + File.separator + prefix + "combined_search_debug.txt";
            String[] combineSearchEnginesInput;

            if (enablePercolator) {
                // User MS:1001493 percolator:PEP
                if (enableMsgf) {
                    combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001493", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001493", "-secondbetterScoresAreLower", "true", "-thirdFile", msgfoutputFile, "-thirdcvTerm", "MS:1002053", "-thirdbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
                } else {
                    combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001493", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001493", "-secondbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
                }
            } else {
                if (enableMsgf) {
                    combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-thirdFile", msgfoutputFile, "-thirdcvTerm", "MS:1002053", "-thirdbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
                } else {
                    combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
                }
            }

            print(out, combineSearchEnginesInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(combineSearchEnginesInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nCombineSearchEngines time " + elapsedTime / 1000 + " Seconds");
                System.out.println("\n\nCombineSearchEngines time " + elapsedTime / 1000 + " Seconds");
            }
            // Global FDR Peptide FDR < 0.01
            out.println("Global FDR Peptide FDR < 0.01");

            String falseDiscoveryRateGlobalOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide.mzid";
            String[] falseDiscoveryRateGlobalInput = {"FalseDiscoveryRateGlobal", combineSearchEnginesOutputFile, falseDiscoveryRateGlobalOutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002356", "-betterScoresAreLower", "true", "-fdrLevel", "Peptide", "-proteinLevel", "PAG", "-compress", "false"};
            print(out, falseDiscoveryRateGlobalInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(falseDiscoveryRateGlobalInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nFalseDiscoveryRateGlobal time " + elapsedTime / 1000 + " Seconds");
            }
            // Threshold 1% Peptide FDR - delete under Threshold
            out.println("Threshold 1% Peptide FDR - delete under Threshold");
            String thresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold.mzid";
            String[] thresholdInput = {"Threshold", falseDiscoveryRateGlobalOutputFile, thresholdOutputFile, "-isPSMThreshold", "true", "-cvAccessionForScoreThreshold", "MS:1002360", "-threshValue", peptideThreshValue, "-betterScoresAreLower", "true", "-deleteUnderThreshold", "true", "-compress", "false"};
            print(out, thresholdInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(thresholdInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nThreshold time " + elapsedTime / 1000 + " Seconds");
            }

            String addGenomeCoordinatesForPeptidesOutputFile1 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff.mzid";
            if (inputGFF_A != null && !inputGFF_A.equals("")) {
                // AddGenomeCoordinatesForPeptides
                out.println("Add genome Coordinates for peptides");

                Path p1 = Paths.get(inputGFF_A);
                String file1 = p1.getFileName().toString();
                String gffOutputFile1 = prefix + file1.substring(0, file1.lastIndexOf(".")) + "_annotated.gff";
                gffOutputFile1 = outputFolder + File.separator + gffOutputFile1;
                String[] addGenomeCoordinatesForPeptidesInput = {"AddGenomeCoordinatesForPeptides", thresholdOutputFile, addGenomeCoordinatesForPeptidesOutputFile1, "-inputGff", inputGFF_A, "-outputGff", gffOutputFile1, "-compress", "false"};
                print(out, addGenomeCoordinatesForPeptidesInput);
                startTime = System.currentTimeMillis();
                mzidLib.init(addGenomeCoordinatesForPeptidesInput);
                stopTime = System.currentTimeMillis();
                elapsedTime = stopTime - startTime;
                if (verbose) {
                    bf.append("\n\nAddGenomeCoordinatesForPeptides time " + elapsedTime / 1000 + " Seconds");
                }
                // AddGenomeCoordinatesForPeptides
                out.println("Add genome Coordinates for peptides");
            } else {
                addGenomeCoordinatesForPeptidesOutputFile1 = thresholdOutputFile;
            }

            String addGenomeCoordinatesForPeptidesOutputFile = null;
            if (inputPredicted != null && !inputPredicted.equals("")) {
                String[] predictedSets = inputPredicted.split("##");
                int i = 2;
                for (String predictedSet : predictedSets) {
                    String[] pairs = predictedSet.split(";");
                    if ((pairs != null && pairs.length == 2) || (pairs != null && pairs.length == 1 && pairs[0].endsWith("gff")) || (pairs != null && pairs.length == 1 && pairs[0].endsWith("gff3"))) {
                        addGenomeCoordinatesForPeptidesOutputFile = outputFolder + File.separator + "combined_fdr_peptide_threshold_mappedGff_" + i + "_.mzid";
                        i = i + 1;
                        Path p2 = Paths.get(pairs[0]);
                        String file2 = p2.getFileName().toString();
                        String gffOutputFile2 = prefix + file2.substring(0, file2.lastIndexOf(".")) + "_annotated.gff";
                        gffOutputFile2 = outputFolder + File.separator + gffOutputFile2;
                        String[] addGenomeCoordinatesForPeptidesInputPredicted = {"AddGenomeCoordinatesForPeptides", addGenomeCoordinatesForPeptidesOutputFile1, addGenomeCoordinatesForPeptidesOutputFile, "-inputGff", pairs[0], "-outputGff", gffOutputFile2, "-compress", "false"};
                        print(out, addGenomeCoordinatesForPeptidesInputPredicted);
                        startTime = System.currentTimeMillis();
                        mzidLib.init(addGenomeCoordinatesForPeptidesInputPredicted);
                        stopTime = System.currentTimeMillis();
                        elapsedTime = stopTime - startTime;
                        if (verbose) {
                            bf.append("\n\nAddGenomeCoordinatesForPeptides time " + elapsedTime / 1000 + " Seconds");
                        }
                    }
                }
            }
            if (addGenomeCoordinatesForPeptidesOutputFile == null) {
                addGenomeCoordinatesForPeptidesOutputFile = addGenomeCoordinatesForPeptidesOutputFile1;
            }

            //  ProteoGrouper 
            out.println("Run ProteoGrouper");

            String proteoGrouperOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper.mzid";
            String[] proteoGrouperInput = {"ProteoGrouper", addGenomeCoordinatesForPeptidesOutputFile, proteoGrouperOutputFile, "-requireSIIsToPassThreshold", "true", "-verboseOutput", "false", "-cvAccForSIIScore", "MS:1002360", "-logTransScore", "true", "-version1_1", "true", "-useProteoAnnotator", "true", "-compress", "false"};
            print(out, proteoGrouperInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nProteoGrouper time " + elapsedTime / 1000 + " Seconds");
            }
            //  FDR Global protein level
            out.println("FDR Global protein level");

            String proteoGrouperFDROutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_fdr.mzid";
            String[] proteoGrouperFDRInput = {"FalseDiscoveryRateGlobal", proteoGrouperOutputFile, proteoGrouperFDROutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002236", "-betterScoresAreLower", "false", "-fdrLevel", "ProteinGroup", "-proteinLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperFDRInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nFalseDiscoveryRateGlobal time " + elapsedTime / 1000 + " Seconds");
            }
            //  Threshold 1% Protein FDR - delete under Threshold
            out.println("Threshold 1% Protein FDR - delete under Threshold");

            String proteoGrouperFDRThresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_fdr_Threshold.mzid";
            String[] proteoGrouperFDRThresholdInput = {"Threshold", proteoGrouperFDROutputFile, proteoGrouperFDRThresholdOutputFile, "-isPSMThreshold", "false", "-cvAccessionForScoreThreshold", "MS:1002373", "-threshValue", proteinThreshValue, "-betterScoresAreLower", "true", "-deleteUnderThreshold", "false", "-scoreLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRThresholdInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperFDRThresholdInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nThreshold time " + elapsedTime / 1000 + " Seconds");
            }
//            //  AddEmpaiToMzid 
//            
//             out.println("AddEmpaiToMzid");
//
//            String proteoGrouperFDRThresholdAddEmpaiToMzidOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_fdr_Threshold_AddEmpaiToMzid.mzid";
//            String[] proteoGrouperFDRThresholdInputAddEmpaiToMzid = {"AddEmpaiToMzid", proteoGrouperFDRThresholdOutputFile, proteoGrouperFDRThresholdAddEmpaiToMzidOutputFile, "-fastaFile", decoyFasta, "-accessionSplitRegex", "/ /" , "-compress", "false"};
//            print(out, proteoGrouperFDRThresholdInputAddEmpaiToMzid);
//            mzidLib.init(proteoGrouperFDRThresholdInputAddEmpaiToMzid);
            out.println("Handling non-A");
            //  Threshold nonAScore < 0.00001 - delete under Threshold
            out.println("Threshold nonAScore < 0.00001 - delete under Threshold");
            String proteoGrouperNonAThresholdFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_nonA_Threshold.mzid";
            String[] proteoGrouperNonAThresholdInput = {"Threshold", proteoGrouperOutputFile, proteoGrouperNonAThresholdFile, "-isPSMThreshold", "false", "-cvAccessionForScoreThreshold", "MS:1002474", "-threshValue", "0.00001", "-betterScoresAreLower", "false", "-deleteUnderThreshold", "true", "-scoreLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperNonAThresholdInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperNonAThresholdInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nThreshold time " + elapsedTime / 1000 + " Seconds");
            }
            //  FDR Global protein level - non-A score
            out.println("FDR Global protein level - non-A score");

            String proteoGrouperFDRNonAOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_nonA_Threshold_FDR.mzid";
            String[] proteoGrouperFDRNonAInput = {"FalseDiscoveryRateGlobal", proteoGrouperNonAThresholdFile, proteoGrouperFDRNonAOutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002474", "-betterScoresAreLower", "false", "-fdrLevel", "ProteinGroup", "-proteinLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRNonAInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperFDRNonAInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nFalseDiscoveryRateGlobal time " + elapsedTime / 1000 + " Seconds");
            }
            //  Threshold 5% Protein FDR - delete under Threshold
            out.println(" Threshold 5% Protein FDR - delete under Threshold");

            String proteoGrouperFDRNonAThresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_nonA_Threshold_FDR_Threshold.mzid";
            String[] proteoGrouperFDRNonAThresholdInput = {"Threshold", proteoGrouperFDRNonAOutputFile, proteoGrouperFDRNonAThresholdOutputFile, "-isPSMThreshold", "false", "-cvAccessionForScoreThreshold", "MS:1002373", "-threshValue", "0.05", "-betterScoresAreLower", "true", "-deleteUnderThreshold", "false", "-scoreLevel", "PAG", "-compress", "false"};
            print(out, proteoGrouperFDRNonAThresholdInput);
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperFDRNonAThresholdInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nThreshold time " + elapsedTime / 1000 + " Seconds");
            }
            // exportProteoAnnotator to csv using Mzid2Csv
            out.println("exportProteoAnnotator to csv using Mzid2Csv");

            String mzid2CsvOutputFile0 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportProteoAnnotator.csv";
            String[] mzid2CsvInput0 = {"Mzid2Csv", proteoGrouperFDRNonAOutputFile, mzid2CsvOutputFile0, "-exportType", "exportProteoAnnotator", "-compress", "false"};
            print(out, mzid2CsvInput0);
            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput0);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            // exportPSMs to csv using Mzid2Csv
            out.println("exportPSMs to csv using Mzid2Csv");

            String mzid2CsvOutputFile1 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportPSMs.csv";
            String[] mzid2CsvInput1 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile1, "-exportType", "exportPSMs", "-compress", "false"};
            print(out, mzid2CsvInput1);
            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput1);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            // exportProteinGroups to csv using Mzid2Csv
            out.println("exportProteinGroups to csv using Mzid2Csv");

            String mzid2CsvOutputFile2 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportProteinGroups.csv";
            String[] mzid2CsvInput2 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile2, "-exportType", "exportProteinGroups", "-compress", "false"};
            print(out, mzid2CsvInput2);
            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput2);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            // exportRepProteinPerPAGOnly to csv using Mzid2Csv
            out.println("exportRepProteinPerPAGOnly to csv using Mzid2Csv");

            String mzid2CsvOutputFile3 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportRepProteinPerPAGOnly.csv";
            String[] mzid2CsvInput3 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile3, "-exportType", "exportRepProteinPerPAGOnly", "-compress", "false"};
            print(out, mzid2CsvInput3);
            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput3);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            // exportProteinsOnly to csv using Mzid2Csv
            out.println("exportProteinsOnly to csv using Mzid2Csv");

            String mzid2CsvOutputFile4 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportProteinsOnly.csv";
            String[] mzid2CsvInput4 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile4, "-exportType", "exportProteinsOnly", "-compress", "false"};
            print(out, mzid2CsvInput4);
            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput4);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            out.println("");
            out.println("");
            out.println("ProteoAnnotator finished on " + date.toString());
            out.close();

            bf.append("\n\n");
            bf.append("\n\n");
            bf.append("\n\nProteoAnnotator finished on " + date.toString());
            out1.println(bf.toString());
            out1.close();
            // executorPool.shutdown();
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
