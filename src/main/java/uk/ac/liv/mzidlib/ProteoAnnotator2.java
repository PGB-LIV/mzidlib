/*
 * Copyright 2016 Fawaz Ghali <fghali@liverpool.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package uk.ac.liv.mzidlib;

/**
 *
 * @author Fawaz Ghali <fghali@liverpool.ac.uk>
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.MzIdentMLObject;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.liv.mzidlib.util.FileHandler;
import uk.ac.liv.mzidlib.util.Utils;

/**
 *
 * @author Fawaz Ghali 04-Apr-2014
 */
public class ProteoAnnotator2 {

    // Calling the MzidLib Main class
    private MzIdentMLLib mzidLib = null;

    private String outputFolder = null;
    private String spectrum_files = null;
    private String inputGFF_A = null;
    private String inputFasta_A = null;
    private String inputPredicted = null;

    private String searchParameters = null;

    private String prefix = null;
    private List<String> fastaFilesList = new ArrayList<>();
    private String debugFile = "ProteoAnnotator.txt";
    private String performanceFile = "Performance.txt";
    private PrintWriter out = null;

//    private boolean useProteoGrouper = true;
    private String peptideThreshValue = "0.01";
    private String proteinThreshValue = "0.01";

    private boolean enableMsgf = true;

    //private boolean deleteFiles = false;
    long startTime, stopTime, elapsedTime;

    private boolean verbose = true;

    private String outputParameterFile = "";
    private SearchGUICLI searchGUICLI;
    private StringBuffer bf;
    private String[] listFiles;
    private File dir;
    private String omssaoutputFile = "";
    private String tandemoutputFile = "";
//    private String msgfoutputFile = "";
    private String decoyFasta = "";

    // Two stage search
    private boolean enableTwoStageSearch = false;
    private File oFolder;
    private String newMGFFolder = "new_mgf_second_stage_search";
    private String newSecondStageSearchFolder = "new_second_stage_search";
    private String firstStageFile;
    private Map<String, String> oldNewIds = new HashMap<>();
    private Map<String, String> oldNewLocations = new HashMap<>();

    //
    private int totalSearches = 0;

    // Constructor 
    public ProteoAnnotator2() {
        mzidLib = new MzIdentMLLib();

    }

    // Constructor 
    public ProteoAnnotator2(String inputGFF, String inputFasta, String spectrum_files, String outputFolder, String inpuPredicted, String searchParameters, String prefix, String peptideThreshValue, String proteinThreshValue, String enableMsgf, String enableTwoStageSearch) {

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
        if (enableMsgf.equals("0")) {
            this.enableMsgf = false;
        }

        if (enableTwoStageSearch != null && enableTwoStageSearch.equals("1")) {
            this.enableTwoStageSearch = true;
        }

        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("========================================");
        System.out.println("No. of availableProcessors: " + cores);
        System.out.println("========================================");

    }

    private void prepareTwoStageSearch() {
        // Check first stage results
        if (outputFolder == null) {
            throw new RuntimeException("First stage search results are missing");
        }

        File f = new File(outputFolder);
        // Unidentified spectra at 1% peptide level
        for (String string : f.list()) {
            if (string.endsWith("combined_fdr_peptide_threshold.mzid")) {
                firstStageFile = string;
                System.out.println("First stage search results: " + firstStageFile);
                break;
            }

        }
        if (firstStageFile == null) {
            throw new RuntimeException("Input file for second stage search is missing");
        } else {

            try {
                firstStageFile = outputFolder + File.separator + firstStageFile;
                MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(firstStageFile));
                System.out.println("Parsing first stage search results: " + firstStageFile);
                Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
                List spectraDataList = inputs.getSpectraData();
                Map<String, List> spectrumIdHashMap = new HashMap<>();
                Iterator<MzIdentMLObject> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
                while (iterSpectrumIdentificationResult.hasNext()) {
                    SpectrumIdentificationResult spectrumIdentificationResult = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.next();
                    String spectrumId = spectrumIdentificationResult.getSpectrumID();
                    String spectrumIndex = spectrumId.substring(6);
                    Integer index1 = Integer.valueOf(spectrumIndex) + 1;

                    String spectraDataRef = spectrumIdentificationResult.getSpectraDataRef();
                    List spectrums = spectrumIdHashMap.get(spectraDataRef);
                    if (spectrums == null) {
                        spectrums = new ArrayList<>();
                        spectrumIdHashMap.put(spectraDataRef, spectrums);
                    }
                    spectrums.add(index1);
                }

                System.out.println("Parsing first stage search results completed ");

                System.out.println("Creating new MGF directory");

                File nMGF = new File(outputFolder + File.separator + newMGFFolder);

                if (!nMGF.exists()) {
                    System.out.println("");
                    System.out.println("Creating directory: " + outputFolder + File.separator + newMGFFolder);
                    System.out.println("");
                    boolean result = nMGF.mkdirs();

                    if (!result) {
                        throw new RuntimeException("Creating the output folder has failed");
                    }
                } else {
                    // check output folder if empty
                    if (nMGF.list().length > 0) {
                        throw new RuntimeException("The output folder is not empty.");
                    }

                }

                for (int i = 0; i < spectraDataList.size(); i++) {
                    try {
                        SpectraData spectraData = (SpectraData) spectraDataList.get(i);
                        String spectraDataID = spectraData.getId();
                        String spectraDataLocation = spectraData.getLocation();
                        PrintWriter out12 = null;
                        Path p1 = Paths.get(spectraDataLocation);
                        String file1 = p1.getFileName().toString();
                        String newMGFFile = nMGF + File.separator + file1;
                        System.out.println("Creating new MGF file: " + newMGFFile);

                        out12 = new PrintWriter(new BufferedWriter(new FileWriter(newMGFFile, true)));

                        MgfFile mgfFile = new MgfFile(new File(spectraDataLocation));
                        List specturms = spectrumIdHashMap.get(spectraDataID);
                        for (int j = 1; j <= mgfFile.getSpectraCount(); j++) {
                            if (!specturms.contains(j)) {
                                // Write to a new MGF file
                                Spectrum specturm = mgfFile.getSpectrumByIndex(j);

                                Map map = specturm.getPeakList();
                                String s = null;
                                if (map == null) {
                                    System.out.println("Spectrum with index: " + j + " has a null PeakList");
                                } else {
                                    s = specturm.toString();
                                    s = s.replace("TITLE=", "TITLE=index=" + j + "@" + spectraDataLocation + ";");
                                    out12.println(s);
                                }

                            }
                        }
                        out12.close();

                    } catch (JMzReaderException ex) {
                        throw new RuntimeException(ex.getMessage());
                    }

                }

            } catch (IOException ex) {
                throw new RuntimeException(ex.getMessage());
            }
        }

        // Create new MGFs based with unidentified spectra at 1% peptide level
    }

    private void runSingleMgf(File mgfFileOrLocation, String string1) {
        try {

            String newOutput = outputFolder + File.separator + totalSearches;
            File f = new File(newOutput);
            boolean b = f.mkdir();

            if (!b) {
                throw new RuntimeException("Creating the output folder has failed");
            }
            String[] searchParamters;
            String inputFilePath = mgfFileOrLocation + File.separator + string1;
            File handledMgfFile = FileHandler.handleFile(inputFilePath, true, true);
            if (enableMsgf) {
                searchParamters = new String[]{"-spectrum_files", handledMgfFile.getAbsolutePath(), "-output_folder", newOutput, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "1", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
            } else {
                searchParamters = new String[]{"-spectrum_files", handledMgfFile.getAbsolutePath(), "-output_folder", newOutput, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
            }
            searchGUICLI.runSearchGUICLI(searchParamters);
            File outputDir = new File(newOutput);
            String[] outputListFiles = outputDir.list();
            String omssaFileName = null, tandemFileName = null, msgfFileName = null;
            for (String string2 : outputListFiles) {
                if (string2.endsWith(".omx")) {
                    omssaFileName = newOutput + File.separator + string2;
                }
                if (string2.endsWith(".t.xml")) {
                    tandemFileName = newOutput + File.separator + string2;
                }
                if (string2.endsWith(".msgf.mzid")) {
                    msgfFileName = newOutput + File.separator + string2;
                }

            }
            // Convert tandem search output file to mzid file
            
            tandemoutputFile = tandemFileName.substring(0, tandemFileName.lastIndexOf(".")) + "_tandem.mzid";
            String[] tandemInput = {"Tandem2mzid", tandemFileName, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
            mzidLib.init(tandemInput);
            
            // Convert omssa search output file to mzid file
            
            omssaoutputFile = omssaFileName.substring(0, omssaFileName.lastIndexOf(".")) + "_omssa.mzid";
            String[] omssaInput = {"Omssa2mzid", omssaFileName, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
            mzidLib.init(omssaInput);
            
            // Combine search engine results
            
            String combineSearchEnginesOutputFile = newOutput + File.separator + prefix + "_combined_search.mzid";
            String combineSearchEnginesOutputdebugFile = newOutput + File.separator + prefix + "combined_search_debug.txt";
            String[] combineSearchEnginesInput;
            if (enableMsgf) {
                combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-thirdFile", msgfFileName, "-thirdcvTerm", "MS:1002053", "-thirdbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
            } else {
                combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
            }
            mzidLib.init(combineSearchEnginesInput);
            
            
            String falseDiscoveryRateGlobalOutputFile = newOutput + File.separator + prefix + "combined_fdr_peptide.mzid";
            String[] falseDiscoveryRateGlobalInput = {"FalseDiscoveryRateGlobal", combineSearchEnginesOutputFile, falseDiscoveryRateGlobalOutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002356", "-betterScoresAreLower", "true", "-fdrLevel", "Peptide", "-proteinLevel", "PAG", "-compress", "false"};
            mzidLib.init(falseDiscoveryRateGlobalInput);
            String thresholdOutputFile = newOutput + File.separator + prefix + "combined_fdr_peptide_threshold.mzid";
            String[] thresholdInput = {"Threshold", falseDiscoveryRateGlobalOutputFile, thresholdOutputFile, "-isPSMThreshold", "true", "-cvAccessionForScoreThreshold", "MS:1002360", "-threshValue", "0.02", "-betterScoresAreLower", "true", "-deleteUnderThreshold", "true", "-compress", "false"};
            mzidLib.init(thresholdInput);
        } catch (Exception ex) {

            System.out.println(ex.toString());
            //ex.printStackTrace();
            //throw new RuntimeException(ex.getMessage());

        }

    }

    // Run ProteoAnnotator
    public void runProteoAnnotator() {
        Date date = new Date();
        String accession_regex = "\\S+";

        try {
            if (enableTwoStageSearch) {
                System.out.println("Prepare Two Stage Search");
                prepareTwoStageSearch();
                System.out.println("Prepare Two Stage Search completed.");
                spectrum_files = outputFolder + File.separator + newMGFFolder;
                outputFolder = outputFolder + File.separator + newSecondStageSearchFolder;
            }
            oFolder = new File(outputFolder);

            if (!oFolder.exists()) {
                System.out.println("");
                System.out.println("Creating directory: " + outputFolder);
                System.out.println("");
                boolean result = oFolder.mkdirs();

                if (!result) {
                    throw new RuntimeException("Creating the output folder has failed");
                }
            } else {
                // check output folder if empty
                if (oFolder.list().length > 0) {
                    throw new RuntimeException("The output folder is not empty.");
                }

            }
            debugFile = prefix + debugFile;
            performanceFile = prefix + performanceFile;
            debugFile = outputFolder + File.separator + debugFile;
            new File(debugFile).delete();

            out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));

            performanceFile = outputFolder + File.separator + performanceFile;
            new File(performanceFile).delete();

            PrintWriter out1 = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile, true)));
            bf = new StringBuffer();

            out.println("ProteoAnnotator");
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
            if (enableTwoStageSearch) {
                outputGenericFastaFile = prefix + "_C_generic.fasta";
            } else {
                outputGenericFastaFile = prefix + "_A_generic.fasta";
            }

            outputGenericFastaFile = outputFolder + File.separator + outputGenericFastaFile;
            if (inputFasta_A == null) {
                inputFasta_A = "";
            }
            if (inputGFF_A == null) {
                inputGFF_A = "";
            }
            String[] genericFastaInput = {"GenericFasta", inputFasta_A, outputGenericFastaFile, "-accession_regex", accession_regex, "-inputGff", inputGFF_A, "-compress", "false"};
            mzidLib.init(genericFastaInput);
            fastaFilesList.add(outputGenericFastaFile);

            if (inputPredicted != null && !inputPredicted.equals("")) {
                out.println("Handling non-canonical gene files");
                String[] predictedSets = inputPredicted.split("##");
                String[] altModels = {"B", "C", "D", "E", "F"};

                int countAlt = 0;
                if (enableTwoStageSearch) {
                    countAlt = 2;
                }
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

                        outputGenericFastaFile = prefix + altModels[countAlt] + "_generic.fasta";
                        countAlt = countAlt + 1;
                        outputGenericFastaFile = outputFolder + File.separator + outputGenericFastaFile;
                        if (inpuPredictedFasta == null) {
                            inpuPredictedFasta = "";
                        }
                        if (inpuPredictedGff == null) {
                            inpuPredictedGff = "";
                        }
                        String[] genericFastaInput1 = {"GenericFasta", inpuPredictedFasta, outputGenericFastaFile, "-accession_regex", accession_regex, "-inputGff", inpuPredictedGff, "-compress", "false"};
                        mzidLib.init(genericFastaInput1);
                        fastaFilesList.add(outputGenericFastaFile);
                    }
                }
            }

            // Combine previous fasta files
            out.println("Combine fasta files");

            dir = new File(outputFolder);
            String fastaFiles = "";
            listFiles = dir.list();

            for (int i = 0; i < fastaFilesList.size(); i++) {
                String object = fastaFilesList.get(i);
                fastaFiles = fastaFiles + object + ";";
            }
            String[] combineFastaFilesInput;
            String outputCombinedFastaFile = outputFolder + File.separator + prefix + "combined.fasta";
            if (enableTwoStageSearch) {
                combineFastaFilesInput = new String[]{"CombineFastaFiles", fastaFiles, outputCombinedFastaFile, "-enableTwoStageSearch", "1", "-compress", "false"};
            } else {
                combineFastaFilesInput = new String[]{"CombineFastaFiles", fastaFiles, outputCombinedFastaFile, "-compress", "false"};
            }

            mzidLib.init(combineFastaFilesInput);
            out.close();
            // Create a decoy database
            //out.println("Create a decoy database");
            //out.println("");
            String[] createDecoyDBInput = {"-in", outputCombinedFastaFile, "-decoy"};
            searchGUICLI = new SearchGUICLI(outputFolder, debugFile);
            decoyFasta = outputCombinedFastaFile.substring(0, outputCombinedFastaFile.lastIndexOf(".")) + "_concatenated_target_decoy.fasta";
            searchGUICLI.runDeocyCLI(outputCombinedFastaFile);

            // Prepare the search identification parameter file
            out.println("Prepare the search identification parameter file");
            out.println("");
            outputParameterFile = outputFolder + File.separator + prefix + "combined.parameters";

            // Run SearchGUI search
            searchGUICLI.runParameterFileCLI(decoyFasta, outputParameterFile, searchParameters);

            File dirMGF = new File(spectrum_files);
            String[] listMGFFiles = dirMGF.list();

            for (String string : listMGFFiles) {
                if (string.toLowerCase().endsWith(".mgf")) {
                    String newMGFLocation = outputFolder + File.separator + "mgf";
                    File n1MGF = new File(newMGFLocation);

                    if (!n1MGF.exists()) {
                        System.out.println("");
                        System.out.println("Creating directory: " + newMGFLocation);
                        System.out.println("");
                        boolean result = n1MGF.mkdirs();

                        if (!result) {
                            throw new RuntimeException("Creating the output folder has failed");
                        }
                    }
                    File mgfFileOrLocation = Utils.splitMGFsOrReturnSame(newMGFLocation, new File(spectrum_files + File.separator + string), (int) Math.pow(1024, 3), 25000);

                    if (mgfFileOrLocation.isDirectory()) {
                        String[] listMGFFiles1 = mgfFileOrLocation.list();
                        for (String string1 : listMGFFiles1) {
                            if (string1.toLowerCase().endsWith(".mgf")) {
                                System.out.println("Running search from: " + mgfFileOrLocation.getName() + " for file " + string1);
                                totalSearches = totalSearches + 1;
                                runSingleMgf(mgfFileOrLocation, string1);

                            }
                        }
                    } else {
                        System.out.println("Running search from: " + mgfFileOrLocation.getParent() + " for file " + string);

                        totalSearches = totalSearches + 1;
                        runSingleMgf(new File(mgfFileOrLocation.getParent()), string);

                    }
                }
            }

            // combine mzid files
            List<String> files2Combine = new ArrayList();
            System.out.println("Files to compile:");
            for (int i = 1; i <= totalSearches; i++) {
                String newOutput = outputFolder + File.separator + i;
                System.out.println("From Dir: " + newOutput);
                File tempFile = new File(newOutput);
                String[] filesList = tempFile.list();
                for (String string : filesList) {
                    if (string.toLowerCase().endsWith("combined_fdr_peptide_threshold.mzid")) {
                        String toAdd = newOutput + File.separator + string;
                        files2Combine.add(toAdd);
                        System.out.println("File: " + toAdd);

                    }
                }
            }

            String combinedPSMfiles = outputFolder + File.separator + prefix + "_combined.mzid";
            String combinedInput = "";
            for (int i = 0; i < files2Combine.size(); i++) {
                String string = files2Combine.get(i);
                File testSize = new File(string);
                if (testSize.length() != 0) {
                    combinedInput = combinedInput + string + ";";
                }
            }
            if (!combinedInput.equals("")) {
                combinedInput = combinedInput.substring(0, combinedInput.length() - 1);
            }
            if (files2Combine.size() == 1) {
                combinedPSMfiles = combinedInput;
            }
            System.out.println("All files to compile: " + combinedInput);
            String[] combinedPSMfilesParams = {"CombinePSMMzidFiles", combinedInput, combinedPSMfiles, "-combineFractions", "true", "-compress", "false"};
            mzidLib.init(combinedPSMfilesParams);

            String falseDiscoveryRateGlobalOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide.mzid";
            String[] falseDiscoveryRateGlobalInput = {"FalseDiscoveryRateGlobal", combinedPSMfiles, falseDiscoveryRateGlobalOutputFile, "-decoyValue", "1", "-decoyRegex", "REVERSED", "-cvTerm", "MS:1002356", "-betterScoresAreLower", "true", "-fdrLevel", "Peptide", "-proteinLevel", "PAG", "-compress", "false"};
            mzidLib.init(falseDiscoveryRateGlobalInput);

            String thresholdOutputFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold.mzid";
            String[] thresholdInput = {"Threshold", falseDiscoveryRateGlobalOutputFile, thresholdOutputFile, "-isPSMThreshold", "true", "-cvAccessionForScoreThreshold", "MS:1002360", "-threshValue", peptideThreshValue, "-betterScoresAreLower", "true", "-deleteUnderThreshold", "true", "-compress", "false"};
            mzidLib.init(thresholdInput);

            if (enableTwoStageSearch) {
                //firstStageFile;
                //thresholdOutputFile;
                //combineBothFiles
                String tempFile = thresholdOutputFile;
                try {
                    MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(tempFile));
                    Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
                    List spectraDataList = inputs.getSpectraData();
                    HashMap<String, List> spectrumIdHashMap = new HashMap<>();
                    Iterator<MzIdentMLObject> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
                    while (iterSpectrumIdentificationResult.hasNext()) {
                        SpectrumIdentificationResult spectrumIdentificationResult = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.next();
                        String spectrumId = spectrumIdentificationResult.getSpectrumID();
                        String spectrumIndex = spectrumId.substring(6);
                        Integer index1 = Integer.valueOf(spectrumIndex) + 1;

                        String spectraDataRef = spectrumIdentificationResult.getSpectraDataRef();
                        List spectrums = spectrumIdHashMap.get(spectraDataRef);
                        if (spectrums == null) {
                            spectrums = new ArrayList<>();
                            spectrumIdHashMap.put(spectraDataRef, spectrums);
                        }
                        spectrums.add(index1);
                    }

                    for (int i = 0; i < spectraDataList.size(); i++) {

                        SpectraData spectraData = (SpectraData) spectraDataList.get(i);
                        String spectraDataID = spectraData.getId();
                        String spectraDataLocation = spectraData.getLocation();

                        MgfFile mgfFile = new MgfFile(new File(spectraDataLocation));
                        List specturms = spectrumIdHashMap.get(spectraDataID);
                        //Fawaz: fix bug spectrums object is null
                        if (specturms != null && specturms.size() > 0) {
                            for (int j = 1; j <= specturms.size(); j++) {
                                Spectrum spectrum = mgfFile.getSpectrumByIndex(j);
                                String spectrumString = spectrum.toString();

                                String result = spectrumString.substring(spectrumString.lastIndexOf("TITLE=") + 6, spectrumString.indexOf(";"));
                                String oldID = result.split("@")[0];
                                String oldLocation = result.split("@")[1];
                                String newID = "index=" + String.valueOf(j);
                                oldNewIds.put(newID, oldID);
                                oldNewLocations.put(spectraDataLocation, oldLocation);

                            }
                        }

                    }

                    String content = FileUtils.readFileToString(new File(tempFile), "UTF-8");
                    for (Map.Entry<String, String> entry : oldNewIds.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        content = content.replaceAll(key, value);
                        System.out.println("Key: " + key + " value: " + value);

                    }

                    for (Map.Entry<String, String> entry : oldNewLocations.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        content = content.replaceAll(key, value);
                        System.out.println("Key: " + key + " value: " + value);

                    }

                    String tempOutName = outputFolder + File.separator + prefix + "OutputFile.mzid";
                    Writer out = new BufferedWriter(new OutputStreamWriter(
                            new FileOutputStream(tempOutName)));
                    try {
                        out.write(content);
                    } finally {
                        out.close();
                    }

                    String inputCombine = tempOutName + ";" + firstStageFile;
                    String outputCombine = outputFolder + File.separator + prefix + "OutputFile2.mzid";

                    startTime = System.currentTimeMillis();
                    String[] stage2CompinedInput = {"CombinePSMMzidFiles", inputCombine, outputCombine, "-combineFractions", "true", "-compress", "false"};
                    mzidLib.init(stage2CompinedInput);
                    stopTime = System.currentTimeMillis();
                    elapsedTime = stopTime - startTime;
                    if (verbose) {
                        bf.append("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                        System.out.println("\n\nCombinePSMMzidFiles time " + elapsedTime / 1000 + " Seconds");
                    }

                    thresholdOutputFile = outputCombine;

                } catch (Exception ex) {
                    throw new RuntimeException(ex.getMessage());
                }

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
            startTime = System.currentTimeMillis();
            mzidLib.init(proteoGrouperFDRThresholdInput);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nThreshold time " + elapsedTime / 1000 + " Seconds");
            }
            out.println("Handling non-A");
            //  Threshold nonAScore < 0.00001 - delete under Threshold
            out.println("Threshold nonAScore < 0.00001 - delete under Threshold");
            String proteoGrouperNonAThresholdFile = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff2_proteoGrouper_nonA_Threshold.mzid";
            String[] proteoGrouperNonAThresholdInput = {"Threshold", proteoGrouperOutputFile, proteoGrouperNonAThresholdFile, "-isPSMThreshold", "false", "-cvAccessionForScoreThreshold", "MS:1002474", "-threshValue", "0.00001", "-betterScoresAreLower", "false", "-deleteUnderThreshold", "true", "-scoreLevel", "PAG", "-compress", "false"};

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

            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput4);
            stopTime = System.currentTimeMillis();
            elapsedTime = stopTime - startTime;
            if (verbose) {
                bf.append("\n\nMzid2Csv time " + elapsedTime / 1000 + " Seconds");
            }
            // Add exportPeptides option
            String mzid2CsvOutputFile5 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportPeptides.csv";
            String[] mzid2CsvInput5 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile5, "-exportType", "exportPeptides", "-compress", "false"};

            startTime = System.currentTimeMillis();
            mzidLib.init(mzid2CsvInput5);
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

}
