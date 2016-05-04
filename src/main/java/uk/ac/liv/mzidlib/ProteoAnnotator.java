package uk.ac.liv.mzidlib;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
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
import uk.ac.liv.mzidlib.util.Utils;

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
    private String msgfoutputFile = "";
    private String decoyFasta = "";

    // Two stage search
    private boolean enableTwoStageSearch = false;
    private File oFolder;
    private String newMGFFolder = "new_mgf_second_stage_search";
    private String newSecondStageSearchFolder = "new_second_stage_search";
    private String firstStageFile;
    private HashMap<String, String> oldNewIds = new HashMap<String, String>();
    private HashMap<String, String> oldNewLocations = new HashMap<String, String>();

    // Constructor 
    public ProteoAnnotator() {
        mzidLib = new MzIdentMLLib();

    }

    // Constructor 
    public ProteoAnnotator(String inputGFF, String inputFasta, String spectrum_files, String outputFolder, String inpuPredicted, String searchParameters, String prefix, String peptideThreshValue, String proteinThreshValue, String enableMsgf, String enableTwoStageSearch) {

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
                HashMap<String, List> spectrumIdHashMap = new HashMap<String, List>();
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
                    spectrums.add(index1.intValue());
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
                                    s = s.replace("TITLE=", "TITLE=" + j + "@" + spectraDataLocation + ";");
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

    // Run ProteoAnnotator
    public void runProteoAnnotator() {
        Date date = new Date();

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

            print(out, combineFastaFilesInput);
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
            // Test if spectrum_files is a single MGF file or a folder containing multiple MGF files
            boolean testDir = new File(spectrum_files).isDirectory();

            runMultipleSearches();
            // Combine search engines
            out.println("Combine search engines.");

            String combineSearchEnginesOutputFile = outputFolder + File.separator + prefix + "combined_search.mzid";

            String combineSearchEnginesOutputdebugFile = outputFolder + File.separator + prefix + "combined_search_debug.txt";
            String[] combineSearchEnginesInput;

            if (enableMsgf) {
                combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-thirdFile", msgfoutputFile, "-thirdcvTerm", "MS:1002053", "-thirdbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
            } else {
                combineSearchEnginesInput = new String[]{"CombineSearchEngines", "-firstFile", omssaoutputFile, "-firstcvTerm", "MS:1001328", "-firstbetterScoresAreLower", "true", "-secondFile", tandemoutputFile, "-secondcvTerm", "MS:1001330", "-secondbetterScoresAreLower", "true", "-decoyRatio", "1", "-rank", "3", "-outputFile", combineSearchEnginesOutputFile, "-debugFile", combineSearchEnginesOutputdebugFile, "-decoyRegex", "REVERSED", "-compress", "false"};
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

            if (enableTwoStageSearch) {
                //firstStageFile;
                //thresholdOutputFile;
                //combineBothFiles
                String tempFile = thresholdOutputFile;
                try {
                    MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(tempFile));
                    Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
                    List spectraDataList = inputs.getSpectraData();
                    HashMap<String, List> spectrumIdHashMap = new HashMap<String, List>();
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
                        spectrums.add(index1.intValue());
                    }

                    for (int i = 0; i < spectraDataList.size(); i++) {

                        SpectraData spectraData = (SpectraData) spectraDataList.get(i);
                        String spectraDataID = spectraData.getId();
                        String spectraDataLocation = spectraData.getLocation();

                        MgfFile mgfFile = new MgfFile(new File(spectraDataLocation));
                        List specturms = spectrumIdHashMap.get(spectraDataID);

                        for (int j = 1; j <= specturms.size(); j++) {
                            Spectrum spectrum = mgfFile.getSpectrumByIndex(j);
                            String spectrumString = spectrum.toString();

                            String result = spectrumString.substring(spectrumString.lastIndexOf("TITLE=") + 6, spectrumString.indexOf(";"));
                            String oldID = result.split("@")[0];

                            String oldLocation = result.split("@")[1];

                            oldNewIds.put(String.valueOf(j), oldID);
                            oldNewLocations.put(spectraDataLocation, oldLocation);

                        }

                    }

                    String content = FileUtils.readFileToString(new File(tempFile), "UTF-8");
                    for (Map.Entry<String, String> entry : oldNewIds.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        content = content.replaceAll(key, value);

                    }

                    for (Map.Entry<String, String> entry : oldNewLocations.entrySet()) {
                        String key = entry.getKey();
                        String value = entry.getValue();
                        content = content.replaceAll(key, value);

                    }

                    String tempOutName = outputFolder + File.separator + prefix + "OutputFile.mzid";
                    File tempOutFile = new File(tempOutName);
                    FileUtils.writeStringToFile(tempOutFile, content, "UTF-8");

                    String inputCombine = tempOutName + ";" + firstStageFile;
                    String outputCombine = outputFolder + File.separator + prefix + "OutputFile2.mzid";

                    startTime = System.currentTimeMillis();
                    String[] stage2CompinedInput = {"CombinePSMMzidFiles", inputCombine, outputCombine, "-combineFractions", "true", "-compress", "false"};
                    print(out, stage2CompinedInput);
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
            // Add exportPeptides option
            String mzid2CsvOutputFile5 = outputFolder + File.separator + prefix + "combined_fdr_peptide_threshold_mappedGff_ProteoGrouper_exportPeptides.csv";
            String[] mzid2CsvInput5 = {"Mzid2Csv", proteoGrouperFDRThresholdOutputFile, mzid2CsvOutputFile5, "-exportType", "exportPeptides", "-compress", "false"};
            print(out, mzid2CsvInput5);
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

    private void runMultipleSearches() throws JMzReaderException, IOException, Exception {
        File dirMGF = new File(spectrum_files);

        String[] listMGFFiles = dirMGF.list();
        // Loop on all MGF files and for each MGF file call SearchGUI
        startTime = System.currentTimeMillis();
        for (String string : listMGFFiles) {
            if (string.toLowerCase().endsWith(".mgf")) {
                String newMGFLocation = outputFolder +File.separator + "mgf";
                File mgfFileOrLocation = Utils.splitMGFsOrReturnSame(newMGFLocation, new File(spectrum_files + File.separator + string), (int) Math.pow(1024, 3), 25000);
                if (mgfFileOrLocation.isDirectory()) {
                    String[] listMGFFiles1 = mgfFileOrLocation.list();
                    for (String string1 : listMGFFiles1) {
                        if (string1.toLowerCase().endsWith(".mgf")) {

                            String[] searchParamters;
                            if (enableMsgf) {
                                searchParamters = new String[]{"-spectrum_files", mgfFileOrLocation + File.separator + string1, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "1", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                            } else {

                                searchParamters = new String[]{"-spectrum_files", mgfFileOrLocation + File.separator + string1, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                            }
                            searchGUICLI.runSearchGUICLI(searchParamters);
                            out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));
                        }
                    }
                } else {
                    String[] searchParamters;
                    if (enableMsgf) {
                        searchParamters = new String[]{"-spectrum_files", spectrum_files + File.separator + string, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "1", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                    } else {
                        searchParamters = new String[]{"-spectrum_files", spectrum_files + File.separator + string, "-output_folder", outputFolder, "-id_params", outputParameterFile, "-comet", "0", "-tide", "0", "-msgf", "0", "-ms_amanda", "0", "-output_option", "3", "-myrimatch", "0"};
                    }
                    searchGUICLI.runSearchGUICLI(searchParamters);
                    out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile, true)));
                }

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

            tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
            String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};
            print(out, tandemInput);
            tandemFDRThresholdoutputFiles.add(tandemoutputFile);
            mzidLib.init(tandemInput);

        }

        stopTime = System.currentTimeMillis();
        elapsedTime = stopTime - startTime;
        if (verbose) {
            bf.append("\n\nTandem2mzid time " + elapsedTime / 1000 + " Seconds");
            System.out.println("\n\nTandem2mzid time " + elapsedTime / 1000 + " Seconds");
        }

        // for every omssa file, convert it to mzid, call fdr and threshold, then combine all mzid files
        startTime = System.currentTimeMillis();

        for (int i = 0; i < omssaFiles.size(); i++) {
            String omssaFile = omssaFiles.get(i);

            omssaoutputFile = omssaFile.substring(0, omssaFile.lastIndexOf(".")) + "_omssa.mzid";
            String[] omssaInput = {"Omssa2mzid", omssaFile, omssaoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};
            print(out, omssaInput);
            omssaFDRThresholdoutputFiles.add(omssaoutputFile);
            mzidLib.init(omssaInput);

        }
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
            for (int i = 0; i < msgfFiles.size(); i++) {
                String string = msgfFiles.get(i);
                msgfInput = msgfInput + string + ";";
            }
            if (!msgfInput.equals("")) {
                msgfInput = msgfInput.substring(0, msgfInput.length() - 1);
            }

            if (msgfFiles.size() == 1) {
                msgfoutputFile = msgfInput;
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

    }

// Run the main function
    public static void main(String[] args) {

    }
}
