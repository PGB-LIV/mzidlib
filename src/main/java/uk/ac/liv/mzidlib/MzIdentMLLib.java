//hello world
package uk.ac.liv.mzidlib;

import uk.ac.liv.mzidlib.fasta.InsertMetaDataFromFasta;
import uk.ac.liv.mzidlib.proteogrouper.ProteoGrouper;
import uk.ac.liv.mzidlib.fdr.FalseDiscoveryRate;
import uk.ac.liv.mzidlib.fdr.FalseDiscoveryRateGlobal;
import uk.ac.liv.mzidlib.converters.MzIdentMLToMzTab;
import uk.ac.liv.mzidlib.converters.MzIdentMLToCSV;
import uk.ac.liv.mzidlib.fasta.CreateRestrictedFASTADatabase;
import bgi.ipeak.percolator.MsgfPercolator;
import bgi.ipeak.percolator.OmssaPercolator;
import bgi.ipeak.percolator.XtandemPercolator;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import uk.ac.liv.mzidlib.converters.Csv2mzid;
import uk.ac.liv.mzidlib.converters.Omssa2mzid;
import uk.ac.liv.mzidlib.converters.Tandem2mzid;
import uk.ac.liv.mzidlib.multiplesearch.CombineSearchEngines;
import uk.ac.liv.mzidlib.util.Gzipper;
import uk.ac.liv.mzidlib.util.Utils;
import uk.ac.liv.mzidlib.experimental.RescoreMods;
import uk.ac.liv.mzidlib.fasta.CombineFastaFiles;
import uk.ac.liv.mzidlib.fasta.AddGenomeCoordinatesForPeptides;
import uk.ac.liv.mzidlib.fasta.GenericFasta;

/**
 *
 * @author Fawaz Ghali, University of Liverpool, 2011
 */
public class MzIdentMLLib {

    public static String fdrParams = "-decoyRegex decoyRegex -decoyValue decoyToTargetRatio -cvTerm cvTerm -betterScoresAreLower true|false [-compress true|false]";
    public static String fdrUsageExample = " -decoyRegex Rev_ -decoyValue 1 -cvTerm MS:1001172 -betterScoresAreLower true -compress true";
    public static String fdrToolDescription = "This tool can be used for mzid files in which a decoy database search has been performed to calculate three new data types for each PSM:"
            + "Local FDR, Q-value and FDRScore [PMID: 19253293] and assign these to every PSM with new CV terms. You must specify the score you wish to order by, using -cvTerm [MS:XXXX] sourced from the PSI-MS CV"
            + "and whether the scores are ordered low to high or vice versa\n. ";
    public static String fdrUsage = "FalseDiscoveryRate input.mzid output.mzid " + fdrParams + " \n\nDescription:\n" + fdrToolDescription;
    public static String fdrGlobalParams = "-decoyValue decoyToTargetRatio -decoyRegex decoyRegex -cvTerm cvTerm -betterScoresAreLower true|false -fdrLevel PSM|Peptide|ProteinGroup -proteinLevel PDH|PAG [-compress true|false]";
    public static String fdrGlobalUsageExample = " -decoyValue 0.01 -decoyRegex REVERSED -cvTerm MS:1002356 -betterScoresAreLower true -fdrLevel Peptide -proteinLevel PAG -compress true";
    public static String fdrGlobalToolDescription = "The Global FDR module calculates the FDR on one of the three levels. 1) PSM, 2) Peptide, 3) ProteinGroup. If ProteinGroup is chosen, there are two options for protein level PAG or PDH.";
    public static String fdrGlobalUsage = "FalseDiscoveryRateGlobal input.mzid output.mzid " + fdrGlobalParams + " \n\nDescription:\n" + fdrGlobalToolDescription;
    
    
    public static String omssa2mzidparams = "[-outputFragmentation true|false] -decoyRegex decoyRegex [-omssaModsFile pathToLocalOmssaModsFile] [-userModsFile pathToLocalUserModsFile] [-compress true|false]";
    public static String omssa2mzidToolDescription = "This tool converts OMSSA omx (XML) files into mzid. It has optional parameters for inserting fragment ions into mzid (much larger files). If a decoy Regex is specified, the mzid attribute isDecoy will be set correctly for peptides."
            + " No protein inference is done by this tool (no protein list produced). To make valid mzid output, OMSSA must have been run with the option \"-w include spectra and search params in search results\"."
            + " Without this option, search paramaters cannot be extracted from OMSSA. In this case, the OMSSA CSV converter should be used. ";
    public static String omssa2mzidUsage = "Omssa2mzid input.omx output.mzid " + omssa2mzidparams + " \n\nDescription:\n" + omssa2mzidToolDescription;
    public static String omssa2mzidUsageExample = " -outputFragmentation false -decoyRegex Rev_ -compress true";
    
    public static String tandem2mzidParams = "[-outputFragmentation (true|false)] [-decoyRegex decoyRegex] [-databaseFileFormatID (e.g. MS:1001348 is FASTA format) \"MS:100blah\"] [-massSpecFileFormatID (e.g. MS:1001062 is MGF) \"MS:100blah\"] [-idsStartAtZero (true for mzML searched, false otherwise) [true|false]] [-compress true|false]";
    public static String tandem2mzidToolDescription = "This tool converts X!Tandem XML results files into mzid. There are several optional parameters: whether to export fragment ions (makes bigger files), "
            + " and include a decoy regular expression to set the isDecoy attribute in mzid. Valid mzid files require several pieces of metadata that are difficult to extract from mzid files, the format of the database searched and the file format of the input spectra. "
            + " If these parameters are not set, the converter attempts to guess these based on the file extension. In X!Tandem, the numbering of spectra differs dependent upon the input spectra type - the IDs start at zero for mzML files, the IDs start at one for other spectra types e.g. MGF. "
            + "This is a command line parameter which should be set to make sure that the mzid file references the correct spectrum in the source spectrum file. ";
    public static String tandem2mzidUsage = "Tandem2mzid input_tandem.xml output.mzid " + tandem2mzidParams + " \n\nDescription:\n" + tandem2mzidToolDescription;
    public static String tandem2mzidUsageExample = "-outputFragmentation false -decoyRegex Rev_ -databaseFileFormatID MS:1001348 -massSpecFileFormatID MS:1001062 -idsStartAtZero false -compress true";
    
    
    public static String csv2mzidParams = "-paramsFile paramsFileLocation -cvAccessionForPSMOrdering (e.g. \"MS:1001328\" is OMSSA:evalue) [-applyFixedMods true|false] [-decoyRegex decoyRegex] [-compress true|false]";
    public static String csv2mzidToolDescription = "This tool is intended for converting OMSSA CSV output in mzid. Since the CSV format does not contain any search metadata, this tool requires a separate "
            + "parameters file containing the search metadata (see example under example_files). By default, applyFixedMods = true, this means that the code attempts to add fixed modifications to every peptide based on the search parameters."
            + " This mode is required since OMSSA does not report fixed mods on peptides, which must be present in mzIdentML. This behaviour can be switched off using -applyFixedMods false. "
            + " Developers can easily adapt this tool to process other types of CSV file into mzid, by altering the file csv_config_file.csv under src/resources and re-building.";
    public static String csv2mzidUsage = "Csv2mzid input_omssa.csv output.mzid " + csv2mzidParams + " \n\nDescription:\n" + csv2mzidToolDescription;
    public static String csv2mzidUsageExample = " -paramsFile example_files/toxo_omssa_params.csv -cvAccessionForPSMOrdering \"MS:1001328\" -decoyRegex Rev_ -compress true";
    
    
    public static String mzid2CsvParams = "-exportType exportProteinGroups|exportPSMs|exportProteinsOnly|exportRepProteinPerPAGOnly|exportProteoAnnotator [-verboseOutput true|false] [-compress true|false]";
    public static String mzid2CsvToolDescription = "This tool can export from an mzid file into CSV, according to one of the four types of export specified as parameters.";
    public static String mzid2CsvUsage = "Mzid2Csv input.mzid output.csv " + mzid2CsvParams + " \n\nDescription:\n" + mzid2CsvToolDescription;
    public static String mzid2CsvUsageExample = " -exportType exportPSMs -verboseOutput false -compress true ";
    
    
    public static String thresholdParams = "-isPSMThreshold true|false -cvAccessionForScoreThreshold \"MS:100blah\" -threshValue doubleValue  -betterScoresAreLower true|false -deleteUnderThreshold true|false [-compress true|false]";
    public static String thresholdToolDescription = "This tool can be used to set the passThreshold parameter for PSMs or proteins in an mzid file, to indicate high-quality identifications that will be used by another tool. "
            + "It can handle any type of score (sourced from the PSI-MS CV) and scores can be ordered low to high or vice versa."
            + "If deleteUnderThreshold is specified, PSMs and referenced proteins under the threshold will be removed from the file.";
    public static String thresholdUsage = "Threshold input.mzid output.mzid " + thresholdParams + " \n\nDescription:\n" + thresholdToolDescription;
    public static String thresholdUsageExample = " -isPSMThreshold true -cvAccessionForScoreThreshold \"MS:1001171\" -threshValue 40 -betterScoresAreLower false -deleteUnderThreshold true -compress true";
    
    
    
    public static String proteogrouperParams = "-requireSIIsToPassThreshold true|false -cvAccForSIIScore cvAccForSIIScore -logTransScore true|false -verboseOutput true|false [-version1_1 true|false] [-useProteoAnnotator true|false] [-compress true|false]";
    public static String proteogrouperToolDescription = "This tool can perform sequence-based protein inference, based on a set of PSMs. It should be parameterized with the CV accession for the PSM score used to create a protein score. "
            + "The tool also needs to know whether the score should be log transformed (true for e/p-values etc) to create a positive protein score. If version1_1 is set to false, export is to the draft mzid 1.2 specification.";
    public static String proteogrouperUsage = "ProteoGrouper input.mzid output.mzid " + proteogrouperParams + " \n\nDescription:\n" + proteogrouperToolDescription;
    public static String proteogrouperUsageExample = " -requireSIIsToPassThreshold true -verboseOutput false -cvAccForSIIScore \"MS:1001171\" -logTransScore false -version1_1 true -compress true";
    
    
    
    public static String insertMetaDataParams = "-fastaFile fastaFilelocation -accessionSplitRegex (regularExpressionTosplitAccSurroundByForwardSlashes e.g. \"/ /\") [-compress true|false]";
    public static String insertMetaDataToolDescription = "This tool can be used to extract the description line from a given FASTA file and insert into an mzid file."
            + "The tool needs a regular expression to split the accession from the description line in the FASTA file.";
    public static String insertMetaDataUsage = "InsertMetaDataFromFasta input.mzid output.mzid " + insertMetaDataParams + " \n\nDescription:\n" + insertMetaDataToolDescription;
    public static String insertMetaDataUsageExample = " -fastaFile example_files/TgondiiME49_ToxoDB-6_2.fasta -accessionSplitRegex \"/ /\" -compress true";
    

    
    public static String emPAIParams = "-fastaFile fastaFilelocation -accessionSplitRegex (regularExpressionTosplitAccSurroundByForwardSlashes e.g. \"/ /\") [-enzymeRegex (enzymeRegex defaults to trypsin - \"(?<=[KR])(?!P)\" )] [-verboseOutput true|false] [-compress true|false]";
    public static String emPAIToolDescription = "This tool applies the emPAI protocol to mzid files [PMID 15958392]. The tool will only work if protein inference has already been performed and there is a protein list in the file. "
            + "It requires the location of the fasta file searches, a regular expression (such as a space) to split the accessions from the description lines and optionally a regular expression of the enzyme used in the search (if not trypsin). ";
    public static String emPAIUsage = "AddEmpaiToMzid input.mzid output.mzid " + emPAIParams + " \n\nDescription:\n" + emPAIToolDescription;
    public static String emPAIUsageExample = " -fastaFile example_files/TgondiiME49_ToxoDB-6_2.fasta -accessionSplitRegex \"/ /\" -verboseOutput false  -compress true";
    
    
    
    public static String combinedSearchParams = "-firstFile firstFile -firstcvTerm firstcvTerm -firstbetterScoresAreLower firstbetterScoresAreLower -secondFile secondFile -secondcvTerm secondcvTerm -secondbetterScoresAreLower secondbetterScoresAreLower -thirdFile thirdFile -thirdcvTerm thirdcvTerm -thirdbetterScoresAreLower thirdbetterScoresAreLower -rank rank -decoyRatio decoyRatio -outputFile outputFile -debugFile debugFile -decoyRegex decoyRegex -compress false";
    public static String combinedSearchUsageExample = "-firstFile iprg_omssa.mzid -firstcvTerm MS:1001328 -firstbetterScoresAreLower true -secondFile iprg_mascot.mzid -secondcvTerm MS:1001172 -secondbetterScoresAreLower true -thirdFile iprg_tandem.mzid -thirdcvTerm MS:1001330 -thirdbetterScoresAreLower true -rank 1 -decoyRatio 3 -outputFile iprg_combined.mzid -debugFile debug.txt -decoyRegex RRRR -compress false";
    public static String combinedSearchDescription = "This tool can be used for combining multiple search engines and can output csv as well as mzid files. This example for three search engines but can be used for two search engines as well\n";
    public static String combinedSearchUsage = "CombineSearchEngines " + combinedSearchParams + " \n\nDescription:\n" + combinedSearchDescription;    
    
    
    public static String createRestrictedFASTADatabaseParams = "[-compress true|false]";
    public static String createRestrictedFASTADatabaseUsageExample = "-compress true";
    public static String createRestrictedFASTADatabaseToolDescription = "read all PDHs with passthreshold=true and create a new FASTA file from these, assuming that InsertMetaDataFromFasta has already been run to insert sequences and descriptions into the file";
    public static String createRestrictedFASTADatabaseUsage = "CreateRestrictedFASTADatabase input.mzid output.fasta " + createRestrictedFASTADatabaseParams + " \n\nDescription:\n" + createRestrictedFASTADatabaseToolDescription;
    
    
    public static String mzIdentMLToMzTabParams = "[-compress true|false]";
    public static String mzIdentMLToMzTabUsageExample = "-compress true";
    public static String mzIdentMLToMzTabToolDescription = "Convert mzidentml file to mztab file";
    public static String mzIdentMLToMzTabUsage = "MzIdentMLToMzTab input.mzid output.mztab " + mzIdentMLToMzTabParams + " \n\nDescription:\n" + mzIdentMLToMzTabToolDescription;

    
    public static String rescoreModsParams = "-cvAccForScoreToAdapt \"MS:100XXXX\"  -logTransformPSMScore true|false -commonModificationWeight X -mediumModificationWeight X -rareModificationWeight X -generalModificationWeight X -pairedModificationAndUnmodWeight X -multipleVariableModWeight X -compress true|false";
    public static String rescoreModsDescription = "Note: !!This tool is experimental at this stage, and hasn't been demonstrated to work effetively yet!! \n"
            + "This tool re-scores modifications identified by the search engine. It creates a new \"PTM score\" for all PSMs identifying modifications, by multiplying a user entered (e-value value) score type (low values are better), by various weighting factors depending on the modification classification.";
    public static String rescoreModsUsage = "RescoreMods input.mzid output.mzid " + rescoreModsParams + " \n\nDescription:\n" + rescoreModsDescription;
    public static String rescoreModsUsageExample = " -cvAccForScoreToAdapt MS:1002053 -commonModificationWeight 0.5 -mediumModificationWeight 2.0 -rareModificationWeight 10.0 -generalModificationWeight 10.0 -pairedModificationAndUnmodWeight 0.1 -multipleVariableModWeight 10.0  -compress false  -compress true";
    
    
    public static String combinePSMMzidFilesParams = "-combineFractions ture|false [-compress true|false]";
    public static String combinePSMMzidFilesUsageExample = "-combineFractions true -compress true";
    public static String combinePSMMzidFilesToolDescription = "Combine multiple mzid files on PSM level. Protein-level results (if any) will be removed from the final output";
    public static String combinePSMMzidFilesUsage = "CombinePSMMzidFiles inputfolder output.mzid " + createRestrictedFASTADatabaseParams + " \n\nDescription:\n" + createRestrictedFASTADatabaseToolDescription;
    
    
    
    public static String genericFastaParams = "-accession_regex accession_regex [-inputGff inputGff] [-compress true|false]";
    public static String genericFastaUsageExample = "-accession_regex \\S+ -inputGff inputGff.gff -compress true";
    public static String genericFastaToolDescription = "Create a generic Fasta file to be used as an input for SearchGUI.";
    public static String genericFastaUsage = "GenericFasta input.fasta output.fasta " + genericFastaParams + " \n\nDescription:\n" + genericFastaToolDescription;
    
    
    
    
    
    public static String addGenomeCoordinatesForPeptidesParams = "-inputMzid inputMzid -outputMzid outputMzid -inputGff inputGff -outputGff outputGff [-compress true|false]";
    public static String addGenomeCoordinatesForPeptidesUsageExample = "-inputGff input.gff -outputGff output.gff -compress true";
    public static String addGenomeCoordinatesForPeptidesToolDescription = "Add genome coordinates for peptides from the gff file to the mzid file";
    public static String addGenomeCoordinatesForPeptidesUsage = "AddGenomeCoordinatesForPeptides input.mzid output.mzid " + addGenomeCoordinatesForPeptidesParams + " \n\nDescription:\n" + addGenomeCoordinatesForPeptidesToolDescription;
    
    
    
    public static String addRetentionTimeToMzidParams = " -inputSourceFile inputSourceFile -compress true|false";
    public static String addRetentionTimeToMzidDescription = "add Retention Time to Mzid";
    public static String addRetentionTimeToMzidUsage = "AddRetentionTimeToMzid input.mzid output.mzid " + addRetentionTimeToMzidParams + " \n\nDescription:\n" + addRetentionTimeToMzidDescription;
    public static String addRetentionTimeToMzidExample = " -inputSourceFile input.mgf -compress false";
    
    
    
    
    
    
    public static String combineFastaFilesParams = " -compress true|false";
    public static String combineFastaFilesDescription = "Combine Fasta Files";
    public static String combineFastaFilesUsage = "CombineFastaFiles 1.fasta;2.fasta output.fasta " + addRetentionTimeToMzidParams + " \n\nDescription:\n" + combineFastaFilesDescription;
    public static String combineFastaFilesExample = " -compress false";

    //XtandemPercolator, OmssaPercolator, MsgfPercolator
    public static String xtandemPercolatorParams = "-decoyRegex decoyregex -compress true|false";
    public static String xtandemPercolatorDescription = "Running XtandemPercolator";
    public static String xtandemPercolatorUsage = "XtandemPercolator input.xml outdir " + xtandemPercolatorParams + " \n\nDescription:\n" + xtandemPercolatorDescription;
    public static String xtandemPercolatorExample = " -decoyRegex decoyregex -compress false";
    //public static String omssaPercolatorParams = "-mod_file mod_file -database database -decoyregex decoyregex -compress true|false";
    public static String omssaPercolatorParams = " -database database -decoyRegex decoyregex -compress true|false";
    public static String omssaPercolatorDescription = "Running OmssaPercolator";
    public static String omssaPercolatorUsage = "OmssaPercolator input.omx outdir " + omssaPercolatorParams + " \n\nDescription:\n" + omssaPercolatorDescription;
    //public static String omssaPercolatorExample = " -mod_file E:\\all\\bo\\IPeak_V1.0\\IPeak_release\\mods.xml -database E:\\all\\bo\\IPeak_V1.0\\IPeak_release\\*.fasta -decoyregex REVERSED -compress false";
    public static String omssaPercolatorExample = " -database E:\\all\\bo\\IPeak_V1.0\\IPeak_release\\*.fasta -decoyRegex REVERSED -compress false";
    public static String msgfPercolatorParams = " -decoyRegex decoyregex -compress true|false";
    public static String msgfPercolatorDescription = "Running MsgfPercolator";
    public static String msgfPercolatorUsage = "MsgfPercolator input.mzid outdir " + msgfPercolatorParams + " \n\nDescription:\n" + msgfPercolatorDescription;
    public static String msgfPercolatorExample = " -decoyRegex decoyregex -compress false";

    


    public static String userFeedback = "java -jar jar-location/mzidentml-lib.jar ";

// Added by Fawaz Ghali to automatically update the MzidLib GUI 
    private Map<String, String> allFunctions;

    /**
     * Init all functions hashmap
     */
    public MzIdentMLLib() {
        allFunctions = new HashMap<String, String>();
        allFunctions.put("FalseDiscoveryRate", fdrParams + ";@;" + fdrUsage + ";@;" + fdrUsageExample);
        allFunctions.put("FalseDiscoveryRateGlobal", fdrGlobalParams + ";@;" + fdrGlobalUsage + ";@;" + fdrGlobalUsageExample);
        allFunctions.put("Omssa2mzid", omssa2mzidparams + ";@;" + omssa2mzidUsage + ";@;" + omssa2mzidUsageExample);
        allFunctions.put("Tandem2mzid", tandem2mzidParams + ";@;" + tandem2mzidUsage + ";@;" + tandem2mzidUsageExample);
        allFunctions.put("Csv2mzid", csv2mzidParams + ";@;" + csv2mzidUsage + ";@;" + csv2mzidUsageExample);
        allFunctions.put("Mzid2Csv", mzid2CsvParams + ";@;" + mzid2CsvUsage + ";@;" + mzid2CsvUsageExample);
        allFunctions.put("Threshold", thresholdParams + ";@;" + thresholdUsage + ";@;" + thresholdUsageExample);
        allFunctions.put("ProteoGrouper", proteogrouperParams + ";@;" + proteogrouperUsage + ";@;" + proteogrouperUsageExample);
        allFunctions.put("InsertMetaDataFromFasta", insertMetaDataParams + ";@;" + insertMetaDataUsage + ";@;" + insertMetaDataUsageExample);
        allFunctions.put("AddEmpaiToMzid", emPAIParams + ";@;" + emPAIUsage + ";@;" + emPAIUsageExample);
        allFunctions.put("CombineSearchEngines", combinedSearchParams + ";@;" + combinedSearchUsage + ";@;" + combinedSearchUsageExample);
        allFunctions.put("CreateRestrictedFASTADatabase", createRestrictedFASTADatabaseParams + ";@;" + createRestrictedFASTADatabaseUsage + ";@;" + createRestrictedFASTADatabaseUsageExample);
        allFunctions.put("MzIdentMLToMzTab", mzIdentMLToMzTabParams + ";@;" + mzIdentMLToMzTabUsage + ";@;" + mzIdentMLToMzTabUsageExample);
        allFunctions.put("RescoreMods", rescoreModsParams + ";@;" + rescoreModsUsage + ";@;" + rescoreModsUsageExample);
        allFunctions.put("CombinePSMMzidFiles", combinePSMMzidFilesParams + ";@;" + combinePSMMzidFilesUsage + ";@;" + combinePSMMzidFilesUsageExample);
        allFunctions.put("GenericFasta", genericFastaParams + ";@;" + genericFastaUsage + ";@;" + genericFastaUsageExample);
        allFunctions.put("AddGenomeCoordinatesForPeptides", addGenomeCoordinatesForPeptidesParams + ";@;" + addGenomeCoordinatesForPeptidesUsage + ";@;" + addGenomeCoordinatesForPeptidesUsageExample);
        allFunctions.put("AddRetentionTimeToMzid", addRetentionTimeToMzidParams + ";@;" + addRetentionTimeToMzidUsage + ";@;" + addRetentionTimeToMzidExample);
        allFunctions.put("CombineFastaFiles", combineFastaFilesParams + ";@;" + combineFastaFilesUsage + ";@;" + combineFastaFilesExample);
        allFunctions.put("XtandemPercolator", xtandemPercolatorParams + ";@;" + xtandemPercolatorUsage + ";@;" + xtandemPercolatorExample);
        allFunctions.put("OmssaPercolator", omssaPercolatorParams + ";@;" + omssaPercolatorUsage + ";@;" + omssaPercolatorExample);
        allFunctions.put("MsgfPercolator", msgfPercolatorParams + ";@;" + msgfPercolatorUsage + ";@;" + msgfPercolatorExample);
        

    }

    /**
     * Getter for all functions to be used in MzidLib GUI
     *
     * @return all functions as hashmap
     */
    public Map<String, String> getAllFunctions() {
        return allFunctions;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        MzIdentMLLib mzidLib = new MzIdentMLLib();
        try {

            mzidLib.init(args);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public String init(String[] args) throws Exception {
        String inputFileName = "";
        String outputFileName = "";

        String guiFeedback = "";     //0 for successful execution

        if (args != null && args.length > 1) {
            System.out.println("");
            System.out.println("About to run " + args[0]);
            System.out.println("Using the following parameters:");
            for (int i = 1; i < args.length; i++) {
                String string = args[i];
                System.out.print(string + " ");
            }
            System.out.println();
        }

        if (args.length > 3) {

            inputFileName = args[1];
            outputFileName = args[2];
            // Added by FG uncompress
            boolean testDir = new File(inputFileName).isDirectory();
            if (inputFileName.equals(outputFileName) && !testDir) {
                guiFeedback = "Error: input and output file names are same";

            } else {
                boolean uncompress = false;
                File inputFile;
                if (inputFileName.endsWith(".gz")) {
                    inputFile = Gzipper.extractFile(new File(inputFileName));
                    uncompress = true;

                } else {
                    inputFile = new File(inputFileName);
                }
                // Added by FG compress file
                String compress = Utils.getCmdParameter(args, "compress", false);
                //set default to "false" if compress parameter is not set:
                if (compress == null) {
                    compress = "false";
                }
                if (Boolean.valueOf(compress)) {
                    if (outputFileName.endsWith(".gz")) {
                        outputFileName = outputFileName.substring(0, outputFileName.length() - 3);
                    }

                }
                // Added by FG check if path is folder

                if (args[0].equals("MzIdentMLToMzTab")) {
                    if (inputFileName != null && outputFileName != null) {
                        MzIdentMLToMzTab mzIdentMLToMzTab = new MzIdentMLToMzTab(inputFileName, outputFileName);

                    } else {
                        guiFeedback = "Error, usage: " + mzIdentMLToMzTabUsage;
                    }

                } else if (args[0].equals("FalseDiscoveryRate")) {

                    if (args.length < 7) {
                        userFeedback += fdrUsage;
                        guiFeedback = fdrUsage;
                    } else {

                        String decoyValue = Utils.getCmdParameter(args, "decoyValue", true);
                        String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);
                        String cvTerm = Utils.getCmdParameter(args, "cvTerm", true);
                        String betterScoresAreLower = Utils.getCmdParameter(args, "betterScoresAreLower", true);

                        //Removed by ARJ - this is not required here
                        //String compress = Utils.getCmdParameter(args, "compress", true);
                        //boolean compressValue = Boolean.valueOf(compress);
                        if (inputFileName != null && decoyRegex != null && decoyValue != null) {
                            FalseDiscoveryRate fdr = new FalseDiscoveryRate(inputFile.getAbsolutePath(), decoyValue, decoyRegex, cvTerm, Boolean.valueOf(betterScoresAreLower));

                            fdr.computeFDRusingJonesMethod();
                            if (outputFileName != null) {
                                fdr.writeToMzIdentMLFile(outputFileName);

                            }
                        }
                    }
                } else if (args[0].equals("FalseDiscoveryRateGlobal")) {

                    String decoyValue = Utils.getCmdParameter(args, "decoyValue", true);
                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);
                    String cvTerm = Utils.getCmdParameter(args, "cvTerm", true);
                    String betterScoresAreLower = Utils.getCmdParameter(args, "betterScoresAreLower", true);

                    String fdrLevel = Utils.getCmdParameter(args, "fdrLevel", true);
                    String proteinLevel = Utils.getCmdParameter(args, "proteinLevel", true);

                    if (inputFileName != null && decoyRegex != null && decoyValue != null && fdrLevel != null) {
                        FalseDiscoveryRateGlobal fdrGlobal = new FalseDiscoveryRateGlobal(inputFile.getAbsolutePath(), decoyValue, decoyRegex, cvTerm, Boolean.valueOf(betterScoresAreLower), fdrLevel, proteinLevel);

                        if (outputFileName != null) {
                            fdrGlobal.computeFDRusingJonesMethod();
                            fdrGlobal.writeToMzIdentMLFile(outputFileName);

                        }
                    }
                } else if (args[0].equals("GenericFasta")) {
                    String accession_regex = Utils.getCmdParameter(args, "accession_regex", true);
                    
                    // Simon: changed last parameter to false for the GenericSearch case, where we do not have a GFF file and therefore can not require it.
                    String inputGff = Utils.getCmdParameter(args, "inputGff", false);
                    
                    if (inputGff==null){
                        inputGff = "";
                    }
                    if (inputFileName != null && accession_regex != null) {
                        GenericFasta genericFasta = new GenericFasta(inputFile.getAbsolutePath(), outputFileName, accession_regex, inputGff);
                    }

                } else if (args[0].equals("CombineFastaFiles")) {

                    if (inputFileName != null) {
                        CombineFastaFiles combineFastaFiles = new CombineFastaFiles();
                        combineFastaFiles.combine(inputFileName, outputFileName);
                    }

                } else if (args[0].equals("AddGenomeCoordinatesForPeptides")) {
                    String inputGff = Utils.getCmdParameter(args, "inputGff", true);
                    String outputGff = Utils.getCmdParameter(args, "outputGff", true);
                    if (inputFileName != null && inputGff != null) {
                        AddGenomeCoordinatesForPeptides addGenomeCoordinatesForPeptides = new AddGenomeCoordinatesForPeptides(inputFile.getAbsolutePath(), outputFileName, inputGff, outputGff);
                        addGenomeCoordinatesForPeptides.writeMappingResults();
                    }

                } else if (args[0].equals("AddRetentionTimeToMzid")) {
                    String inputSourceFile = Utils.getCmdParameter(args, "inputSourceFile", true);
                    if (inputFileName != null && inputSourceFile != null) {
                        AddRetentionTimeToMzid addRetentionTimeToMzid = new AddRetentionTimeToMzid(inputFile.getAbsolutePath(), inputSourceFile, outputFileName);
                    }
                } //XtandemPercolator, OmssaPercolator, MsgfPercolator 
                else if (args[0].equals("XtandemPercolator")) {
//                    String mod_file = Utils.getCmdParameter(args, "mod_file", true);
                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);
                    String proteinCodeRegex = Utils.getCmdParameter(args, "proteinCodeRegex", false);
                    if (proteinCodeRegex != null) {
                        System.out.println("Using following regex for extracting the protein codes: " + proteinCodeRegex);
                    }
                    if (inputFileName != null && decoyRegex != null) {
                        XtandemPercolator.xtandem_percolator(inputFileName, outputFileName, decoyRegex, proteinCodeRegex);
                    }

                } else if (args[0].equals("OmssaPercolator")) {
                    // String mod_file = Utils.getCmdParameter(args, "mod_file", true);
                    String database = Utils.getCmdParameter(args, "database", true);
                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);

                    if (inputFileName != null && decoyRegex != null) {
                        OmssaPercolator.omssa_percolator(inputFileName, outputFileName, database, decoyRegex);
                    }

                } else if (args[0].equals("MsgfPercolator")) {

                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);

                    if (inputFileName != null && decoyRegex != null) {
                        MsgfPercolator.msgf_percolaotr(inputFileName, outputFileName, decoyRegex);
                    }

                } else if (args[0].equals("Mzid2Csv")) {

                    String options = Utils.getCmdParameter(args, "exportType", true);
                    String verboseOutputParam = Utils.getCmdParameter(args, "verboseOutput", false);
                    Boolean verboseOutput = false;

                    if (verboseOutputParam != null) {
                        verboseOutput = Boolean.valueOf(verboseOutputParam);
                    }

                    if (inputFileName != null && outputFileName != null && options != null) {
                        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
                        mzidToCsv.useMzIdentMLToCSV(inputFile.getAbsolutePath(), outputFileName, options, verboseOutput);
                    } else {
                        guiFeedback = "Error, usage: " + mzid2CsvUsage;
                    }

                } else if (args[0].equals("CreateRestrictedFASTADatabase")) {

                    if (inputFileName != null && outputFileName != null) {
                        CreateRestrictedFASTADatabase createRestrictedFASTADatabase = new CreateRestrictedFASTADatabase();
                        createRestrictedFASTADatabase.use(inputFile.getAbsolutePath(), outputFileName);
                    } else {
                        guiFeedback = "Error, usage: " + createRestrictedFASTADatabaseUsage;
                    }

                } else if (args[0].equals("CombinePSMMzidFiles")) {
                    String combineFractions = Utils.getCmdParameter(args, "combineFractions", true);
                    CombinePSMMzidFiles combinePSMMzidFiles = new CombinePSMMzidFiles(inputFileName, outputFileName, Boolean.valueOf(combineFractions));
                } else if (args[0].equals("ProteoGrouper")) {
                    String requireSIIsToPassThresholdParam = Utils.getCmdParameter(args, "requireSIIsToPassThreshold", true);
                    String verboseOutputParam = Utils.getCmdParameter(args, "verboseOutput", true);
                    String cvAccForSIIScoreParam = Utils.getCmdParameter(args, "cvAccForSIIScore", true);
                    String logTransScoreParam = Utils.getCmdParameter(args, "logTransScore", true);

                    String v11String = Utils.getCmdParameter(args, "version1_1", false);
                    String useProteoAnnotatorValue = Utils.getCmdParameter(args, "useProteoAnnotator", false);

                    Boolean version1_1 = true;
                    if (v11String != null) {
                        version1_1 = Boolean.parseBoolean(v11String);
                    }

                    Boolean requireSIIsToPassThreshold = Boolean.valueOf(requireSIIsToPassThresholdParam);
                    Boolean verboseOutput = Boolean.valueOf(verboseOutputParam);
                    String cvAccForSIIScore = cvAccForSIIScoreParam;
                    Boolean logTransScore = Boolean.valueOf(logTransScoreParam);
                    Boolean useProteoAnnotator = Boolean.valueOf(useProteoAnnotatorValue);

                    //ARJ 15/03/2013 Removed Boolean bestScoreForPep - not implemented properly
                    if (inputFileName != null && outputFileName != null && requireSIIsToPassThreshold != null && verboseOutput != null && cvAccForSIIScore != null
                            && logTransScore != null) {
                        if (version1_1) {
                            System.out.println("About to run protein inference and export to version 1.1.0");
                        } else {
                            System.out.println("About to run protein inference and export to version 1.2.0 draft");
                        }
                        ProteoGrouper proteinInference = new ProteoGrouper(inputFile.getAbsolutePath(), outputFileName, requireSIIsToPassThreshold, verboseOutput, cvAccForSIIScore, logTransScore, version1_1, useProteoAnnotator);

                    } else {

                        System.out.println("Error, usage: " + userFeedback + proteogrouperUsage);
                        guiFeedback = "Error, usage: " + proteogrouperUsage;
                    }

                } else if (args[0].equals("Threshold")) {
                    Boolean isPSMThreshold = Boolean.valueOf(Utils.getCmdParameter(args, "isPSMThreshold", true));
                    //String inFile, String outFile, Boolean psmThresh, String cvAccScoreThreshold, Double threshValue

                    String cvAccForScoreThreshold = Utils.getCmdParameter(args, "cvAccessionForScoreThreshold", true);
                    Double threshValue = Double.parseDouble(Utils.getCmdParameter(args, "threshValue", true));
                    Boolean scoreLowToHigh = Boolean.valueOf(Utils.getCmdParameter(args, "betterScoresAreLower", true));
                    Boolean deleteUnderThreshold = Boolean.valueOf(Utils.getCmdParameter(args, "deleteUnderThreshold", true));
                    String scoreLevel = Utils.getCmdParameter(args, "scoreLevel", false);

                    if (inputFileName != null && outputFileName != null && isPSMThreshold != null && cvAccForScoreThreshold != null && threshValue != null && scoreLowToHigh != null) {
                        ThresholdMzid thresholdMzid = new ThresholdMzid(inputFile.getAbsolutePath(), outputFileName, isPSMThreshold, cvAccForScoreThreshold, threshValue, scoreLowToHigh, deleteUnderThreshold, scoreLevel);
                    } else {
                        System.out.println("Error in parameters\n" + "Usage: " + userFeedback + thresholdUsage);
                        guiFeedback = "Error in parameters\n" + "Usage: " + thresholdUsage;
                    }
                } else if (args[0].equals("InsertMetaDataFromFasta")) {
                    String fastaFile = Utils.getCmdParameter(args, "fastaFile", true);
                    String accSplitRegex = Utils.getCmdParameter(args, "accessionSplitRegex", true);
                    accSplitRegex = accSplitRegex.replaceAll("/", "");

                    if (inputFileName != null && outputFileName != null && fastaFile != null && accSplitRegex != null) {
                        InsertMetaDataFromFasta insertMD = new InsertMetaDataFromFasta(inputFile.getAbsolutePath(), outputFileName, fastaFile, accSplitRegex);
                    } else {
                        System.out.println("Error in parameters\n" + "Usage: " + userFeedback + insertMetaDataUsage);
                        guiFeedback = "Error in parameters\n" + "Usage: " + insertMetaDataUsage;
                    }
                } else if (args[0].equals("Omssa2mzid")) {

                    Boolean outputFrags = Boolean.valueOf(Utils.getCmdParameter(args, "outputFragmentation", false));
                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", true);

                    String omssaModsFile = Utils.getCmdParameter(args, "omssaModsFile", false);
                    String userModsFile = Utils.getCmdParameter(args, "userModsFile", false);

                    if (decoyRegex != null && outputFrags != null) {
                        // Multiple conversion
                        if (inputFile.isDirectory()) {
                            File[] listOfFiles = inputFile.listFiles();
                            for (int i = 0; i < listOfFiles.length; i++) {
                                File file = listOfFiles[i];

                                if (file.getAbsolutePath().endsWith("omx")) {
                                    String outputName = file.getName().substring(0, file.getName().lastIndexOf('.'));
                                    String outputFileName1 = outputFileName + "\\" + outputName + ".mzid";
                                    Omssa2mzid omssa2mzid = new Omssa2mzid(file.getAbsolutePath(), outputFileName1, outputFrags, decoyRegex, omssaModsFile, userModsFile);
                                    System.out.println("Searching " + file.getName() + " for decoys containing " + decoyRegex);
                                }

                            }
                        } else {
                            new Omssa2mzid(inputFile.getAbsolutePath(), outputFileName, outputFrags, decoyRegex, omssaModsFile, userModsFile);
                            System.out.println("Searching for decoys containing " + decoyRegex);
                        }
                    } /*
                     * ARJ - only one constructor to maintain now else if
                     * (outputFrags != null) { new
                     * Omssa2mzid(inputFile.getAbsolutePath(), outputFileName,
                     * outputFrags); }
                     *
                     */ else {
                        System.out.println("Error in parameters\n" + "Usage: " + userFeedback + omssa2mzidUsage);
                        guiFeedback = "Error in parameters\n" + "Usage: " + omssa2mzidUsage;
                    }

                } else if (args[0].equals("Tandem2mzid")) {

                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", false);
                    if (decoyRegex != null) {
                        System.out.println("Searching for decoys containing " + decoyRegex);
                    }
                    String proteinCodeRegex = Utils.getCmdParameter(args, "proteinCodeRegex", false);
                    if (proteinCodeRegex != null) {
                        System.out.println("Using following regex for extracting the protein codes: " + proteinCodeRegex);
                    }
                    String outputFragsValue = Utils.getCmdParameter(args, "outputFragmentation", false);
                    boolean outputFrags = true; //true is the default 

                    if (outputFragsValue != null) {
                        outputFrags = (outputFragsValue.equalsIgnoreCase("true") ? true : false);
                    }

                    //the parameter idsStartAtZero can be used to overrule the default behavior of 
                    //the tool, which is to assume that spectrumIds start at zero if the spectra 
                    //file format is mzML, start at one otherwise. 
                    //Now in some cases it could be that your specific spectra format also starts at
                    //zero. Then set this to TRUE. Anyway, for mzML spectra file this parameter 
                    //should be TRUE, for other formats (like mzData) the user should check whether the 
                    //spectrum ids start at zero or at 1.
                    Boolean idsStartAtZero = null; //null is the default    
                    String idsStartAtZeroValue = Utils.getCmdParameter(args, "idsStartAtZero", false);
                    if (idsStartAtZeroValue != null) {
                        idsStartAtZero = Boolean.valueOf(idsStartAtZeroValue.toLowerCase());
                    }

                    String databaseFileFormatID = Utils.getCmdParameter(args, "databaseFileFormatID", false);
                    String massSpecFileFormatID = Utils.getCmdParameter(args, "massSpecFileFormatID", false);
                    try {
                        // Multiple conversion
                        if (inputFile.isDirectory()) {
                            File[] listOfFiles = inputFile.listFiles();
                            for (int i = 0; i < listOfFiles.length; i++) {
                                File file = listOfFiles[i];

                                if (file.getAbsolutePath().endsWith(".t.xml")) {
                                    String outputName = file.getName().substring(0, file.getName().lastIndexOf('.'));
                                    String outputFileName1 = outputFileName + "\\" + outputName + ".mzid";
                                    Tandem2mzid tandem2mzid = new Tandem2mzid(file.getAbsolutePath(), outputFileName1, databaseFileFormatID, massSpecFileFormatID, idsStartAtZero, decoyRegex, proteinCodeRegex, outputFrags);
                                    System.out.println("Searching " + file.getName() + " for decoys containing " + decoyRegex);

                                }

                            }
                        } else {
                            Tandem2mzid tandem2mzid = new Tandem2mzid(inputFile.getAbsolutePath(), outputFileName, databaseFileFormatID, massSpecFileFormatID, idsStartAtZero, decoyRegex, proteinCodeRegex, outputFrags);
                        }
                    } catch (Exception e) {
                        System.out.println("Error running Tandem2mzid:" + userFeedback + tandem2mzidUsage);
                        e.printStackTrace();
                        guiFeedback = "Error running Tandem2mzid:" + tandem2mzidUsage;
                    }

                } else if (args[0].equals("Csv2mzid")) {

                    String paramsFile = Utils.getCmdParameter(args, "paramsFile", true);
                    String decoyRegex = Utils.getCmdParameter(args, "decoyRegex", false);
                    String cvAccForPSMOrdering = Utils.getCmdParameter(args, "cvAccessionForPSMOrdering", true);

                    Boolean applyFixedMods = true;
                    String applyFixedModsString = Utils.getCmdParameter(args, "applyFixedMods", false);
                    if (applyFixedModsString != null) {
                        applyFixedMods = Boolean.parseBoolean(applyFixedModsString);
                    }

                    if (inputFileName != null && outputFileName != null && paramsFile != null && cvAccForPSMOrdering != null && decoyRegex != null) {
                        new Csv2mzid(inputFile.getAbsolutePath(), outputFileName, paramsFile, cvAccForPSMOrdering, decoyRegex, applyFixedMods);
                    } else if (inputFileName != null && outputFileName != null && paramsFile != null && cvAccForPSMOrdering != null) {
                        new Csv2mzid(inputFile.getAbsolutePath(), outputFileName, paramsFile, cvAccForPSMOrdering, applyFixedMods);
                    } else {
                        System.out.println("Error, usage: " + userFeedback + csv2mzidUsage);
                        guiFeedback = "Error, usage: " + csv2mzidUsage;
                    }
                } else if (args[0].equals("CombineSearchEngines") && args.length == 25) {
                    //String[] newArgs = new String[11];
                    String[] newArgs = new String[13];
                    newArgs[0] = Utils.getCmdParameter(args, "firstFile", true);
                    newArgs[1] = "s1";
                    newArgs[2] = Utils.getCmdParameter(args, "firstcvTerm", true);
                    newArgs[3] = Utils.getCmdParameter(args, "firstbetterScoresAreLower", true);
                    newArgs[4] = Utils.getCmdParameter(args, "secondFile", true);
                    newArgs[5] = "s2";
                    newArgs[6] = Utils.getCmdParameter(args, "secondcvTerm", true);
                    newArgs[7] = Utils.getCmdParameter(args, "secondbetterScoresAreLower", true);
                    newArgs[8] = Utils.getCmdParameter(args, "rank", true);
                    newArgs[9] = Utils.getCmdParameter(args, "decoyRatio", true);
                    newArgs[10] = Utils.getCmdParameter(args, "outputFile", true);
                    newArgs[11] = Utils.getCmdParameter(args, "debugFile", true);
                    newArgs[12] = Utils.getCmdParameter(args, "decoyRegex", true);

                    CombineSearchEngines.runTwoSearchEngines(newArgs);

                } else if (args[0].equals("CombineSearchEngines") && args.length == 31) {
                    //String[] newArgs = new String[14];
                    String[] newArgs = new String[17];
                    newArgs[0] = Utils.getCmdParameter(args, "firstFile", true);
                    newArgs[1] = "s1";
                    newArgs[2] = Utils.getCmdParameter(args, "firstcvTerm", true);
                    newArgs[3] = Utils.getCmdParameter(args, "firstbetterScoresAreLower", true);
                    newArgs[4] = Utils.getCmdParameter(args, "secondFile", true);
                    newArgs[5] = "s2";
                    newArgs[6] = Utils.getCmdParameter(args, "secondcvTerm", true);
                    newArgs[7] = Utils.getCmdParameter(args, "secondbetterScoresAreLower", true);
                    newArgs[8] = Utils.getCmdParameter(args, "thirdFile", true);
                    newArgs[9] = "s3";
                    newArgs[10] = Utils.getCmdParameter(args, "thirdcvTerm", true);
                    newArgs[11] = Utils.getCmdParameter(args, "thirdbetterScoresAreLower", true);
                    newArgs[12] = Utils.getCmdParameter(args, "rank", true);
                    newArgs[13] = Utils.getCmdParameter(args, "decoyRatio", true);
                    newArgs[14] = Utils.getCmdParameter(args, "outputFile", true);
                    newArgs[15] = Utils.getCmdParameter(args, "debugFile", true);
                    newArgs[16] = Utils.getCmdParameter(args, "decoyRegex", true);
                    CombineSearchEngines.runThreeSearchEngines(newArgs);

                } else if (args[0].equals("AddEmpaiToMzid")) {

                    //AddEmpaiToMzid(String mzidIn, String mzidOut, String fastaIn, String accRegex, String enzRegex){
                    String fastaFile = Utils.getCmdParameter(args, "fastaFile", true);
                    String accSplitRegex = Utils.getCmdParameter(args, "accessionSplitRegex", true);
                    accSplitRegex = accSplitRegex.replaceAll("/", "");
                    String enzRegex = Utils.getCmdParameter(args, "enzymeRegex", false);

                    String verboseMode = Utils.getCmdParameter(args, "verboseOutput", false);

                    Boolean verbose = false;
                    if (verboseMode != null) {
                        verbose = Boolean.parseBoolean(verboseMode);
                    }

                    if (inputFileName != null && outputFileName != null && fastaFile != null && accSplitRegex != null) {
                        AddEmpaiToMzid empai = new AddEmpaiToMzid(inputFile.getAbsolutePath(), outputFileName, fastaFile, accSplitRegex, enzRegex, verbose);
                    } else {
                        System.out.println("Error in parameters\n" + "Usage: " + userFeedback + emPAIUsage);
                        guiFeedback = "Error, usage: " + emPAIUsage;
                    }
                } else if (args[0].equals("ProteoAnnotator")) {
                    //Calling ProteoAnnotator
                    String inputGFF = Utils.getCmdParameter(args, "inputGFF", false);
                    String inputFasta = Utils.getCmdParameter(args, "inputFasta", false);
                    String spectrum_files = Utils.getCmdParameter(args, "spectrum_files", true);
                    String outputFolder = Utils.getCmdParameter(args, "outputFolder", true);
                    String inputPredicted = Utils.getCmdParameter(args, "inputPredicted", false);
                    String searchParameters = Utils.getCmdParameter(args, "searchParameters", true);
                    String prefix = Utils.getCmdParameter(args, "prefix", false);
                    String peptideThreshValue = Utils.getCmdParameter(args, "peptideThreshValue", true);
                    String proteinThreshValue = Utils.getCmdParameter(args, "proteinThreshValue", true);

                    ProteoAnnotator proteoAnnotator = new ProteoAnnotator(inputGFF, inputFasta, spectrum_files, outputFolder, inputPredicted, searchParameters, prefix, peptideThreshValue, proteinThreshValue);
                    proteoAnnotator.runProteoAnnotator();

                } else if (args[0].equals("GenericSearch")) {
                    //Calling GenericSearch
                    String inputFasta = Utils.getCmdParameter(args, "inputFasta", false);
                    String spectrum_files = Utils.getCmdParameter(args, "spectrum_files", true);
                    String outputFolder = Utils.getCmdParameter(args, "outputFolder", true);
                    String searchParameters = Utils.getCmdParameter(args, "searchParameters", true);
                    String prefix = Utils.getCmdParameter(args, "prefix", false);
                    String peptideThreshValue = Utils.getCmdParameter(args, "peptideThreshValue", true);
                    String proteinThreshValue = Utils.getCmdParameter(args, "proteinThreshValue", true);

                    GenericSearch genericSearch = new GenericSearch(inputFasta, spectrum_files, outputFolder, searchParameters, prefix, peptideThreshValue, proteinThreshValue);
                    genericSearch.runGenericSearch();

                } else if (args[0].equals("RescoreMods")) {
                    //-cvAccForScoreToAdapt MS:1002053  -logTransformPSMScore true|false -commonModificationWeight 2.0 -mediumModificationWeight 0.5 -rareModificationWeight 0.1 -generalModificationWeight 0.5 -pairedModificationAndUnmodWeight 3.0 -multipleVariableModWeight 0.1 
                    //AddEmpaiToMzid(String mzidIn, String mzidOut, String fastaIn, String accRegex, String enzRegex){
                    String cvForScore = Utils.getCmdParameter(args, "cvAccForScoreToAdapt", true);
                    //Boolean logTransform = Boolean.parseBoolean(Utils.getCmdParameter(args, "logTransformPSMScore", true));
                    Double commonModWeight = Double.parseDouble(Utils.getCmdParameter(args, "commonModificationWeight", true));
                    Double mediumModWeight = Double.parseDouble(Utils.getCmdParameter(args, "mediumModificationWeight", true));
                    Double rareModWeight = Double.parseDouble(Utils.getCmdParameter(args, "rareModificationWeight", true));
                    Double generalModWeight = Double.parseDouble(Utils.getCmdParameter(args, "generalModificationWeight", true));
                    Double pairedModWeight = Double.parseDouble(Utils.getCmdParameter(args, "pairedModificationAndUnmodWeight", true));
                    Double multipleModWeight = Double.parseDouble(Utils.getCmdParameter(args, "multipleVariableModWeight", true));

                    if (inputFileName != null && outputFileName != null) {
                        System.out.println("Input:" + inputFile.getAbsolutePath());
                        System.out.println("Output:" + outputFileName);
                        RescoreMods rescore = new RescoreMods(inputFile.getAbsolutePath(), outputFileName, cvForScore, commonModWeight, mediumModWeight, rareModWeight, generalModWeight, pairedModWeight, multipleModWeight);
                    } else {
                        System.out.println("Error in parameters\n" + "Usage: " + userFeedback + emPAIUsage);
                        guiFeedback = "Error, usage: " + rescoreModsUsage;
                    }
                } else {
                    String tempFeedback = "Program within MzidLib not recognized: " + args[0] + omssa2mzidUsage + "\n" + fdrUsage + "\n" + mzid2CsvUsage + "\n" + proteogrouperUsage + "\n" + thresholdUsage + "\n" + insertMetaDataUsage
                            + "\n" + tandem2mzidUsage + "\n" + csv2mzidUsage + "\n" + emPAIUsage;

                    System.out.println(tempFeedback);
                    guiFeedback = "Error, usage: " + tempFeedback;
                    userFeedback = tempFeedback;
                }
                if (userFeedback.equals("")) {
                    userFeedback = "Completed successfully, output written to " + outputFileName;
                }
                System.out.println(userFeedback);
                // Added by FG delete tmp file
                if (uncompress) {
//                    Gzipper.deleteFile(inputFile);
                    inputFile.deleteOnExit();
                }

                if (Boolean.valueOf(compress)) {
                    Gzipper.compressFile(new File(outputFileName));
                }
            }

        } else {
            String tempFeedback = "Error insufficient arguments entered, options: " + userFeedback + " toolname options\n"
                    + "\n\nTools:\n\n***************\n"
                    + omssa2mzidUsage
                    + "\n" + "***************\n" + fdrUsage
                    + "\n" + "***************\n" + mzid2CsvUsage
                    + "\n" + "***************\n" + proteogrouperUsage
                    + "\n" + "***************\n" + thresholdUsage
                    + "\n" + "***************\n" + insertMetaDataUsage
                    + "\n" + "***************\n" + tandem2mzidUsage
                    + "\n" + "***************\n" + csv2mzidUsage
                    + "\n" + "***************\n" + emPAIUsage
                    + "\n" + "***************\n" + omssaPercolatorUsage
                    + "\n" + "***************\n" + xtandemPercolatorUsage
                    + "\n" + "***************\n" + combinedSearchUsage
                    + "\n" + "***************\n" + createRestrictedFASTADatabaseUsage
                    + "\n" + "***************\n" + mzIdentMLToMzTabUsage
                    + "\n" + "***************\n" + rescoreModsUsage
                    + "\n" + "***************\n" + fdrGlobalUsage
                    + "\n" + "***************\n" + combinePSMMzidFilesUsage
                    + "\n" + "***************\n" + genericFastaUsage
                    + "\n" + "***************\n" + addGenomeCoordinatesForPeptidesUsage
                    + "\n" + "***************\n" + addRetentionTimeToMzidUsage
                    + "\n" + "***************\n" + combineFastaFilesUsage
                    + "\n" + "***************\n" + msgfPercolatorUsage + "\n";

            System.out.println(tempFeedback);
            guiFeedback = tempFeedback;
        }

        return guiFeedback;

    }

}
