
package uk.ac.liv.mzidlib.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.ParamList;
import uk.ac.liv.mzidlib.constants.CvConstants;

/**
 * General utilities class
 *
 * @author lukas007
 *
 */
public class Utils {

    public static double round(double value, int numberOfDecimalPlaces) {
        double multipicationFactor = Math.pow(10, numberOfDecimalPlaces);
        return Math.round(value * multipicationFactor) / multipicationFactor;
    }

    public static Map<String, String> getInitializedCVMap()
            throws IOException {
        //Read resource file and build up map:
        BufferedReader in = null;
        Map<String, String> resultMap = new HashMap<>();
        try {
            //Use the getResourceAsStream trick to read the file packaged in
            //the .jar .  This simplifies usage of the solution as no extra 
            //classpath or path configurations are needed: 
            InputStream resourceAsStream = ClassLoader.getSystemClassLoader().
                    getResourceAsStream("CV_psi-ms.obo.txt");
            Reader reader = new InputStreamReader(resourceAsStream);
            in = new BufferedReader(reader);
            String inputLine;
            String key = "";
            String value = "";

            while ((inputLine = in.readLine()) != null) {
                if (inputLine.startsWith("id:")) {
                    key = inputLine.split("id:")[1].trim();
                }
                if (inputLine.startsWith("name:")) {
                    //validate:
                    if (key.equals("")) {
                        throw new RuntimeException(
                                "Unexpected name: preceding id: entry in CV file");
                    }
                    value = inputLine.split("name:")[1].trim();
                    resultMap.put(key, value);
                    //reset:
                    key = "";
                    value = "";
                }
            }
            return resultMap;

        } finally {
            if (in != null) {
                in.close();
            }
        }

    }

    public static String getCmdParameter(String[] args, String name,
                                         boolean required) {
        for (int i = 0; i < args.length; i++) {
            String argName = args[i];
            if (argName.equals("-" + name)) {
                String argValue = "";
                if (i + 1 < args.length) {
                    argValue = args[i + 1];
                }
                if (required && (argValue.trim().length() == 0 || argValue.
                        startsWith("-"))) {
                    System.err.
                            println("Parameter value expected for " + argName);
                    throw new RuntimeException(
                            "Expected parameter value not found: " + argName);
                } else if (argValue.trim().length() == 0 || argValue.startsWith(
                        "-")) {
                    return "";
                } else {
                    return argValue;
                }
            }
        }
        //Nothing found, if required, throw error, else return "";
        if (required) {
            System.err.println("Parameter -" + name + " expected ");
            throw new RuntimeException("Expected parameter not found: " + name);
        }

        return null;
    }

    /**
     * Split the input mgf file according to fileUpperLimit and
     * spectraUpperLimit.
     * If the input file is below the size limit, then the original file is
     * output.
     * Otherwise, the file is split within the limit into several files and
     * saved in the designated path.
     *
     * @param path              the folder path for the new files
     * @param masterMgf         the original mgf file
     * @param fileUpperLimit    file size upper limit
     * @param spectraUpperLimit spectra size upper limit
     *
     * @return the same file if it is below the fileUpperLimit and
     *         spectraUpperLimit, otherwise the folder holding the split files.
     *
     * @throws IOException
     */
    public static File splitMGFsOrReturnSame(String path, File masterMgf,
                                             int fileUpperLimit,
                                             int spectraUpperLimit)
            throws IOException {
        System.out.println("splitMGFsOrReturnSame: ");

        System.out.println("path: " + path);
        System.out.println("masterMgf: " + masterMgf);
        int entries = countMgfEntries(masterMgf);
        if (entries < spectraUpperLimit && masterMgf.length() < fileUpperLimit) {
            return masterMgf;
        }

        double entriesPerFile = entries;
        int filesToSplitInto = 1;
        while (entriesPerFile > spectraUpperLimit) {
            filesToSplitInto++;
            entriesPerFile = ((double) entries) / ((double) filesToSplitInto);
            entriesPerFile = Math.ceil(entriesPerFile);
        }

        Path temporaryMgfFolderPath = Paths.get(path);
        Set<BufferedWriter> writers = new HashSet<>();
        for (int i = 0; i < filesToSplitInto; i++) {
            Path newMgfPath = temporaryMgfFolderPath.resolve("peaks_" + i
                    + ".mgf");
            writers.add(Files.newBufferedWriter(newMgfPath,
                                                StandardCharsets.UTF_8));
        }

        BufferedReader reader = Files.newBufferedReader(Paths.get(masterMgf.
                getAbsolutePath()), StandardCharsets.UTF_8);
        String line = null;

        Iterator<BufferedWriter> writerIterator = writers.iterator();
        BufferedWriter writer = null;
        int entriesInThisFile = 0;
        while ((line = reader.readLine()) != null) {
            if (writer == null) {
                if (writerIterator.hasNext()) {
                    writer = writerIterator.next();
                } else {
                    break;
                }
            }

            writer.write(line);
            writer.newLine();

            if (line.trim().equalsIgnoreCase("END IONS")) {
                entriesInThisFile++;
                if (entriesInThisFile == entriesPerFile) {
                    writer = null;
                    entriesInThisFile = 0;
                }
            }
        }

        for (BufferedWriter writerToClose : writers) {
            writerToClose.close();
        }

        return temporaryMgfFolderPath.toFile();
    }

    public static Enzyme getXtandemEnzyme(String tandemRegex, int missedCleavage) {

        Enzyme enzyme = new Enzyme();
        //[KR]|{P}

        //TODO only trypsin implemented - this is difficult to convert from Regex used in X!Tandem
        enzyme.setId("Enz1");
        enzyme.setCTermGain("OH");
        enzyme.setNTermGain("H");
        enzyme.setMissedCleavages(missedCleavage);
        enzyme.setSemiSpecific(false);

        if (tandemRegex.equalsIgnoreCase("[KR]|{P}")) {

            ParamList paramList = enzyme.getEnzymeName();
            if (paramList == null) {
                paramList = new ParamList();
                enzyme.setEnzymeName(paramList);
            }
            List<CvParam> cvParamList = paramList.getCvParam();
            cvParamList.add(MzidLibUtils.makeCvParam("MS:1001251", "Trypsin",
                                                     CvConstants.PSI_CV));
        } else {
            //TODO
            /*
             *
             * [Term] id: MS:1001303 name: Arg-C is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001272 ! (?<=R)(?!P)
             *
             * [Term] id: MS:1001304 name: Asp-N is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001273 ! (?=[BD])
             *
             * [Term] id: MS:1001305 name: Asp-N_ambic is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001274 !
             * (?=[DE])
             *
             * [Term] id: MS:1001306 name: Chymotrypsin is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001332 !
             * (?<=[FYWL])(?!P)
             *
             * [Term] id: MS:1001307 name: CNBr is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001333 ! (?<=M)
             *
             * [Term] id: MS:1001308 name: Formic_acid is_a: MS:1001045 !
             * cleavage agent name relationship: has_regexp MS:1001334 !
             * ((?<=D))|((?=D))
             *
             * [Term] id: MS:1001309 name: Lys-C is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001335 ! (?<=K)(?!P)
             *
             * [Term] id: MS:1001310 name: Lys-C/P is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001336 ! (?<=K)
             *
             * [Term] id: MS:1001311 name: PepsinA is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001337 ! (?<=[FL])
             *
             * [Term] id: MS:1001312 name: TrypChymo is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001338 !
             * (?<=[FYWLKR])(?!P)
             *
             * [Term] id: MS:1001313 name: Trypsin/P is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001339 ! (?<=[KR])
             *
             * [Term] id: MS:1001314 name: V8-DE is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001340 ! (?<=[BDEZ])(?!P)
             *
             * [Term] id: MS:1001315 name: V8-E is_a: MS:1001045 ! cleavage
             * agent name relationship: has_regexp MS:1001341 ! (?<=[EZ])(?!P)
             */
        }

        return enzyme;

    }

    private static int countMgfEntries(File mgf)
            throws IOException {
        int ionCount;
        try (BufferedReader reader
                = Files.newBufferedReader(Paths.get(mgf.getAbsolutePath()),
                                          StandardCharsets.UTF_8)) {
            String line = null;
            ionCount = 0;
            while ((line = reader.readLine()) != null) {
                if (line.trim().equalsIgnoreCase("BEGIN IONS")) {
                    ionCount++;
                }
            }
        }
        return ionCount;
    }

}
