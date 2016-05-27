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
import java.util.Map;
import java.util.Set;

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

    public static Map<String, String> getInitializedCVMap() throws IOException {
        //Read resource file and build up map:
        BufferedReader in = null;
        Map<String, String> resultMap = new HashMap<String, String>();
        try {
            //Use the getResourceAsStream trick to read the file packaged in
            //the .jar .  This simplifies usage of the solution as no extra 
            //classpath or path configurations are needed: 
            InputStream resourceAsStream = ClassLoader.getSystemClassLoader().getResourceAsStream("CV_psi-ms.obo.txt");
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
                        throw new RuntimeException("Unexpected name: preceding id: entry in CV file");
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

    public static String getCmdParameter(String[] args, String name, boolean required) {
        for (int i = 0; i < args.length; i++) {
            String argName = args[i];
            if (argName.equals("-" + name)) {
                String argValue = "";
                if (i + 1 < args.length) {
                    argValue = args[i + 1];
                }
                if (required && (argValue.trim().length() == 0 || argValue.startsWith("-"))) {
                    System.err.println("Parameter value expected for " + argName);
                    throw new RuntimeException("Expected parameter value not found: " + argName);
                } else if (argValue.trim().length() == 0 || argValue.startsWith("-")) {
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

    public static File splitMGFsOrReturnSame(String path, File masterMgf, int fileUpperLimit, int spectraUpperLimit) throws IOException {
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
            Path newMgfPath = temporaryMgfFolderPath.resolve("peaks_" + i + ".mgf");
            writers.add(Files.newBufferedWriter(newMgfPath, StandardCharsets.UTF_8));
        }

        BufferedReader reader = Files.newBufferedReader(Paths.get(masterMgf.getAbsolutePath()), StandardCharsets.UTF_8);
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

    private static int countMgfEntries(File mgf) throws IOException {
        BufferedReader reader = Files.newBufferedReader(Paths.get(mgf.getAbsolutePath()), StandardCharsets.UTF_8);
        String line = null;
        int ionCount = 0;
        while ((line = reader.readLine()) != null) {
            if (line.trim().equalsIgnoreCase("BEGIN IONS")) {
                ionCount++;
            }
        }

        reader.close();
        return ionCount;
    }

}
