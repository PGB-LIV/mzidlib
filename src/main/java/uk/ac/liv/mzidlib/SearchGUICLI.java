
package uk.ac.liv.mzidlib;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/**
 * A class to call SearchGUI from the command line
 *
 * @author Fawaz Ghali 09-Apr-2014
 */
public class SearchGUICLI {

    private String searchGUIPath = "";
    private String outputFolder = "";
    private PrintWriter out = null;
    private String debugFile = "";

    public SearchGUICLI(String outputFolder, String debugFile) {
        try {
            this.outputFolder = outputFolder;

            searchGUIPath = getSearchGUIFile();
            // make sure that the SearchGUI jar file is executable
            File searchGUIExecutable = new File(searchGUIPath);
            searchGUIExecutable.setExecutable(true);
            out = new PrintWriter(new BufferedWriter(new FileWriter(debugFile,
                                                                    true)));
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Returns the path to the jar file's parent folder.
     *
     * @return the path to the jar file's parent folder
     */
    public String getSearchGUIPath() {
        return new File(SearchGUICLI.class.getProtectionDomain().getCodeSource()
                .getLocation().getPath()).getParent() + File.separator
                + "searchgui";
    }

    /**
     * Returns SearchGUI jar file.
     *
     * @return SearchGUI jar file.
     */
    public String getSearchGUIFile() {
        String path = null;
        try {
            File dir = new File(URLDecoder.decode(getSearchGUIPath(), "UTF-8"));
            System.out.println("SearchGUI dir: " + dir);

            File[] matches = dir.listFiles(new FilenameFilter() {

                public boolean accept(File dir, String name) {
                    return name.endsWith(".jar");
                }

            });
            path = matches[0].getAbsolutePath();

        } catch (UnsupportedEncodingException ex) {
            ex.printStackTrace();
        }
        return path;
    }

    /**
     * Returns output file name from the decoy function
     *
     * @param in input fasta file
     *
     * @return output file name from the decoy function
     */
    public String runDeocyCLI(String in) {
        String outputFileName = "";
        try {
            System.out.println(
                    "Running SearchGUI decoy function using the following command:");
            //out.println("Running SearchGUI decoy function using the following command:");
            String args[] = {"java", "-Djava.awt.headless=true", "-cp",
                searchGUIPath, "eu.isas.searchgui.cmd.FastaCLI", "-in", in,
                "-decoy"};
            for (int i = 0; i < args.length; i++) {
                System.out.print(args[i] + " ");
                out.print(args[i] + " ");
            }
            System.out.println();
            out.println();
            ProcessBuilder pb = new ProcessBuilder(args);
            pb.redirectErrorStream(true);
            final Process process = pb.start();
            InputStream is = process.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) {
                System.out.println(line);
                out.println(line);
            }
            try {
                process.waitFor();

            } catch (InterruptedException ex) {
                if (process != null) {
                    process.destroy();
                }

                ex.printStackTrace();
            }
            System.out.println("SearchGUI decoy function finished.");

            //out.println("SearchGUI decoy function finished.");
        } catch (IOException ex) {

            ex.printStackTrace();
        }
        String fileNameWithOutExt = FilenameUtils.removeExtension(in);
        outputFileName = fileNameWithOutExt + "_concatenated_target_decoy.fasta";
        return outputFileName;
    }

    void runParameterFileCLI(String decoyFasta, String outputParameterFile,
                             String paramters) {
        try {
            System.out.println(
                    "Running SearchGUI ParameterFileCLI function using the following command:");
            //out.println("Running SearchGUI ParameterFileCLI function using the following command:");
            String javaCLI[] = {"java", "-Djava.awt.headless=true", "-cp",
                searchGUIPath,
                "eu.isas.searchgui.cmd.IdentificationParametersCLI", "-db",
                decoyFasta, "-out", outputParameterFile};

            BufferedReader br1 = new BufferedReader(new FileReader(paramters));
            String everything = "";
            try {
                StringBuilder sb = new StringBuilder();
                String line = br1.readLine();

                while (line != null) {
                    sb.append(line);
                    sb.append(System.lineSeparator());
                    line = br1.readLine();
                }
                everything = sb.toString();
            } finally {
                br1.close();
            }
            System.out.println("Search parameters are: " + everything);
            List<String> list = new ArrayList<>();
            Matcher m = Pattern.compile("([^\"]\\S*|\".+?\")\\s*").matcher(
                    everything);
            while (m.find()) {
                list.add(m.group(1)); // Add .replace("\"", "") to remove surrounding quotes.
            }

            System.out.println(list);
            String[] args = (String[]) ArrayUtils
                    .addAll(javaCLI, list.toArray());
            for (int i = 0; i < args.length; i++) {
                System.out.print(args[i] + " ");
                //out.print(args[i] + " ");
            }
            System.out.println();
            //out.println();
            ProcessBuilder pb = new ProcessBuilder(args);
            pb.redirectErrorStream(true);
            final Process process = pb.start();
            InputStream is = process.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }
            try {
                process.waitFor();

            } catch (InterruptedException ex) {
                if (process != null) {
                    process.destroy();
                }

                ex.printStackTrace();
            }
            System.out.println("SearchGUI ParameterFileCLI function finished.");
            //out.println("SearchGUI ParameterFileCLI function finished.");
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    public void runSearchGUICLI(String[] searchParamters) {
        try {
            System.out.println(
                    "Running SearchGUI search function using the following command:");
            out.println(
                    "Running SearchGUI search function using the following command:");
            String javaCLI[] = {"java", "-Djava.awt.headless=true", "-cp",
                searchGUIPath, "eu.isas.searchgui.cmd.SearchCLI"};
            String[] args = (String[]) ArrayUtils.addAll(javaCLI,
                                                         searchParamters);

            for (int i = 0; i < args.length; i++) {
                System.out.print(args[i] + " ");
                //out.print(args[i] + " ");
            }
            System.out.println();
            //out.println();

            ProcessBuilder pb = new ProcessBuilder(args);
            pb.redirectErrorStream(true);
            final Process process = pb.start();
            InputStream is = process.getInputStream();
            InputStreamReader isr = new InputStreamReader(is,
                                                          StandardCharsets.UTF_8);
            try (BufferedReader br = new BufferedReader(isr)) {
                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                    out.println(line);
                }
                try {
                    process.waitFor();

                } catch (InterruptedException ex) {
                    process.destroy();
                    is.close();
                }

            }
            System.out.println("SearchGUI search function finished.");
            out.println("SearchGUI search function finished.");
            out.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

}
