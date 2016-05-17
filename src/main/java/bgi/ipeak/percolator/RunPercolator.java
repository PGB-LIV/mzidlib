/**
 *
 */
package bgi.ipeak.percolator;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.io.StreamHog;
import bgi.ipeak.util.Properties;

/**
 *
 * @author Administrator
 *
 */
public class RunPercolator {

    @Option(name = "-f", required = true, usage = "(required) Path and file name of tab delimited features file")
    private File feature = null;
    @Option(name = "-i", required = false, usage = "(required) Results path and file name (without extension)")
    private String outFilesName = null;
    @Option(name = "-x", required = false, usage = "(optional) Write supplemental XML output as defined here: http://noble.gs.washington.edu/proj/percolator/model/percolator_out.xsd")
    private String temp_xml_file = null;
    @Option(name = "-u", required = true, usage = "Activates -U switch when running percolator (older percolator version it switches peptide grouping on newer versions it switches it off) default off")
    private boolean us = false;

    /**
     * Constructors for RunPercolator.
     *
     * @param features Features file.
     * @param out Output file.
     * @param tempXML Temp XML file.
     * @param useU Use 'U' switch.
     * @throws IOException IOException thrown my constructor.
     * @throws InterruptedException InterruptedException thrown by constructor.
     */
    public RunPercolator(String features, String out, String tempXML, boolean useU) throws IOException, InterruptedException {
        this.feature = new File(features);
        this.temp_xml_file = tempXML;
        this.outFilesName = out;
        this.us = useU;
    }

    public RunPercolator() {

    }

    public static void main(String[] args) throws IOException, InterruptedException {
        RunPercolator rp = new RunPercolator();
        CmdLineParser parser = new CmdLineParser(rp);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.RunPercolator [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        rp.execute();
    }

    public int execute() throws IOException, InterruptedException {
        FileWriter scoreFile = new FileWriter(new File(this.outFilesName + ".target.txt"));
        FileWriter logFile = new FileWriter(new File(this.outFilesName + ".log.txt"));

        int ret = execute(Properties.getPercolator_path(), this.us, this.temp_xml_file, this.outFilesName + ".decoy.txt", this.feature.getAbsolutePath(), scoreFile, logFile);
        if (ret != 0) {
            System.err.println("Error exists while running percolator: " + ret);
            System.exit(0);
        }
        scoreFile.flush();
        scoreFile.close();
        logFile.flush();
        logFile.close();
        return ret;
    }

    public static int execute(String percolatorPath, boolean uSwitch, String tempxml, String decoyResults, String features,
            FileWriter stdoutStream, FileWriter percolatorLogStream)
            throws IOException, InterruptedException {
        String uFlag = "";
        if (uSwitch) {
            uFlag = " -U ";
        }
        String xmlFlag = "";
        if ((tempxml != null) && (tempxml.trim().length() > 0)) {
            xmlFlag = " -X " + tempxml + " ";
        }

        String decoyResultsFlag = "";
        if ((decoyResults != null) && (decoyResults.trim().length() > 0)) {
            decoyResultsFlag = " -B " + decoyResults + " ";
        }

        String cmd = percolatorPath + " -j " + features + xmlFlag + decoyResultsFlag + uFlag;
        System.out.println("Run percolator: " + cmd);
        Runtime rt = Runtime.getRuntime();
        Process p = null;
        try {
            p = rt.exec(cmd);
            StreamHog errorHog = new StreamHog(p.getErrorStream(), percolatorLogStream);
            StreamHog outputHog = new StreamHog(p.getInputStream(), stdoutStream);
            p.waitFor();
            errorHog.join();
            outputHog.join();
        } catch (IOException e) {
            System.err.println("Can not run this command: '" + cmd + "'");
            throw e;
        }
        int exitValue;
        if ((exitValue = p.exitValue()) != 0) {
            System.err.println("Error with Run Percolator: " + p.exitValue());
            System.exit(p.exitValue());
        }

        return exitValue;
    }
}
