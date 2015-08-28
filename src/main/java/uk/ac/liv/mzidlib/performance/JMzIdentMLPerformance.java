package uk.ac.liv.mzidlib.performance;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import static uk.ac.ebi.jmzidml.MzIdentMLElement.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.Callable.MzidCallable;
import uk.ac.liv.mzidlib.MzIdentMLLib;

/**
 *
 * @author Fawaz Ghali 03-Feb-2015
 */
public class JMzIdentMLPerformance {

    private static MzIdentMLUnmarshaller unmarshaller;
    private long startTime;
    private long stopTime;

    private long difference;

    private String inputDir = "C:\\Users\\fghali\\Desktop\\performance";
    private String outDir = "C:\\Users\\fghali\\Desktop\\performance";

    List<MzidCallable> collection = new ArrayList();
    int cores = 1;

    private ExecutorService executorPool;
    private StringBuffer buffer = new StringBuffer();

    // Constructor 
    public JMzIdentMLPerformance(String in, String out) {
        inputDir = in;
        outDir = out;

    }

    public String runXtandem() throws InterruptedException {

        List<String> tandemFiles = new ArrayList();

        MzIdentMLLib mzidLib = new MzIdentMLLib();
        File dir1 = new File(inputDir);
        String[] listFiles = dir1.list();

        for (String string : listFiles) {

            if (string.endsWith(".t.xml")) {
                tandemFiles.add(inputDir + File.separator + string);
            }

        }

        // for every tandem file, convert it to mzid, call fdr and threshold, then combine all mzid files
        startTime = System.currentTimeMillis();
        for (int i = 0; i < tandemFiles.size(); i++) {
            try {
                String tandemFile = tandemFiles.get(i);
                

                String tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};

                mzidLib.init(tandemInput);

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        stopTime = System.currentTimeMillis();
        difference = stopTime - startTime;

        String s = "\n\nTandem2mzid time without invokeAll " + difference / 1000 + " Seconds";

        cores = Runtime.getRuntime().availableProcessors();
        System.out.println("cores: " + cores);
        startTime = System.currentTimeMillis();
        for (int i = 0; i < tandemFiles.size(); i++) {
            try {
                String tandemFile = tandemFiles.get(i);

                String tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_tandem.mzid";
                String[] tandemInput = {"Tandem2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-databaseFileFormatID", "MS:1001348", "-massSpecFileFormatID", "MS:1001062", "-idsStartAtZero", "false", "-proteinCodeRegex", "\\S+", "-compress", "false"};

                //mzidLib.init(tandemInput);
                MzidCallable task = new MzidCallable(tandemInput);
                collection.add(task);

            } catch (Exception e) {
                e.printStackTrace();
            }

        }
        executorPool = Executors.newFixedThreadPool(cores);
        executorPool.invokeAll(collection);
        executorPool.shutdown();
        stopTime = System.currentTimeMillis();
        difference = stopTime - startTime;
        s = s + "\n\nTandem2mzid time with invokeAll " + difference / 1000 + " Seconds";
        return s;

    }

    public String runOmssa() throws InterruptedException {

        List<String> omssaFiles = new ArrayList();

        MzIdentMLLib mzidLib = new MzIdentMLLib();
        File dir1 = new File(inputDir);
        String[] listFiles = dir1.list();

        for (String string : listFiles) {

            if (string.endsWith(".omx")) {
                omssaFiles.add(inputDir + File.separator + string);
            }

        }

        // for every tandem file, convert it to mzid, call fdr and threshold, then combine all mzid files
        startTime = System.currentTimeMillis();
        for (int i = 0; i < omssaFiles.size(); i++) {
            try {
                String tandemFile = omssaFiles.get(i);

                String tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_omssa.mzid";
                String[] tandemInput = {"Omssa2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};

                mzidLib.init(tandemInput);

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        stopTime = System.currentTimeMillis();
        difference = stopTime - startTime;

        String s = "\n\nomssaFiles time without invokeAll " + difference / 1000 + " Seconds";

        cores = Runtime.getRuntime().availableProcessors();
        System.out.println("cores: " + cores);
        startTime = System.currentTimeMillis();
        for (int i = 0; i < omssaFiles.size(); i++) {
            try {
                String tandemFile = omssaFiles.get(i);

                String tandemoutputFile = tandemFile.substring(0, tandemFile.lastIndexOf(".")) + "_omssa.mzid";
                String[] tandemInput = {"Omssa2mzid", tandemFile, tandemoutputFile, "-outputFragmentation", "false", "-decoyRegex", "REVERSED", "-compress", "false"};

                //mzidLib.init(tandemInput);
                MzidCallable task = new MzidCallable(tandemInput);
                collection.add(task);

            } catch (Exception e) {
                e.printStackTrace();
            }

        }

        executorPool = Executors.newFixedThreadPool(cores);
        executorPool.invokeAll(collection);
        executorPool.shutdown();

        stopTime = System.currentTimeMillis();
        difference = stopTime - startTime;

        s = s + "\n\nomssaFiles time with invokeAll " + difference / 1000 + " Seconds";

        return s;

    }

    public static void main(String args[]) throws InterruptedException {

        JMzIdentMLPerformance jMzIdentML = new JMzIdentMLPerformance(args[0], args[1]);

        String s = jMzIdentML.runXtandem() + jMzIdentML.runOmssa();

        System.out.println(s);
    }

}
