package uk.ac.liv.mzidlib.test;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import uk.ac.liv.mzidlib.MzIdentMLLib;

/**
 * 
 * @author Fawaz Ghali 28-Jul-2015 
 */


public class Performance {
    
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

    private boolean enablePercolator = true;
    private boolean enableMsgf = true;

    //private boolean deleteFiles = false;
    long startTime, stopTime, elapsedTime;

    private boolean verbose = true;
    
    public Performance(){
          mzidLib = new MzIdentMLLib();
    }
    
    // Run ProteoAnnotator
    public void runPerformance() {
        
    }

}
