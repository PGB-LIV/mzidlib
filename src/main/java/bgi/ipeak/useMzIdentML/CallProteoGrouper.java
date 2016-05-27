package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.proteogrouper.ProteoGrouper;

/*
 import java.io.File;
 import java.io.FileWriter;
 import java.io.IOException;

 import org.kohsuke.rngom.parse.compact.UCode_UCodeESC_CharStream;
 import bgi.ipeak.io.StreamHog;
 import bgi.ipeak.util.Properties;
 */
public class CallProteoGrouper {

    private String mzid2group;
    private String output_mzid;
    private String cvAccForSIIScore;
    private Boolean requireSIIsToPassThreshold;
    private Boolean verboseOutput;
    private Boolean logTransScore;

    public static void main(String[] args) {
        // TODO 自动生成的方法存根

    }

    public CallProteoGrouper(String mzid2group, String output_mzid, Boolean requireSIIsToPassThreshold,
            Boolean verboseOutput, String cvAccForSIIScore, Boolean logTransScore) {
        // TODO 
        this.mzid2group = mzid2group;
        this.output_mzid = output_mzid;
        this.requireSIIsToPassThreshold = requireSIIsToPassThreshold;
        this.verboseOutput = verboseOutput;
        this.logTransScore = logTransScore;
        this.cvAccForSIIScore = cvAccForSIIScore;
    }

    @SuppressWarnings("unused")
    public void Use_mzidlib2ProteinGroup() {
        ProteoGrouper proteinInference = new ProteoGrouper(mzid2group, output_mzid,
                requireSIIsToPassThreshold, verboseOutput, cvAccForSIIScore, logTransScore, Boolean.valueOf("false"), Boolean.valueOf("false"));
    }
}
