package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.fdr.FalseDiscoveryRate;

public class CallFalseDiscoveryRate {

    private String mzid2fdr;
    private String output_mzid;
    private String decoyregex;
    private Integer decoyvalue;
    private String cvterm;
    private Boolean betterScoresAreLower;

    public static void main(String[] args) {

    }

    public CallFalseDiscoveryRate(String mzid2fdr, String output_mzid, String decoyregex,
            Integer decoyvalue, String cvterm, Boolean betterScoresAreLower) {
        this.mzid2fdr = mzid2fdr;
        this.output_mzid = output_mzid;
        this.decoyregex = decoyregex;
        this.decoyvalue = decoyvalue;
        this.cvterm = cvterm;
        this.betterScoresAreLower = betterScoresAreLower;
    }

    public void UseMzid2GetFalseDiscoverRate() {
        FalseDiscoveryRate fdr = new FalseDiscoveryRate(mzid2fdr,
                String.valueOf(decoyvalue), decoyregex, cvterm, betterScoresAreLower);

        fdr.computeFDRusingJonesMethod();
        fdr.writeToMzIdentMLFile(output_mzid);
    }
}
