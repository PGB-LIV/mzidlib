package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.fdr.FalseDiscoveryRateGlobal;


public class CallFalseDiscoveryRateGlobal {
	private String mzid;
	private String output_mzid;
	private String decoyregex;
	private String decoy_ratio;
	private String cvterm;
	private Boolean betterScoresAreLower;
	private String fdr_levle;
	private String protein_levle;
	public CallFalseDiscoveryRateGlobal() {
		
	}
	/**
	 * 
	 * @param mzid Mzid file.
	 * @param output_mzid Output mzid.
	 * @param decoyRatio Decoy ratio.
	 * @param decoy Decoy regex.
	 * @param cvTerm CV term.
	 * @param betterScoresAreLower Whether better scores are lower.
	 * @param fdrLevel FDR level.
	 * @param proteinLevel Protein level.
	 */
	public CallFalseDiscoveryRateGlobal(String mzid, String output_mzid, String decoyRatio, String decoy,
			String cvTerm, boolean betterScoresAreLower, String fdrLevel, String proteinLevel ) {
		this.mzid=mzid;
		this.output_mzid=output_mzid;
		this.decoyregex=decoy;
		this.decoy_ratio=decoyRatio;
		this.cvterm=cvTerm;
		this.betterScoresAreLower=betterScoresAreLower;
		this.fdr_levle=fdrLevel;
		this.protein_levle=proteinLevel;		
	}
	public void run_FdrAnalyse() {
        if (mzid != null && decoyregex != null && decoy_ratio != null && fdr_levle != null) {
            System.out.println("FalseDiscoveryRateGlobal");
            FalseDiscoveryRateGlobal fdrGlobal = new FalseDiscoveryRateGlobal(mzid, decoy_ratio, 
            		decoyregex, cvterm, betterScoresAreLower, fdr_levle, protein_levle);

            if (output_mzid != null) {
                fdrGlobal.computeFDRusingJonesMethod();
                fdrGlobal.writeToMzIdentMLFile(output_mzid);
            }
        }
	}
}
