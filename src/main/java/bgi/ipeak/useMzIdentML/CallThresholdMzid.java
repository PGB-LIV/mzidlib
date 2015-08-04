package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.ThresholdMzid;

public class CallThresholdMzid {
	private String mzid2threshold;
	private String output_mzid;
	private String cvAccessionForScoreThreshold;
	private Double threshValue;
	
	private Boolean isPSMThreshold;
	private Boolean betterScoresAreLower;
	private Boolean deleteUnderThreshold;
	
	public static void main(String[] args) {
	
	}
	/**
	 * 
	 * @param mzid2threshold
	 * @param output_mzid
	 * @param cvAccessionForScoreThreshold
	 * @param threshValue
	 * @param isPSMThreshold
	 * @param betterScoreAreLower
	 * @param deleteUnderThreshold
	 */
	public CallThresholdMzid(String mzid2threshold,String output_mzid,String cvAccessionForScoreThreshold,Double threshValue,
			Boolean isPSMThreshold,Boolean betterScoreAreLower,Boolean deleteUnderThreshold) {
		this.mzid2threshold=mzid2threshold;
		this.output_mzid=output_mzid;
		this.threshValue=threshValue;
		this.cvAccessionForScoreThreshold=cvAccessionForScoreThreshold;
		this.isPSMThreshold=isPSMThreshold;
		this.betterScoresAreLower=betterScoreAreLower;
		this.deleteUnderThreshold=deleteUnderThreshold;
	}
	
	public void Use_mzidlib2threshold() {
		@SuppressWarnings("unused")
		ThresholdMzid thresholdMzid = new ThresholdMzid(mzid2threshold, output_mzid, isPSMThreshold, cvAccessionForScoreThreshold, threshValue, betterScoresAreLower, deleteUnderThreshold, "PDH");
	}
}
