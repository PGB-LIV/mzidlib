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
     * @param mzid2threshold Mzid File
     * @param output_mzid Output mzid file.
     * @param cvAccessionForScoreThreshold CV accession for score threshold.
     * @param threshValue Threshold value.
     * @param isPSMThreshold Whether to use PSM threshold.
     * @param betterScoreAreLower Whether better scores are lower.
     * @param deleteUnderThreshold Whether to delete under threshold.
     */
    public CallThresholdMzid(String mzid2threshold, String output_mzid, String cvAccessionForScoreThreshold, Double threshValue,
            Boolean isPSMThreshold, Boolean betterScoreAreLower, Boolean deleteUnderThreshold) {
        this.mzid2threshold = mzid2threshold;
        this.output_mzid = output_mzid;
        this.threshValue = threshValue;
        this.cvAccessionForScoreThreshold = cvAccessionForScoreThreshold;
        this.isPSMThreshold = isPSMThreshold;
        this.betterScoresAreLower = betterScoreAreLower;
        this.deleteUnderThreshold = deleteUnderThreshold;
    }

    public void Use_mzidlib2threshold() {
        @SuppressWarnings("unused")
        ThresholdMzid thresholdMzid = new ThresholdMzid(mzid2threshold, output_mzid, isPSMThreshold, cvAccessionForScoreThreshold, threshValue, betterScoresAreLower, deleteUnderThreshold, "PDH");
    }
}
