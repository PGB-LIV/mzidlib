package bgi.ipeak.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;

public class PSMSummaryInformation {

    private String psm_id = "";
    private String spectrum_index = "";
    private String peptide_sequence = "";
    private Double mz = 0.0;
    private Double mz_error = 0.0;
    private Double theory_mz = 0.0;
    private Vector<String> mods = new Vector<String>();
    private HashMap<String, String> Scores = new HashMap<String, String>();
    private Vector<String> proteins = new Vector<String>();
    private Integer charge = 0;
    private Integer rank = 0;
    private Boolean is_decoy = false;
    private Boolean pass_threshold = false;
    private String main_threshold_score;
    private String main_scoreName = "";
    private String other_scoreName = "";

    public PSMSummaryInformation() {

    }

    public void print_psm_infor(BufferedWriter writed_file) throws IOException {
        String tem = "\t";
        String infor = psm_id + tem + spectrum_index + tem + peptide_sequence + tem + charge + tem + main_threshold_score + tem
                + Scores.values() + tem + mz + tem + mz_error + tem + theory_mz + tem + mods + tem + rank + tem + is_decoy + tem + pass_threshold + tem;
        String proString = "";
        for (int i = 0; i < proteins.size() - 1; i++) {
            proString += proteins.get(i) + ";";
        }
        proString += proteins.get(proteins.size() - 1);
        infor += proString;
        writed_file.write(infor);
        writed_file.newLine();
    }

    public void print_tile(BufferedWriter writed_file) throws IOException {
        String title = "PSMid\tSpectrum_index\tSequence\tCharge\t" + main_scoreName + "\t" + other_scoreName + "\tmz\tmz_error\ttheory_mz\t"
                + "Mods\trank\tis_decoy\tpass_threshold\tproteins";
        writed_file.write(title);
        writed_file.newLine();
    }

    public void setCharge(Integer charge) {
        this.charge = charge;
    }

    public void setIs_decoy(Boolean is_decoy) {
        this.is_decoy = is_decoy;
    }

    public void setMods(Vector<String> mods) {
        this.mods = mods;
    }

    public void setMz(Double mz) {
        this.mz = mz;
    }

    public void setMz_error(Double mz_error) {
        this.mz_error = mz_error;
    }

    public void setTheory_mz(Double theory_mz) {
        this.theory_mz = theory_mz;
    }

    public void setPass_threshold(Boolean pass_threshold) {
        this.pass_threshold = pass_threshold;
    }

    public void setPeptide_sequence(String peptide_sequence) {
        this.peptide_sequence = peptide_sequence;
    }

    public void setProteins(Vector<String> proteins) {
        this.proteins = proteins;
    }

    public void setPsm_id(String psm_id) {
        this.psm_id = psm_id;
    }

    public void setRank(Integer rank) {
        this.rank = rank;
    }

    public void setScores(HashMap<String, String> scores) {
        if (scores.isEmpty()) {
            other_scoreName = "NoOtherScore";
        }
        Scores = scores;
        other_scoreName = scores.keySet().toString();
    }

    public void setSpectrum_index(String spectrum_index) {
        this.spectrum_index = spectrum_index;
    }

    public void setMain_threshold_score(HashMap<String, String> main_threshold_score) {
        this.main_threshold_score = main_threshold_score.values().iterator().next();
        this.main_scoreName = main_threshold_score.keySet().iterator().next();
    }

    public Boolean getPass_threshold() {
        return pass_threshold;
    }

    public Boolean getIs_decoy() {
        return is_decoy;
    }

    public String getPeptide_sequence() {
        return peptide_sequence;
    }

    public Vector<String> getProteins() {
        return proteins;
    }

    public String getSpectrum_index() {
        return spectrum_index;
    }

    public String getMain_scoreName() {
        return main_scoreName;
    }

}
