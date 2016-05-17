package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.multiplesearch.CombineSearchEngines;

public class CallCombineSearchEngines {
    // the files add percolator score: 

    private String msgf_mzidaddPS;
    private String omssa_mzidaddPS;
    private String xtandem_mzidaddPS;
    private String output_file;
    private String debug_file;
    private String decoy_regex;
    private int Score_type;

    public static void main(String[] args) {

    }

    public CallCombineSearchEngines(String msgf_mzid, String omssa_mzid, String xtandem_mzid, String output_file,
            String decoy_regex, String debug_file, int Score_type) {
        this.msgf_mzidaddPS = msgf_mzid;
        this.xtandem_mzidaddPS = xtandem_mzid;
        this.omssa_mzidaddPS = omssa_mzid;
        this.output_file = output_file;
        this.debug_file = debug_file;
        this.decoy_regex = decoy_regex;
        this.Score_type = Score_type;
    }

    public void Comined_Search() {
        if (this.Score_type == 1) {
            try {
                Use_MzidlibeCombineEvalue();
            } catch (Exception e) {

                e.printStackTrace();
            }
        } else {
            try {
                Use_MzidlibeCombine();
            } catch (Exception e) {

                e.printStackTrace();
            }
        }
    }

    private void Use_MzidlibeCombine() throws Exception {
//		String parames="java  -jar "+Properties.getMzidLib_path()+" CombineSearchEngines -firstFile  "+this.omssa_mzidaddPS+ " -firstSearchEngine s1 "
//				+"-firstcvTerm MS:1001492  -firstbetterScoresAreLower false -secondFile "+this.xtandem_mzidaddPS
//				+" -secondSearchEngine s2 -secondcvTerm MS:1001492 -secondbetterScoresAreLower false -thirdFile "
//				+this.msgf_mzidaddPS + " -thirdSearchEngine s3 -thirdcvTerm MS:1001492 -thirdbetterScoresAreLower false "
//				+ "-rank 1 -decoyRatio 1 -outputFile " + this.output_file +" -debugFile " +this.debug_file 
//				+" -decoyRegex " + this.decoy_regex  + " -compress false";
        String[] newArgs = new String[17];
        newArgs[0] = this.omssa_mzidaddPS;
        newArgs[1] = "s1";
        newArgs[2] = "MS:1001492";
        newArgs[3] = "false";
        newArgs[4] = this.xtandem_mzidaddPS;
        newArgs[5] = "s2";
        newArgs[6] = "MS:1001492";
        newArgs[7] = "false";
        newArgs[8] = this.msgf_mzidaddPS;
        newArgs[9] = "s3";
        newArgs[10] = "MS:1001492";;
        newArgs[11] = "false";
        newArgs[12] = "1";
        newArgs[13] = "1";
        newArgs[14] = output_file;
        newArgs[15] = this.debug_file;
        newArgs[16] = this.decoy_regex;
        CombineSearchEngines.runThreeSearchEngines(newArgs);
    }

    private void Use_MzidlibeCombineEvalue() throws Exception {
//		String parames="java  -jar "+Properties.getMzidLib_path()+" CombineSearchEngines -firstFile  "+this.omssa_mzidaddPS+ " -firstSearchEngine s1 "
//				+"-firstcvTerm MS:1001328  -firstbetterScoresAreLower true -secondFile "+this.xtandem_mzidaddPS
//				+" -secondSearchEngine s2 -secondcvTerm MS:1001330 -secondbetterScoresAreLower true -thirdFile "
//				+this.msgf_mzidaddPS + " -thirdSearchEngine s3 -thirdcvTerm MS:1002053 -thirdbetterScoresAreLower true "
//				+ "-rank 1 -decoyRatio 1 -outputFile " + this.output_file +" -debugFile " +this.debug_file 
//				+" -decoyRegex " + this.decoy_regex  + " -compress false";
        String[] newArgs = new String[17];
        newArgs[0] = this.omssa_mzidaddPS;
        newArgs[1] = "s1";
        newArgs[2] = "MS:1001328";
        newArgs[3] = "true";
        newArgs[4] = this.xtandem_mzidaddPS;
        newArgs[5] = "s2";
        newArgs[6] = "MS:1001330";
        newArgs[7] = "true";
        newArgs[8] = this.msgf_mzidaddPS;
        newArgs[9] = "s3";
        newArgs[10] = "MS:1002053";
        newArgs[11] = "true";
        newArgs[12] = "1";
        newArgs[13] = "1";
        newArgs[14] = output_file;
        newArgs[15] = this.debug_file;
        newArgs[16] = this.decoy_regex;
        CombineSearchEngines.runThreeSearchEngines(newArgs);
    }
}
