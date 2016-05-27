package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.converters.MzIdentMLToCSV;

public class CallMzIdentMLToCSV {

    private String mzid2trans;
    private String out_put_csv;
    private String verboseOutput;
    private String exportType;

    public CallMzIdentMLToCSV(String mzid2trans, String out_put_csv, String verboseOutput, String exportType) {
        this.mzid2trans = mzid2trans;
        this.out_put_csv = out_put_csv;
        this.verboseOutput = verboseOutput;
        this.exportType = exportType;
    }

    public void Use_mzidlib2trans() {
        MzIdentMLToCSV mzidToCsv = new MzIdentMLToCSV();
        mzidToCsv.useMzIdentMLToCSV(mzid2trans, out_put_csv, exportType, Boolean.valueOf(verboseOutput));
    }
}
