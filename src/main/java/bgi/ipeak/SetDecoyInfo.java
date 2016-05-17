package bgi.ipeak;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.HashMap;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.MzIdentML;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class SetDecoyInfo {

    @Option(name = "-f", required = true, usage = "(required) input mzid file path")
    String mzid_path;
    @Option(name = "-o", required = true, usage = "(required) output mzid file path")
    String out_path;
    @Option(name = "-dr", required = false, usage = "(options) decoy proteins tag")
    String decoyRegex = "###REV###";
    private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;

    public static void main(String[] args) throws FileNotFoundException {

        SetDecoyInfo set = new SetDecoyInfo();
        CmdLineParser parser = new CmdLineParser(set);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-11)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.CombineMzidFiles [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        set.SetAndPrint();
    }

    public SetDecoyInfo(String mzid_file_path, String out_path, String decoyRegex) {

        this.mzid_path = mzid_file_path;
        this.out_path = out_path;
        this.decoyRegex = decoyRegex;
    }

    public SetDecoyInfo() {

    }

    public void SetAndPrint() throws FileNotFoundException {
        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid_path));
        mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
        SequenceCollection aCollection = mzIdentML.getSequenceCollection();
        HashMap<String, DBSequence> id_dbseq = new HashMap<String, DBSequence>();
        for (DBSequence db : aCollection.getDBSequence()) {
            id_dbseq.put(db.getId(), db);
        }
        List<PeptideEvidence> iterPeptideEvidence = aCollection.getPeptideEvidence();
        for (PeptideEvidence peptideEvidence : iterPeptideEvidence) {
//			System.out.println("db_ref:"+peptideEvidence.getDBSequenceRef());			
            String db_acc = id_dbseq.get(peptideEvidence.getDBSequenceRef()).getAccession();
//			System.out.println("db_acc:"+db_acc);
            if (db_acc.contains(decoyRegex)) {
                peptideEvidence.setIsDecoy(true);
            } else {
                peptideEvidence.setIsDecoy(false);
            }
        }

        MzIdentMLMarshaller m = new MzIdentMLMarshaller();

        m.marshal(mzIdentML, new FileOutputStream(out_path));
    }
}
