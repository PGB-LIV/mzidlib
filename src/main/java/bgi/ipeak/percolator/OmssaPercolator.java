package bgi.ipeak.percolator;

import java.io.File;
import java.io.IOException;
import java.util.Set;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.percolator.paser.OmssaParser;
import bgi.ipeak.useMzIdentML.CallExportAsMzid;
import bgi.ipeak.util.Properties;

public class OmssaPercolator {

    @Option(name = "-f", required = true, usage = "(required) The input search result. OMSSA in omx")
    private String searchResult = "";
    @Option(name = "-database", required = true, usage = "(required) database path")
    private String database;
	//@Option(name="-pp",required=true,usage="(required) percolator path. ")
    //private String percolator_path="";	
    //@Option(name="-mod",required=true,usage="(required) percolator path. ")
    private String mod_file = "";

    @Option(name = "-out_dir", required = false, usage = "(optional) output directory")
    private String out_dir;
    @Option(name = "-decoyregex", required = false, usage = "(optional) Decoy regular expression, default ###\\w+###")
    private String decoyregex = "###\\w+###";
    @Option(name = "-ms2", required = false, usage = "(optional) use the MS/MS  features[use[true],not use[false](default)]")
    private Boolean ms2 = false;

    private String usermodfile = "";
    private Integer searchtype;
    private String feature_file;
    private String mzidfile;
    private File[] result_list;
    private Vector<String> features_file_list = new Vector<String>();
    private Vector<String> mzid_list = new Vector<String>();
    private Vector<String> mzid_addScore_list = new Vector<String>();

    public static void main(String[] args) throws Exception {
        OmssaPercolator ip = new OmssaPercolator();
        CmdLineParser parser = new CmdLineParser(ip);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp IPeak.jar bgi.ipeak.percolator.OmssaPercolator [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        ip.analyse();
    }

    public OmssaPercolator(String input_file, String output_dir, String database, String decoyregex) {
        //public OmssaPercolator(String input_file,String output_dir,String mod_file,String database,String decoyregex, String percolator_path) {
        this.searchResult = input_file;
        this.out_dir = output_dir;
        this.decoyregex = decoyregex;
        
        this.database = database;
        this.mod_file = Properties.getModFile_path();
        //this.percolator_path=percolator_path;
        //Properties.Properties_set(;this.percolator_path, "");

    }

    public OmssaPercolator() {
        this.mod_file = Properties.getModFile_path();
    }

    public void analyse() throws Exception {
        searchtype = check_filePath();
        trans_file2mzid();
        get_Features();
        runPercolator();
        add_scores2mzid();
        System.out.println("Omssa Percolator Done!");
    }

    private void trans_file2mzid() throws Exception {
        if (searchtype == 0) {
            File omssa_file = new File(searchResult);
            String name = omssa_file.getName().substring(0, omssa_file.getName().lastIndexOf("."));
            mzidfile = out_dir + name + ".mzid";
            CallExportAsMzid omx2mzid = new CallExportAsMzid(searchResult, mzidfile, usermodfile, decoyregex);
            omx2mzid.convert();
        } else if (searchtype == 1) {
            for (File result_file : result_list) {
                String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
                String out_file = out_dir + name + ".mzid";
                CallExportAsMzid omx2mzid = new CallExportAsMzid(result_file.getAbsolutePath(), out_file, usermodfile, decoyregex);
                omx2mzid.convert();
                mzid_list.add(out_file);
            }
        }
        System.out.println("convert prosperity\n\n");
    }

    private void get_Features() throws Exception {
        if (searchtype == 0) {
            File result_file = new File(searchResult);
            String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
            feature_file = out_dir + name + ".features";
            OmssaParser paser = new OmssaParser(this.searchResult, feature_file, usermodfile, database, decoyregex);
            paser.get_feature_file();
        } else if (searchtype == 1) {
            for (File result_file : result_list) {
                String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
                String feature = out_dir + name + ".features";
                OmssaParser paser = new OmssaParser(result_file.getAbsolutePath(), feature, usermodfile, database, decoyregex);
                paser.get_feature_file();
                features_file_list.add(feature);
            }
        }
    }

    private void runPercolator() throws IOException, InterruptedException {
        if (searchtype == 0) {
            File result_file = new File(searchResult);
            String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
            String tempXML = out_dir + name + ".tempXML";
            RunPercolator run_percolator = new RunPercolator(feature_file, out_dir + name, tempXML, true);
            run_percolator.execute();
            IPeakPercolator.catTargetAndDecoy(out_dir + name + ".target.txt", out_dir + name + ".decoy.txt", out_dir + name + ".per.txt");
        } else if (searchtype == 1) {
            feature_file = out_dir + "combined.features";
            IPeakPercolator.combined_features(features_file_list, feature_file);
            String name = "combined";
            String tempXML = out_dir + name + ".tempXML";
            RunPercolator run_percolator = new RunPercolator(feature_file, out_dir + name, tempXML, true);
            run_percolator.execute();
            IPeakPercolator.catTargetAndDecoy(out_dir + name + ".target.txt", out_dir + name + ".decoy.txt", out_dir + name + ".per.txt");
        }
    }

    private void add_scores2mzid() throws IOException {
        if (searchtype == 0) {
            String pertxtfile = feature_file.substring(0, feature_file.lastIndexOf(".")) + ".per.txt";
            PlusPropertiesToMzid pop = new PlusPropertiesToMzid(mzidfile, pertxtfile);
            pop.export();
        } else if (searchtype == 1) {
            String pertxtfile = feature_file.substring(0, feature_file.lastIndexOf(".")) + ".per.txt";
            Set<String> per_name_list = IPeakPercolator.split_percolatorResult(pertxtfile);
            for (String per_name : per_name_list) {
                String mzid_name = per_name.substring(0, per_name.indexOf(".per.txt")) + ".mzid";
                PlusPropertiesToMzid pop = new PlusPropertiesToMzid(mzid_name, per_name);
                pop.export();
            }
        }
    }

    private int check_filePath() {
        File out_file = new File(out_dir);
        if (out_file != null && !out_file.exists()) {
            out_file.mkdir();
        }
        out_dir = out_file.getAbsolutePath() + File.separator;
        File result_file = new File(searchResult);
        int type = -1;

        if (result_file.isFile()) {
            type = 0;
        } else if (result_file.isDirectory()) {
            result_list = result_file.listFiles();
            type = 1;
        } else {
            type = 2;
        }
        return type;
    }

    public static void omssa_percolator(String input_file, String output_dir, String database, String decoyregex) throws Exception {
        //OmssaPercolator run_percolator=new OmssaPercolator(input_file,output_dir,mod_file,database,decoyregex,percolator_path);
        OmssaPercolator run_percolator = new OmssaPercolator(input_file, output_dir, database, decoyregex);
        run_percolator.analyse();
    }
}
