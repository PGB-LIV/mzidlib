package bgi.ipeak.percolator;

import java.io.File;
import java.io.IOException;
import java.util.Set;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.percolator.paser.XTandemParser;
import bgi.ipeak.useMzIdentML.CallExportAsMzid;
import bgi.ipeak.util.Properties;

public class XtandemPercolator {

    @Option(name = "-f", required = true, usage = "(required) The input search result. Xtandem in xml")
    private String searchResult = "";
    @Option(name = "-out_dir", required = true, usage = "(optional) output directory")
    private String out_dir;
    //@Option(name="-pp",required=true,usage="(required) percolator path. ")
    //private String percolator_path="";	

    @Option(name = "-databaseformat", required = false, usage = "(optional) The database file format, default MS:1001348(FASTA).")
    private String databaseformat = "MS:1001348";
    @Option(name = "-spectrumformat", required = false, usage = "(optional) The spectrum file format, default MS:1001062(MGF)")
    private String spectrumformat = "MS:1001062";
    @Option(name = "-accessionSplitRegex", required = false, usage = "(optional) Regular expression to split protein ID, surrounded by forward slashes e.g. \"/ /\"")
    private String accessionSplitRegex = "/ /";
    @Option(name = "-bc", required = false, usage = "The coordinates for spectra. Turn on for 0-based, mzML searched. Otherwise, off.[false](default)")
    private boolean bc = false;
    @Option(name = "-decoyregex", required = false, usage = "(optional) Decoy regular expression, default ###\\w+###")
    private String decoyregex = "###\\w+###";
    @Option(name = "-ms2", required = false, usage = "(optional) use the MS/MS  features[use[true],not use[false](default)]")
    private Boolean ms2 = false;
    private Integer searchtype;
    private String feature_file;
    private String mzidfile;
    private File[] result_list;
    private Vector<String> features_file_list = new Vector<String>();
    private Vector<String> mzid_list = new Vector<String>();
    private Vector<String> mzid_addScore_list = new Vector<String>();

    public static void main(String[] args) throws Exception {
        XtandemPercolator ip = new XtandemPercolator();
        CmdLineParser parser = new CmdLineParser(ip);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp IPeak.jar bgi.ipeak.percolator.XtandemPercolator [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        ip.analyse();
    }

    //public XtandemPercolator(String input_file,String output_file,String decoyregex,String proteinSplitRegex, String percolator_path) {
    public XtandemPercolator(String input_file, String output_file, String decoyregex, String proteinSplitRegex) {
        searchResult = input_file;
        out_dir = output_file;
        this.decoyregex = decoyregex;
        //this.percolator_path=percolator_path;
        this.accessionSplitRegex = proteinSplitRegex;
        //Properties.Properties_set(this.percolator_path, "");
    }

    public XtandemPercolator() {
    }

    public void analyse() throws Exception {
        searchtype = check_filePath();
        trans_file2mzid();
        get_Features();
        runPercolator();
        add_scores2mzid();
        System.out.println("Xtandem Percolator Done!");
    }

    private void trans_file2mzid() throws Exception {
        if (searchtype == 0) {
            File result_file = new File(searchResult);
            String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
            mzidfile = out_dir + name + ".mzid";
            CallExportAsMzid xtandem2mzid = new CallExportAsMzid(result_file.getAbsolutePath(), mzidfile, decoyregex, databaseformat,
                    accessionSplitRegex, spectrumformat, bc);
            xtandem2mzid.convert();
        } else if (searchtype == 1) {
            for (File result_file : result_list) {
                String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
                String out_file = out_dir + name + ".mzid";
                CallExportAsMzid xtandem2mzid = new CallExportAsMzid(result_file.getAbsolutePath(), out_file, decoyregex,
                        databaseformat, accessionSplitRegex, spectrumformat, bc);
                xtandem2mzid.convert();
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
            XTandemParser paser = new XTandemParser(result_file.getAbsolutePath(), feature_file, this.decoyregex);
            paser.get_feature_file();
        } else if (searchtype == 1) {
            for (File result_file : result_list) {
                String name = result_file.getName().substring(0, result_file.getName().lastIndexOf("."));
                String feature = out_dir + name + ".features";
                XTandemParser paser = new XTandemParser(result_file.getAbsolutePath(), feature_file, decoyregex);
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

    public static void xtandem_percolator(String input_file, String output_dir, String decoyregex, String proteinSplitRegex) throws Exception {
        //XtandemPercolator run_percolator=new XtandemPercolator(input_file,output_dir,decoyregex,proteinSplitRegex,percolator_path);
        XtandemPercolator run_percolator = new XtandemPercolator(input_file, output_dir, decoyregex, proteinSplitRegex);
        run_percolator.analyse();
    }
}
