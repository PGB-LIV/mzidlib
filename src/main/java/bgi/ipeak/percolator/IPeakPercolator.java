/**
 *
 */
package bgi.ipeak.percolator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.percolator.paser.MSGFParser;
import bgi.ipeak.percolator.paser.OmssaParser;
import bgi.ipeak.percolator.paser.XTandemParser;
import bgi.ipeak.util.Properties;

/**
 * @author Administrator
 *
 */
public class IPeakPercolator {

    @Option(name = "-s", required = true, usage = "(required) The search engine code.\n\t\t[0:features,1:MSGF+(default),2:OMSSA,3:X!Tandem]\n\t\tif -s 0,the option -feature is needed")
    private int searchEngine = 1;
    @Option(name = "-r", required = false, usage = "(option) The input search result. MSGF+ in mzid, OMSSA in omx and X!Tandem in pepXml")
    private String searchResult = "";
    @Option(name = "-out", required = false, usage = "(optional) Path and file name for output.")
    private String out;
    @Option(name = "-decoyregex", required = false, usage = "(optional) Decoy regular expression, default ###\\w+###")
    private String decoyregex = "###\\w+###";
    @Option(name = "-ms2", required = false, usage = "(optional) use the MS/MS  features[use[true],not use[false](default)]")
    private Boolean ms2 = false;
    @Option(name = "-feature", required = false, usage = "(optional) features file")
    private String feature_file = "";
    @Option(name = "-database", required = true, usage = "(required) database path")
    private String database;
    private String modfile = Properties.getModFile_path();
    private String usermodfile = "";

    /**
     * Constructors for IPeakPercolator.
     *
     * @param searchenginID Search engine ID.
     * @param searchResult Search result.
     * @param output_name Output name.
     */
    public IPeakPercolator(int searchenginID, String searchResult, String output_name) {
        this.searchEngine = searchenginID;
        this.searchResult = searchResult;
        this.out = output_name;
    }

    /**
     *
     * @param searchenginID Search engine ID.
     * @param searchResult Search result.
     * @param output_name Output name.
     * @param decoyregex Decoy regex.
     */
    public IPeakPercolator(int searchenginID, String searchResult, String output_name, String decoyregex) {
        this.searchEngine = searchenginID;
        this.searchResult = searchResult;
        this.decoyregex = decoyregex;
        this.out = output_name;
    }

    /**
     *
     * @param searchenginID Search engine ID.
     * @param searchResult Search result.
     * @param output_name Output name.
     * @param Database The database.
     * @param usermod_fiel User mod field.
     * @param decoyregex Decoy regex.
     */
    public IPeakPercolator(int searchenginID, String searchResult, String output_name, String Database,
            String usermod_fiel, String decoyregex) {
        this.searchEngine = searchenginID;
        this.searchResult = searchResult;
        this.decoyregex = decoyregex;
        this.out = output_name;
        this.database = Database;
        //this.modfile=mod_file;
        this.usermodfile = usermod_fiel;
    }

    public IPeakPercolator(String feature_file, String out_file) throws Exception {
        this.feature_file = feature_file;
        this.out = out_file;
    }

    public IPeakPercolator() {

    }

    /**
     * @param args Arguments to the main method.
     * @throws Exception Exception thrown by the main method.
     */
    public static void main(String[] args) throws Exception {
        IPeakPercolator ip = new IPeakPercolator();
        CmdLineParser parser = new CmdLineParser(ip);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Guilin Li (liguilin@genomics.org.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.org.cn.ipeak.core.IPeakPercolator [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        ip.get_Features();
        ip.runPercolator();
    }

    public void get_Features() throws Exception {
        if (this.out == null) {
            this.out = this.searchResult.substring(0, this.searchResult.lastIndexOf("."));
        }
        feature_file = this.out + ".features";
        if (this.searchEngine == 1) {
            MSGFParser.setUse_ms2feature(ms2);
            MSGFParser paser = new MSGFParser(this.searchResult, feature_file, this.decoyregex);
            paser.get_feature_file();
        } else if (this.searchEngine == 2) {
            OmssaParser paser = new OmssaParser(this.searchResult, feature_file, usermodfile, database, decoyregex);
            paser.get_feature_file();
        } else if (this.searchEngine == 3) {
            XTandemParser paser = new XTandemParser(this.searchResult, feature_file, this.decoyregex);
            paser.get_feature_file();
        } else if (this.searchEngine != 0 || feature_file.equals("")) {
            System.err.println("the optins -s 0, and -feature is need");
        }
    }

    public void runPercolator() throws Exception {
        String tempXML = this.out + ".temp.XML";
        try {
            RunPercolator run_percolator = new RunPercolator(feature_file, this.out, tempXML, true);
            run_percolator.execute();
            IPeakPercolator.catTargetAndDecoy(this.out + ".target.txt", this.out + ".decoy.txt", this.out + ".per.txt");
        } catch (IOException e) {
            throw e;
        }
    }

    /**
     * Cat percolator target result and decoy result to one file "*.per.txt".
     *
     * @param target Target result.
     * @param decoy Decoy result.
     * @param per Concatenated result.
     */
    public static void catTargetAndDecoy(String target, String decoy, String per) {
        try {
            BufferedReader fr = new BufferedReader(new FileReader(new File(target)));
            FileWriter fw = new FileWriter(new File(per));
            String line = fr.readLine();
            String[] head = line.split("\t");
            String newHead = head[0] + "\tlabel";
            for (int i = 1; i < head.length; i++) {
                newHead += "\t" + head[i];
            }
            fw.append(newHead + "\n");
            fw.flush();
            while ((line = fr.readLine()) != null) {
                String[] items = line.split("\t");
                String newLine = items[0] + "\t1";
                for (int i = 1; i < items.length; i++) {
                    newLine += "\t" + items[i];
                }
                fw.append(newLine + "\n");
                fw.flush();
            }
            BufferedReader bf = new BufferedReader(new FileReader(new File(decoy)));
            line = bf.readLine();
            while ((line = bf.readLine()) != null) {
                String[] items = line.split("\t");
                String newLine = items[0] + "\t-1";
                for (int i = 1; i < items.length; i++) {
                    newLine += "\t" + items[i];
                }
                fw.append(newLine + "\n");
                fw.flush();
            }
            fr.close();
            bf.close();
            fw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void combined_features(Vector<String> feature_file_list, String all_feature_file) throws IOException {

        File test_file = new File(all_feature_file);
        if (test_file.exists()) {
            test_file.delete();
        }
        BufferedReader readtitle = new BufferedReader(new FileReader(feature_file_list.get(0)));
        String title = readtitle.readLine();
        readtitle.close();
        Vector<String> featureinfors = new Vector<String>();
        for (String feature_file : feature_file_list) {
            System.out.println(feature_file);
            BufferedReader read = new BufferedReader(new FileReader(feature_file));
            String line = read.readLine();
            while ((line = read.readLine()) != null) {
                featureinfors.add(line);
            }
            read.close();
        }

        BufferedWriter writeFile = new BufferedWriter(new FileWriter(all_feature_file));
        writeFile.write(title);
        for (String infor : featureinfors) {
            writeFile.append(infor);
            writeFile.newLine();
        }
        writeFile.close();
    }

    public static Set<String> split_percolatorResult(String pertxtfile) throws IOException {
        File per_file = new File(pertxtfile);
        String out_dir = per_file.getParent() + "/";
        BufferedReader read = new BufferedReader(new FileReader(pertxtfile));
        String line = read.readLine();
        String title = line;
        HashMap<String, Vector<String>> per_file_list = new HashMap<String, Vector<String>>();
        while ((line = read.readLine()) != null) {
            String[] items = line.split(":");
            String id = out_dir + items[0] + ".per.txt";
            if (per_file_list.containsKey(id)) {
                Vector<String> infor = per_file_list.get(id);
                infor.add(line);
                per_file_list.put(id, infor);
            } else {
                Vector<String> infor = new Vector<String>();
                infor.add(line);
                per_file_list.put(id, infor);
            }
        }
        read.close();

        Set<String> out_names = per_file_list.keySet();
        for (String out : out_names) {
            BufferedWriter writer = new BufferedWriter(new FileWriter(out));
            writer.append(title);
            writer.newLine();
            for (String line_infor : per_file_list.get(out)) {
                writer.append(line_infor);
                writer.newLine();
            }
            writer.close();
        }
        return out_names;
    }

    public static Vector<String> addPercolator2mzid(Set<String> per_name, String mziddir, String out_dir) throws IOException {
        if (mziddir.lastIndexOf("/") != mziddir.length() - 1) {
            mziddir = mziddir + "/";
        }
        if (out_dir.lastIndexOf("/") != out_dir.length() - 1) {
            out_dir = out_dir + "/";
        }
        String tag = ".mzid";
        String tag_out = "AddP.mzid";

        Vector<String> addp_mzdi_list = new Vector<String>();
        for (String file_name : per_name) {
            File per_file = new File(file_name);
            String name = per_file.getName().substring(0, per_file.getName().lastIndexOf(".per.txt"));
            String mzid_path = mziddir + name + tag;
            String mzidP_path = out_dir + name + tag_out;
            File mz_file = new File(mzid_path);
            if (mz_file.exists()) {
                PlusPropertiesToMzid pmp = new PlusPropertiesToMzid(mzid_path, file_name, mzidP_path);
                pmp.export();
            }
            addp_mzdi_list.add(mzidP_path);
        }
        return addp_mzdi_list;
    }

    public String getFeature_file() {
        File file = new File(feature_file);
        return file.getAbsolutePath();
    }
}
