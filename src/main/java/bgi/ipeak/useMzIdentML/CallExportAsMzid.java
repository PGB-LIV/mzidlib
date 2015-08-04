/**
 *
 */
package bgi.ipeak.useMzIdentML;

import bgi.ipeak.util.Properties;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.kohsuke.args4j.Option;
import uk.ac.liv.mzidlib.converters.Omssa2mzid;
import uk.ac.liv.mzidlib.converters.Tandem2mzid;


/**
 * @author Administrator
 *
 */
class ConvertLog extends Thread {

    InputStream is;

    ConvertLog(InputStream is) {
        this.is = is;
        start();
    }

    public void run() {
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String line = null;
            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }
        } catch (IOException e) {
            e.printStackTrace();

        }
    }
}

public class CallExportAsMzid {

    @Option(name = "-f", required = true, usage = "(required) The input search result to convertion to mzid.")
    private String input_result_file;
    @Option(name = "-s", required = true, usage = "(required) The type of search engine.[0:OMSSA,1:X!Tandem]")
    private int search_engineID;
    @Option(name = "-o", required = false, usage = "(optional) The output path and file name.")
    private String out_mzidfilepath;
    @Option(name = "-dr", required = true, usage = "(required) The decoy regrex.")
    private String decoy_regrex;
    @Option(name = "-accessionSplitRegex", required = false, usage = "(optional) Regular expression to split protein ID, surrounded by forward slashes e.g. \"/ /\"")
    private String accessionSplitRegex = "/ /";
    @Option(name = "-usermod", required = false, usage = "(optional) The user modification file.")
    private String umod;
    @Option(name = "-frag", required = false, usage = "Output the fragmentation ions.")
    private boolean output_fraginfor;
    @Option(name = "-df", required = false, usage = "(optional) The database file format, default MS:1001348(FASTA).")
    private String database_format;
    @Option(name = "-mf", required = false, usage = "(optional) The spectrum file format, default MS:1001062(MGF)")
    private String spectrum_fileformat;
    @Option(name = "-bc", required = false, usage = "The coordinates for spectra. Turn on for 0-based, mzML searched. Otherwise, off.")
    private boolean bc;
	//private String mods=getClass().getClassLoader().getResource("resources/mods.xml").getPath();;

    private String mod_file = Properties.getModFile_path();

    /**
     *
     * @param aFile
     * @param out_mzid_file
     * @param usermod
     * @param reg
     */
    public CallExportAsMzid(String aFile, String out_mzid_file, String usermod, String reg) {
        this.input_result_file = aFile;
        this.search_engineID = 0;
        this.out_mzidfilepath = out_mzid_file;        
        this.umod = usermod;
        this.decoy_regrex = reg;
        //this.mod_file=mod;
        this.output_fraginfor = false;
    }

    /**
     *
     * @param aFile
     * @param out_mzid_file
     * @param reg
     * @param df
     * @param accessionSplitRegex
     * @param mf
     * @param b0
     */
    public CallExportAsMzid(String aFile, String out_mzid_file, String reg, String df, String accessionSplitRegex, String mf, boolean b0) {
        this.input_result_file = aFile;
        this.search_engineID = 1;
        this.out_mzidfilepath = out_mzid_file;
        this.decoy_regrex = reg;
        this.output_fraginfor = false;
        this.database_format = df;
        this.accessionSplitRegex = accessionSplitRegex;
        this.spectrum_fileformat = mf;
        this.bc = b0;
    }

    /**
     *
     * @throws Exception
     */
    public void convert() throws Exception {
        if (search_engineID == 1) {
            Tandem2mzid tandem2mzid = new Tandem2mzid(input_result_file, out_mzidfilepath, database_format, spectrum_fileformat,
                    bc, decoy_regrex, accessionSplitRegex, output_fraginfor);
        } else if (search_engineID == 0) {
            //System.out.println("mod file="+mod_file+"\t"+input_result_file+"\t");
            Omssa2mzid omssa2mzid = new Omssa2mzid(input_result_file, 
                                                   out_mzidfilepath, 
                                                   output_fraginfor, 
                                                   decoy_regrex, 
                                                   mod_file, 
                                                   umod);
        }
    }

    public String getOut() {
        return out_mzidfilepath;
    }
}
