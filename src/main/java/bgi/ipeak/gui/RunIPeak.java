package bgi.ipeak.gui;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

import javax.swing.JOptionPane;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

import bgi.ipeak.CombineMzidFiles;
import bgi.ipeak.Mzid2Summary;
import bgi.ipeak.SetDecoyInfo;
import bgi.ipeak.io.WritePSM2Summary;
import bgi.ipeak.io.WriteProteins2Summary;
import bgi.ipeak.percolator.PlusPropertiesToMzid;
import bgi.ipeak.percolator.IPeakPercolator;
import bgi.ipeak.test.ThresholdPSM;
import bgi.ipeak.useMzIdentML.CallCombineSearchEngines;
import bgi.ipeak.useMzIdentML.CallExportAsMzid;
import bgi.ipeak.useMzIdentML.CallFalseDiscoveryRate;
import bgi.ipeak.useMzIdentML.CallFalseDiscoveryRateGlobal;
import bgi.ipeak.useMzIdentML.CallMzIdentMLToCSV;
import bgi.ipeak.useMzIdentML.CallInsertMetaDataFromFasta;
import bgi.ipeak.useMzIdentML.CallProteoGrouper;
import bgi.ipeak.useMzIdentML.CallThresholdMzid;
import bgi.ipeak.util.IPeakFileList;
import bgi.ipeak.util.Properties;

public class RunIPeak extends SwingWorker<Void, Void> {

    public static JTextArea log = null;
    private String msgf;
    private String omssa;
    private String xtandem;
    private String database;
    private String out_file_regex;
    private String out_file_dir;
    private String percolator_path;
    private String decoyregrex = "###REV###";
    private String mod = "";
    private String umod = "";
    private String databaseformat = "MS:1001348";
    private String spectrumformat = "MS:1001062";
    private boolean bc = false;
    private String accessionSplitRegex = "/ /";
    private String debug = "false";
    private Integer analyse_step = 0;
    private String use_peptidefdr="true";
    private Double maxfdr = 0.01;
    private Boolean delete = false;
    private Integer searchtype = 0;
    private Integer group_type = 0;
    private String threshold_befor_combine = "false";
    private String thresholdScoreAcc = "MS:1001491";
    private Double thresholdScoreValue = 0.2;
    private IPeakFileList files_list;
    
    private File[] msgf_list;
    private File[] xtandem_list;
    private File[] omssa_list;

    /**
     * @param log the log to set
     */
    public void setLog(JTextArea log) {
        RunIPeak.log = log;
    }

    /**
     * @return the log
     */
    public JTextArea getLog() {
        return log;
    }

    public void setMsgf(String msgf) {
        this.msgf = msgf;
    }

    public void setUse_peptidefdr(String use_peptidefdr) {
		this.use_peptidefdr = use_peptidefdr;
	}
    public void setOmssa(String omssa) {
        this.omssa = omssa;
    }

    public void setXtandem(String xtandem) {
        this.xtandem = xtandem;
    }

    public void setDatabase(String database) {
        this.database = database;
    }

    public void setSpectrum_file_regex(String out_file_regex) {
        this.out_file_regex = out_file_regex;
    }

    public void setOut_file_dir(String out_file_dir) {
        this.out_file_dir = out_file_dir;
    }

    public void setDecoyregrex(String decoyregrex) {
        this.decoyregrex = decoyregrex;
    }

    public void setMod(String mod) {
        this.mod = mod;
    }

    public void setUmod(String umod) {
        this.umod = umod;
    }

    public void setAccessionSplitRegex(String accessionSplitRegex) {
        this.accessionSplitRegex = accessionSplitRegex;
    }

    public void setDebug(String debug) {
        this.debug = debug;
    }

    public void setMaxfdr(Double maxfdr) {
        this.maxfdr = maxfdr;
    }

    public void setDelete(Boolean delete) {
        this.delete = delete;
    }

    public void setBc(boolean bc) {
        this.bc = bc;
    }

    public void setPercolator_path(String percolator_path) {
        this.percolator_path = percolator_path;
    }

    private int set_filelist_and_get_analyse_method() {
        log.append("Check the percolator ...\n");
        log.setCaretPosition(log.getDocument().getLength());

        File msgf_file = new File(msgf);
        File xtandem_file = new File(xtandem);
        File omssa_file = new File(omssa);
        int type = -1;

        if (msgf_file.isFile() && xtandem_file.isFile() && omssa_file.isFile()) {
            //System.out.println(out_file_dir);
            files_list = new IPeakFileList(msgf, omssa, xtandem, out_file_dir, out_file_regex);
            type = 0;
        } else if (msgf_file.isDirectory() && xtandem_file.isDirectory() && omssa_file.isDirectory()) {
            msgf_list = msgf_file.listFiles();
            xtandem_list = xtandem_file.listFiles();
            omssa_list = omssa_file.listFiles();
            files_list = new IPeakFileList(msgf_file.getAbsolutePath() + "all.mzid", omssa_file.getAbsolutePath()
                    + "all.omx", xtandem_file.getAbsolutePath() + "all.xml", out_file_dir, out_file_regex);
            type = 1;
        } else {
            type = 2;
        }
        return type;
    }

    private void basic_analyse() throws Exception {

        if (!isCancelled()) {
            basic_analyse_step_1();
        } else {
            return;
        }

        if (!isCancelled()) {
            add_percolatorScore2Mzid();
        } else {
            return;
        }

        if (searchtype == 1) {
            if (!isCancelled()) {
                CombinedMzid();
            } else {
                return;
            }
        }

        if (!isCancelled()) {
            basic_analyse_step_2();
        } else {
            return;
        }
    }

    private void basic_analyse_step_1() throws Exception {

        if (!isCancelled()) {
            trans_file2mzid();
        } else {
            return;
        }

        if (!isCancelled()) {
            percolator_analyse();
        } else {
            return;
        }
    }

    private void basic_analyse_step_2() throws Exception {

        if (!isCancelled()) {
            AddProtein2Mzid();
        } else {
            return;
        }

        if (!isCancelled()) {
            Combined_ThreeMethod();
            Get_PeptideLevalFDR();
        } else {
            return;
        }

        if (group_type == 1) {
            if (!isCancelled()) {
                Threshold_Result();
            } else {
                return;
            }

            if (!isCancelled()) {
                ProteinGroup();
            } else {
                return;
            }

            if (!isCancelled()) {
                Mzid2csv();
            } else {
                return;
            }

            if (!isCancelled()) {
                trans2summary_report();
            } else {
                return;
            }
        } else {
            if (!isCancelled()) {
                Mzid2Summary();
            } else {
                return;
            }
        }
    }

    private void single_step_run() throws Exception {
        if (analyse_step == 1) {
            basic_analyse_step_1();
            System.exit(0);
        }
        if (group_type == 0) {
            if (analyse_step == 2) {
                basic_analyse_step_2();
            } else if (analyse_step == 3) {
                Combined_ThreeMethod();
                Mzid2Summary();
            } else if (analyse_step > 3) {
                Mzid2Summary();
            }
        } else {
            if (analyse_step == 2) {
                basic_analyse_step_2();
                System.exit(0);
            } else if (analyse_step == 3) {
                Combined_ThreeMethod();
                Threshold_Result();
                ProteinGroup();
                Mzid2csv();
                trans2summary_report();
                System.exit(0);
            } else if (analyse_step == 4) {
                Threshold_Result();
                ProteinGroup();
                Mzid2csv();
                trans2summary_report();
                System.exit(0);
            } else if (analyse_step == 5) {
                Mzid2csv();
                trans2summary_report();
                System.exit(0);
            } else if (analyse_step == 6) {
                trans2summary_report();
                System.exit(0);
            }
        }
    }

    private void trans_file2mzid() throws Exception {

        log.append("Convert *.omx/*.xml to mzid format file(s)...\n");
        log.setCaretPosition(log.getDocument().getLength());
        if (searchtype == 0) {
            String xtandem_mzidfile = files_list.getXtandem_mzid();
            String omssa_mzidfile = files_list.getOmssa_mzid();
            String msgf_mzidfile = files_list.getMsgf_mzid();
            CallExportAsMzid omx2mzid = new CallExportAsMzid(omssa, omssa_mzidfile,  umod, decoyregrex);
            CallExportAsMzid xtd2mzid = new CallExportAsMzid(xtandem, xtandem_mzidfile, decoyregrex, databaseformat,accessionSplitRegex, spectrumformat, bc);
            SetDecoyInfo set_decoy = new SetDecoyInfo(msgf, msgf_mzidfile, decoyregrex);
            set_decoy.SetAndPrint();
            omx2mzid.convert();
            xtd2mzid.convert();
        } else {
            for (File xtandem_file : xtandem_list) {
                String out_dir = files_list.getXtandem_out_dir();
                String name = xtandem_file.getName().substring(0, xtandem_file.getName().lastIndexOf("."));
                CallExportAsMzid xtd2mzid = new CallExportAsMzid(xtandem_file.getAbsolutePath(),
                        out_dir + name + ".mzid", decoyregrex, databaseformat, accessionSplitRegex,spectrumformat, bc);
                xtd2mzid.convert();
            }
            for (File omssa_file : omssa_list) {
                String out_dir = files_list.getOmssa_out_dir();
                String name = omssa_file.getName().substring(0, omssa_file.getName().lastIndexOf("."));
                CallExportAsMzid omx2mzid = new CallExportAsMzid(omssa_file.getAbsolutePath(),
                        out_dir + name + ".mzid",  umod, decoyregrex);
                omx2mzid.convert();
            }
            for (File msgf_file : msgf_list) {
                String out_dir = files_list.getMsgf_out_dir();
                String name = msgf_file.getName().substring(0, msgf_file.getName().lastIndexOf("."));
                SetDecoyInfo set_decoy = new SetDecoyInfo(msgf_file.getAbsolutePath(), out_dir + name + ".mzid", decoyregrex);
                set_decoy.SetAndPrint();
            }
        }
        log.append("Format conversion done!\n");
        log.setCaretPosition(log.getDocument().getLength());
    }

    private void percolator_analyse() throws Exception {
        log.append("Run percolator...\n");
        log.setCaretPosition(log.getDocument().getLength());
        if (searchtype == 0) {
            String xtandem_primaryfilename = files_list.getXtandem_primaryfilename();
            String omssa_primaryfilename = files_list.getOmssa_primaryfilename();
            String msgf_primaryfilename = files_list.getMsgf_primaryfilename();
            log.append("\nCollect features for X!Tandem...");
            IPeakPercolator xp = new IPeakPercolator(3, xtandem, xtandem_primaryfilename, decoyregrex);
            xp.get_Features();
            log.append("\nCollect features for OMSSA...");
            IPeakPercolator op = new IPeakPercolator(2, omssa, omssa_primaryfilename, database,  umod, decoyregrex);
            op.get_Features();
            log.append("\nCollect features for MS-GF+...");
            IPeakPercolator mp = new IPeakPercolator(1, msgf, msgf_primaryfilename, decoyregrex);
            mp.get_Features();
            log.append("Feature collection done\n\n");
            if (analyse_step != 1) {
                log.append("\nRun percolator for X!Tandem...");
                xp.runPercolator();
                log.append("\nRun percolator for MS-GF+...");
                mp.runPercolator();
                log.append("\nRun percolator for OMSSA...");
                op.runPercolator();
                log.append("Percolator running done\n\n");
            }
            log.setCaretPosition(log.getDocument().getLength());
        } else {

            Vector<String> xtandem_feature_files = new Vector<String>();
            Vector<String> omssa_feature_files = new Vector<String>();
            Vector<String> msgf_feature_files = new Vector<String>();
            xtandem_feature_files = get_featureFromSplitResult("xtandem", xtandem_list);
            omssa_feature_files = get_featureFromSplitResult("omssa", omssa_list);
            msgf_feature_files = get_featureFromSplitResult("msgf", msgf_list);
            log.append("Collect features done!\n\n");

            String xtandem_all_feature_file = files_list.getXtandem_feature_file();
            String omssa_all_feature_file = files_list.getOmssa_feature_file();
            String msgf_all_feature_file = files_list.getMsgf_feature_file();
            combined_features(xtandem_feature_files, xtandem_all_feature_file);
            combined_features(omssa_feature_files, omssa_all_feature_file);
            combined_features(msgf_feature_files, msgf_all_feature_file);
            log.append("Merging Feature files done\n\n");
            log.setCaretPosition(log.getDocument().getLength());

            IPeakPercolator xtandem_percolator = new IPeakPercolator(xtandem_all_feature_file,
                    xtandem_all_feature_file.substring(0, xtandem_all_feature_file.lastIndexOf(".")));
            xtandem_percolator.runPercolator();
            IPeakPercolator omssa_percolator = new IPeakPercolator(omssa_all_feature_file,
                    omssa_all_feature_file.substring(0, omssa_all_feature_file.lastIndexOf(".")));
            omssa_percolator.runPercolator();
            IPeakPercolator msgf_percolator = new IPeakPercolator(msgf_all_feature_file,
                    msgf_all_feature_file.substring(0, msgf_all_feature_file.lastIndexOf(".")));
            msgf_percolator.runPercolator();

            String xtandem_per = xtandem_all_feature_file.substring(0, xtandem_all_feature_file.lastIndexOf(".")) + ".per.txt";
            String omssa_per = omssa_all_feature_file.substring(0, omssa_all_feature_file.lastIndexOf(".")) + ".per.txt";
            String msgf_per = msgf_all_feature_file.substring(0, msgf_all_feature_file.lastIndexOf(".")) + ".per.txt";
            Set<String> xtandem_pertxt = split_percolatorResult(xtandem_per);
            Set<String> omssa_pertxt = split_percolatorResult(omssa_per);
            Set<String> msgf_pertxt = split_percolatorResult(msgf_per);
            files_list.setMsgf_per_list(msgf_pertxt);
            files_list.setOmssa_per_list(omssa_pertxt);
            files_list.setXtandem_per_list(xtandem_pertxt);
        }
    }

    private void add_percolatorScore2Mzid() throws InterruptedException, IOException {
        log.append("Adding percolator scores to mzid file...\n");
        if (searchtype == 0) {
            String omssa_plused_file = files_list.getOmssa_mzidAddPer();
            String msgf_plused_file = files_list.getMsgf_mzidAddPer();
            String xtandem_plused_file = files_list.getXtandem_mzidAddPer();
            PlusPropertiesToMzid pmp = new PlusPropertiesToMzid(files_list.getMsgf_mzid(), files_list.getMsgf_percolator_result_file(), msgf_plused_file);
            PlusPropertiesToMzid pop = new PlusPropertiesToMzid(files_list.getOmssa_mzid(), files_list.getOmssa_percolator_result_file(), omssa_plused_file);
            PlusPropertiesToMzid pxp = new PlusPropertiesToMzid(files_list.getXtandem_mzid(), files_list.getXtandem_percolator_result_file(), xtandem_plused_file);
            pmp.export();
            pop.export();
            pxp.export();
            log.append("Adding percolator scores done\n\n");
        } else {
            Vector<String> msgf_addPmzid = addPercolator2mzid(files_list.getMsgf_per_list(),
                    files_list.getMsgf_out_dir(), files_list.getMsgf_out_dir());
            Vector<String> omssa_addPmzid = addPercolator2mzid(files_list.getOmssa_per_list(),
                    files_list.getOmssa_out_dir(), files_list.getOmssa_out_dir());
            Vector<String> xtandem_addPmzid = addPercolator2mzid(files_list.getXtandem_per_list(),
                    files_list.getXtandem_out_dir(), files_list.getXtandem_out_dir());
            files_list.setMsgf_addP_list(msgf_addPmzid);
            files_list.setOmssa_addP_list(omssa_addPmzid);
            files_list.setXtandem_addP_list(xtandem_addPmzid);
        }
        log.setCaretPosition(log.getDocument().getLength());
    }

    private void CombinedMzid() {
        Vector<String> list = files_list.getMsgf_addP_list();
        Integer kInteger = 1;
        HashMap<String, Integer> name_id = new HashMap<String, Integer>();
        for (String j : list) {
            File file = new File(j);
            String name = file.getName();
            name_id.put(name, kInteger);
            kInteger++;
        }
        CombineMzidFiles combine_msgf = new CombineMzidFiles(files_list.getMsgf_addP_list(), files_list.getMsgf_mzidAddPer(), name_id);
        combine_msgf.combined_mzid();
        CombineMzidFiles combine_omssa = new CombineMzidFiles(files_list.getOmssa_addP_list(), files_list.getOmssa_mzidAddPer(), name_id);
        combine_omssa.combined_mzid();
        CombineMzidFiles combine_xtandem = new CombineMzidFiles(files_list.getXtandem_addP_list(), files_list.getXtandem_mzidAddPer(), name_id);
        combine_xtandem.combined_mzid();
    }

    private void AddProtein2Mzid() throws InterruptedException, IOException {

        log.append("Merge proteins from OMSSA/MSGF+/X!Tandem\n");
        System.out.println("accessionSplitRegex:" + accessionSplitRegex);
        CallInsertMetaDataFromFasta plus_pro_msgf = new CallInsertMetaDataFromFasta(files_list.getMsgf_mzidAddPer(), files_list.getMsgf_mzidAddProPer(), database, accessionSplitRegex);
        plus_pro_msgf.AddFasta_UseMzidlib();
        CallInsertMetaDataFromFasta plus_pro_omssa = new CallInsertMetaDataFromFasta(files_list.getOmssa_mzidAddPer(), files_list.getOmssa_mzidAddProPer(), database, accessionSplitRegex);
        plus_pro_omssa.AddFasta_UseMzidlib();
        CallInsertMetaDataFromFasta plus_pro_xtandem = new CallInsertMetaDataFromFasta(files_list.getXtandem_mzidAddPer(), files_list.getXtandem_mzidAddProPer(), database, accessionSplitRegex);
        plus_pro_xtandem.AddFasta_UseMzidlib();
        log.append("Adding Proteins information done\n\n");
        log.setCaretPosition(log.getDocument().getLength());
    }

    private void Combined_ThreeMethod() throws IOException, InterruptedException {

        log.append("Combine percolator results...\n");

        String omssa_plused_file = files_list.getOmssa_mzidAddProPer();
        String msgf_plused_file = files_list.getMsgf_mzidAddProPer();
        String xtandem_plused_file = files_list.getXtandem_mzidAddProPer();
        if (Boolean.valueOf(threshold_befor_combine)) {
            ThresholdPSM thresholdOmssa = new ThresholdPSM(omssa_plused_file, omssa_plused_file + "T", thresholdScoreValue, thresholdScoreAcc);
            thresholdOmssa.ThresholdAndPrint();
            ThresholdPSM thresholdMsgf = new ThresholdPSM(msgf_plused_file, msgf_plused_file + "T", thresholdScoreValue, thresholdScoreAcc);
            thresholdMsgf.ThresholdAndPrint();
            ThresholdPSM thresholdXtandem = new ThresholdPSM(xtandem_plused_file, xtandem_plused_file + "T", thresholdScoreValue, thresholdScoreAcc);
            thresholdXtandem.ThresholdAndPrint();
            omssa_plused_file = files_list.getOmssa_mzidAddProPer() + "T";
            msgf_plused_file = files_list.getMsgf_mzidAddProPer() + "T";
            xtandem_plused_file = files_list.getXtandem_mzidAddProPer() + "T";
        }

        String mxo_combined_fiel = files_list.getCombined_resultUsePerScore_primaryfilename();
        String debug_file = files_list.getCombined_result_debugfile();
        CallCombineSearchEngines combine = new CallCombineSearchEngines(msgf_plused_file, omssa_plused_file, xtandem_plused_file,
                mxo_combined_fiel, decoyregrex, debug_file, 0);
        combine.Comined_Search();
        log.append("Combining results done\n\n");
        log.setCaretPosition(log.getDocument().getLength());

        if (Boolean.valueOf(debug)) {

            String Emxo_combined_fiel = files_list.getEcombined_resultUseEvalue_primaryfilename();
            String Edebug_file = files_list.getEcombined_result_debugfile();
            CallCombineSearchEngines Ecombine = new CallCombineSearchEngines(msgf_plused_file, omssa_plused_file, xtandem_plused_file,
                    Emxo_combined_fiel, decoyregrex, Edebug_file, 1);
            Ecombine.Comined_Search();
            log.append(new Date() + " DEBUG Combine result prosperity\n\n");
            log.setCaretPosition(log.getDocument().getLength());

            String mzid2fdr_m = files_list.getMsgf_mzidAddProPer();
            String mzid2fdr_x = files_list.getXtandem_mzidAddProPer();
            String mzid2fdr_o = files_list.getOmssa_mzidAddProPer();

            String mzid2threshold_me = files_list.getMsgfPP_efdr();
            String mzid2threshold_xe = files_list.getXtandemPP_efdr();
            String mzid2threshold_oe = files_list.getOmssaPP_efdr();
            String mzid2threshold_mp = files_list.getMsgfPP_pfdr();
            String mzid2threshold_xp = files_list.getXtandemPP_pfdr();
            String mzid2threshold_op = files_list.getOmssaPP_pfdr();

            CallFalseDiscoveryRate get_fdr_me = new CallFalseDiscoveryRate(mzid2fdr_m, mzid2threshold_me, decoyregrex, 1, "MS:1002053", true);
            CallFalseDiscoveryRate get_fdr_xe = new CallFalseDiscoveryRate(mzid2fdr_x, mzid2threshold_xe, decoyregrex, 1, "MS:1001330", true);
            CallFalseDiscoveryRate get_fdr_oe = new CallFalseDiscoveryRate(mzid2fdr_o, mzid2threshold_oe, decoyregrex, 1, "MS:1001328", true);
            CallFalseDiscoveryRate get_fdr_mp = new CallFalseDiscoveryRate(mzid2fdr_m, mzid2threshold_mp, decoyregrex, 1, "MS:1001492", false);
            CallFalseDiscoveryRate get_fdr_xp = new CallFalseDiscoveryRate(mzid2fdr_x, mzid2threshold_xp, decoyregrex, 1, "MS:1001492", false);
            CallFalseDiscoveryRate get_fdr_op = new CallFalseDiscoveryRate(mzid2fdr_o, mzid2threshold_op, decoyregrex, 1, "MS:1001492", false);
            get_fdr_me.UseMzid2GetFalseDiscoverRate();
            get_fdr_xe.UseMzid2GetFalseDiscoverRate();
            get_fdr_oe.UseMzid2GetFalseDiscoverRate();
            get_fdr_mp.UseMzid2GetFalseDiscoverRate();
            get_fdr_xp.UseMzid2GetFalseDiscoverRate();
            get_fdr_op.UseMzid2GetFalseDiscoverRate();
        }
    }

   /**
     * 
     */
    private void Get_PeptideLevalFDR(){
    	/*
    	 * peptide leval fdr options;
    	 */
    	String decoy_rate="1";
    	String FDRScoreCV="MS:1002355";
    	String CombinedFDRScore="MS:1002356";
    	Boolean better_score_are_lower=true;
    	String fdr_leval="Peptide";
    	String protein_leval="PAG";
    	
        String mxo_percolator = files_list.getCombined_resultUsePerScore();
        String outpep_mxo_percolator=files_list.getPepFDR_combined_resultUsePerScore();
        System.out.println("Combine result prosperity\n\n");
        CallFalseDiscoveryRateGlobal mxo_pep =new CallFalseDiscoveryRateGlobal(mxo_percolator,outpep_mxo_percolator,
        		decoy_rate,decoyregrex,CombinedFDRScore,better_score_are_lower,fdr_leval,protein_leval);
        mxo_pep.run_FdrAnalyse();
        if (Boolean.valueOf(debug)) {
            String mxo_evalue = files_list.getEcombined_resultUseEvalue_primaryfilename();
            String outpep_mxo_evalue = files_list.getPepFDR_Ecombined_resultUseEvalue();
            CallFalseDiscoveryRateGlobal mxo_evalue_pep =new CallFalseDiscoveryRateGlobal(mxo_evalue,outpep_mxo_evalue,
            		decoy_rate,decoyregrex,CombinedFDRScore,better_score_are_lower,fdr_leval,protein_leval);
            mxo_evalue_pep.run_FdrAnalyse();
            System.out.println("Combine result prosperity\n\n");

            String mzid2pepFDR_me = files_list.getMsgfPP_efdr();
            String mzid2pepFDR_xe = files_list.getXtandemPP_efdr();
            String mzid2pepFDR_oe = files_list.getOmssaPP_efdr();
            String mzid2pepFDR_mp = files_list.getMsgfPP_pfdr();
            String mzid2pepFDR_xp = files_list.getXtandemPP_pfdr();
            String mzid2pepFDR_op = files_list.getOmssaPP_pfdr();
            String outpep_mzid2pepFDR_me = files_list.getPepFDR_msgfPP_efdr();
            String outpep_mzid2pepFDR_xe = files_list.getPepFDR_xtandemPP_efdr();
            String outpep_mzid2pepFDR_oe = files_list.getPepFDR_omssaPP_efdr();
            String outpep_mzid2pepFDR_mp = files_list.getPepFDR_msgfPP_pfdr();
            String outpep_mzid2pepFDR_xp = files_list.getPepFDR_xtandemPP_pfdr();
            String outpep_mzid2pepFDR_op = files_list.getPepFDR_omssaPP_pfdr();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_me =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_me,outpep_mzid2pepFDR_me,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_me.run_FdrAnalyse();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_xe =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_xe,outpep_mzid2pepFDR_xe,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_xe.run_FdrAnalyse();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_oe =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_oe,outpep_mzid2pepFDR_oe,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_oe.run_FdrAnalyse();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_mp =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_mp,outpep_mzid2pepFDR_mp,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_mp.run_FdrAnalyse();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_xp =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_xp,outpep_mzid2pepFDR_xp,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_xp.run_FdrAnalyse();
            CallFalseDiscoveryRateGlobal getpep_mzid2pepFDR_op =new CallFalseDiscoveryRateGlobal(mzid2pepFDR_op,outpep_mzid2pepFDR_op,
            		decoy_rate,decoyregrex,FDRScoreCV,better_score_are_lower,fdr_leval,protein_leval);
            getpep_mzid2pepFDR_op.run_FdrAnalyse();
            
        }
    }
    /**
     * 
     * @throws IOException
     * @throws InterruptedException
     */
    private void Mzid2Summary() throws IOException, InterruptedException {  
    	String Combined_FDRScoreCV="MS:1002356";
    	String FDRScoreCV="MS:1002355";
    	if(Boolean.valueOf(use_peptidefdr)){
    		Combined_FDRScoreCV="MS:1002360";
    		FDRScoreCV="MS:1002360";
    	}    	
        String mzid2threshold = files_list.getPepFDR_combined_resultUsePerScore();
        String output_summary = files_list.getCombined_resultUsePerScore_primaryfilename();
        Mzid2Summary mxo_mzid2Summary = new Mzid2Summary(mzid2threshold, output_summary, database,
                maxfdr, true, Combined_FDRScoreCV, decoyregrex, true);
        mxo_mzid2Summary.trans_mzid2txt();
        Integer mxo_protein_number = mxo_mzid2Summary.getProtein_number();
        Integer mxo_peptide_number = mxo_mzid2Summary.getPeptide_number();
        Integer mxo_psm_number = mxo_mzid2Summary.getPSM_number();
        System.out.println("trans combined percolator result file to summary result prosperity\n\n");

        FileWriter writeFile = new FileWriter(files_list.getStatistic_analysis_file());
        writeFile.write("the number of psm,peptide,protein\n");
        writeFile.write("\tPSM\tPeptide\tProtein\n");
        writeFile.write("IPeak\t" + mxo_psm_number + "\t" + mxo_peptide_number + "\t" + mxo_protein_number + "\n");

        if (Boolean.valueOf(debug)) {

            String Emzid2threshold = files_list.getPepFDR_Ecombined_resultUseEvalue();
            String Eoutput_summary = files_list.getEcombined_resultUseEvalue_primaryfilename();
            Mzid2Summary Emxo_mzid2Summary = new Mzid2Summary(Emzid2threshold, Eoutput_summary, database,
                    maxfdr, true, Combined_FDRScoreCV, decoyregrex, true);
            Emxo_mzid2Summary.trans_mzid2txt();
            Integer Emxo_protein_number = Emxo_mzid2Summary.getProtein_number();
            Integer Emxo_peptide_number = Emxo_mzid2Summary.getPeptide_number();
            Integer Emxo_psm_number = Emxo_mzid2Summary.getPSM_number();
            System.out.println("trans combined evalue result file to summary result prosperity\n\n");

            String mzid2summary_me = files_list.getPepFDR_msgfPP_efdr();
            String mzid2summary_xe = files_list.getPepFDR_xtandemPP_efdr();
            String mzid2summary_oe = files_list.getPepFDR_omssaPP_efdr();
            String mzid2summary_mp = files_list.getPepFDR_msgfPP_pfdr();
            String mzid2summary_xp = files_list.getPepFDR_xtandemPP_pfdr();
            String mzid2summary_op = files_list.getPepFDR_omssaPP_pfdr();
            String summary_me = files_list.getMsgf_primaryfilename() + "efdr_";
            String summary_xe = files_list.getXtandem_primaryfilename() + "efdr";
            String summary_oe = files_list.getOmssa_primaryfilename() + "efdr";
            String summary_mp = files_list.getMsgf_primaryfilename() + "pfdr";
            String summary_xp = files_list.getXtandem_primaryfilename() + "pfdr";
            String summary_op = files_list.getOmssa_primaryfilename() + "pfdr";
            Mzid2Summary me_Mzid2Summary = new Mzid2Summary(mzid2summary_me, summary_me, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            Mzid2Summary xe_Mzid2Summary = new Mzid2Summary(mzid2summary_xe, summary_xe, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            Mzid2Summary oe_Mzid2Summary = new Mzid2Summary(mzid2summary_oe, summary_oe, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            Mzid2Summary mp_Mzid2Summary = new Mzid2Summary(mzid2summary_mp, summary_mp, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            Mzid2Summary xp_Mzid2Summary = new Mzid2Summary(mzid2summary_xp, summary_xp, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            Mzid2Summary op_Mzid2Summary = new Mzid2Summary(mzid2summary_op, summary_op, database, maxfdr, true, FDRScoreCV, decoyregrex, true);
            me_Mzid2Summary.trans_mzid2txt();
            Integer me_protein_number = me_Mzid2Summary.getProtein_number();
            Integer me_peptide_number = me_Mzid2Summary.getPeptide_number();
            Integer me_psm_number = me_Mzid2Summary.getPSM_number();
            System.out.println("trans msgf evalue result file to summary result prosperity\n\n");

            mp_Mzid2Summary.trans_mzid2txt();
            Integer mp_protein_number = mp_Mzid2Summary.getProtein_number();
            Integer mp_peptide_number = mp_Mzid2Summary.getPeptide_number();
            Integer mp_psm_number = mp_Mzid2Summary.getPSM_number();
            System.out.println("trans msgf percolator result file to summary result prosperity\n\n");

            xe_Mzid2Summary.trans_mzid2txt();
            Integer xe_protein_number = xe_Mzid2Summary.getProtein_number();
            Integer xe_peptide_number = xe_Mzid2Summary.getPeptide_number();
            Integer xe_psm_number = xe_Mzid2Summary.getPSM_number();
            System.out.println("trans xtandem evalue result file to summary result prosperity\n\n");

            xp_Mzid2Summary.trans_mzid2txt();
            Integer xp_protein_number = xp_Mzid2Summary.getProtein_number();
            Integer xp_peptide_number = xp_Mzid2Summary.getPeptide_number();
            Integer xp_psm_number = xp_Mzid2Summary.getPSM_number();
            System.out.println("trans xtandem percolator result file to summary result prosperity\n\n");

            oe_Mzid2Summary.trans_mzid2txt();
            Integer oe_protein_number = oe_Mzid2Summary.getProtein_number();
            Integer oe_peptide_number = oe_Mzid2Summary.getPeptide_number();
            Integer oe_psm_number = oe_Mzid2Summary.getPSM_number();
            System.out.println("trans omssa evalue result file to summary result prosperity\n\n");

            op_Mzid2Summary.trans_mzid2txt();
            Integer op_protein_number = op_Mzid2Summary.getProtein_number();
            Integer op_peptide_number = op_Mzid2Summary.getPeptide_number();
            Integer op_psm_number = op_Mzid2Summary.getPSM_number();
            System.out.println("trans omssa percolaotor result file to summary result prosperity\n\n");

            writeFile.write("FDRAnalyse\t" + Emxo_psm_number + "\t" + Emxo_peptide_number + "\t" + Emxo_protein_number + "\n");
            writeFile.write("MP\t" + mp_psm_number + "\t" + mp_peptide_number + "\t" + mp_protein_number + "\n");
            writeFile.write("MS-GF+\t" + me_psm_number + "\t" + me_peptide_number + "\t" + me_protein_number + "\n");
            writeFile.write("OP\t" + op_psm_number + "\t" + op_peptide_number + "\t" + op_protein_number + "\n");
            writeFile.write("Omssa\t" + oe_psm_number + "\t" + oe_peptide_number + "\t" + oe_protein_number + "\n");
            writeFile.write("XP\t" + xp_psm_number + "\t" + xp_peptide_number + "\t" + xp_protein_number + "\n");
            writeFile.write("X!Tandem\t" + xe_psm_number + "\t" + xe_peptide_number + "\t" + xe_protein_number + "\n");
        }
        writeFile.close();
    }
    private void Threshold_Result() throws IOException, InterruptedException {
        log.append("Filtering...\n");
        String mzid2threshold = files_list.getCombined_resultUsePerScore();
        String output_mzid_t = files_list.getCombined_threshold_file();
        CallThresholdMzid mxo_threshold = new CallThresholdMzid(mzid2threshold, output_mzid_t, "MS:1002356", maxfdr, true, true, delete);
        mxo_threshold.Use_mzidlib2threshold();
        log.append("threshold the combined result prosperity\n\n");
        log.setCaretPosition(log.getDocument().getLength());

        if (Boolean.valueOf(debug)) {

            String Emzid2threshold = files_list.getEcombined_resultUseEvalue();
            String Eoutput_mzid_t = files_list.getEcombined_threshold_file();
            CallThresholdMzid Emxo_threshold = new CallThresholdMzid(Emzid2threshold, Eoutput_mzid_t, "MS:1002356", maxfdr, true, true, delete);
            Emxo_threshold.Use_mzidlib2threshold();
            log.append(new Date() + " DEBUG threshold the combined result prosperity\n\n");
            log.setCaretPosition(log.getDocument().getLength());

            String mzid2threshold_me = files_list.getMsgfPP_efdr();
            String mzid2threshold_xe = files_list.getXtandemPP_efdr();
            String mzid2threshold_oe = files_list.getOmssaPP_efdr();
            String mzid2threshold_mp = files_list.getMsgfPP_pfdr();
            String mzid2threshold_xp = files_list.getXtandemPP_pfdr();
            String mzid2threshold_op = files_list.getOmssaPP_pfdr();
            String mzid2group_me = files_list.getMsgfPP_efdr_threshold();
            String mzid2group_xe = files_list.getXtandemPP_efdr_threshold();
            String mzid2group_oe = files_list.getOmssaPP_efdr_threshold();
            String mzid2group_mp = files_list.getMsgfPP_pfdr_threshold();
            String mzid2group_xp = files_list.getXtandemPP_pfdr_threshold();
            String mzid2group_op = files_list.getOmssaPP_pfdr_threshold();
            CallThresholdMzid me_threshold = new CallThresholdMzid(mzid2threshold_me, mzid2group_me, "MS:1001874", maxfdr, true, true, delete);
            CallThresholdMzid xe_threshold = new CallThresholdMzid(mzid2threshold_xe, mzid2group_xe, "MS:1001874", maxfdr, true, true, delete);
            CallThresholdMzid oe_threshold = new CallThresholdMzid(mzid2threshold_oe, mzid2group_oe, "MS:1001874", maxfdr, true, true, delete);
            CallThresholdMzid mp_threshold = new CallThresholdMzid(mzid2threshold_mp, mzid2group_mp, "MS:1001874", maxfdr, true, true, delete);
            CallThresholdMzid xp_threshold = new CallThresholdMzid(mzid2threshold_xp, mzid2group_xp, "MS:1001874", maxfdr, true, true, delete);
            CallThresholdMzid op_threshold = new CallThresholdMzid(mzid2threshold_op, mzid2group_op, "MS:1001874", maxfdr, true, true, delete);
            mp_threshold.Use_mzidlib2threshold();
            xp_threshold.Use_mzidlib2threshold();
            op_threshold.Use_mzidlib2threshold();
            me_threshold.Use_mzidlib2threshold();
            xe_threshold.Use_mzidlib2threshold();
            oe_threshold.Use_mzidlib2threshold();
        }

    }

    private void ProteinGroup() throws IOException, InterruptedException {

        log.append("Protein Group Mergement...\n");
        String mzid2group = files_list.getCombined_threshold_file();
        String output_mzid_g = files_list.getGroup_combined_threshold_file();
        CallProteoGrouper mxo_pro_group = new CallProteoGrouper(mzid2group, output_mzid_g, true, false, "MS:1002356", true);
        mxo_pro_group.Use_mzidlib2ProteinGroup();
        log.append("group proteins prosperity\n\n");
        log.setCaretPosition(log.getDocument().getLength());

        if (Boolean.valueOf(debug)) {

            String Emzid2group = files_list.getEcombined_threshold_file();
            String Eoutput_mzid_g = files_list.getEgroup_combined_threshold_file();
            CallProteoGrouper Emxo_pro_group = new CallProteoGrouper(Emzid2group, Eoutput_mzid_g, true, false, "MS:1002356", true);
            Emxo_pro_group.Use_mzidlib2ProteinGroup();
            log.append(new Date() + " DEBUG group proteins prosperity\n\n");
            log.setCaretPosition(log.getDocument().getLength());

            String mzid2group_me = files_list.getMsgfPP_efdr_threshold();
            String mzid2group_xe = files_list.getXtandemPP_efdr_threshold();
            String mzid2group_oe = files_list.getOmssaPP_efdr_threshold();
            String mzid2group_mp = files_list.getMsgfPP_pfdr_threshold();
            String mzid2group_xp = files_list.getXtandemPP_pfdr_threshold();
            String mzid2group_op = files_list.getOmssaPP_pfdr_threshold();
            String grouped_mzid_me = files_list.getMsgfPP_efdr_threshold_group();
            String grouped_mzid_xe = files_list.getXtandemPP_efdr_threshold_group();
            String grouped_mzid_oe = files_list.getOmssaPP_efdr_threshold_group();
            String grouped_mzid_mp = files_list.getMsgfPP_pfdr_threshold_group();
            String grouped_mzid_xp = files_list.getXtandemPP_pfdr_threshold_group();
            String grouped_mzid_op = files_list.getOmssaPP_pfdr_threshold_group();
            CallProteoGrouper me_pro_group = new CallProteoGrouper(mzid2group_me, grouped_mzid_me, true, false, "MS:1001874", true);
            CallProteoGrouper xe_pro_group = new CallProteoGrouper(mzid2group_xe, grouped_mzid_xe, true, false, "MS:1001874", true);
            CallProteoGrouper oe_pro_group = new CallProteoGrouper(mzid2group_oe, grouped_mzid_oe, true, false, "MS:1001874", true);
            CallProteoGrouper mp_pro_group = new CallProteoGrouper(mzid2group_mp, grouped_mzid_mp, true, false, "MS:1001874", true);
            CallProteoGrouper xp_pro_group = new CallProteoGrouper(mzid2group_xp, grouped_mzid_xp, true, false, "MS:1001874", true);
            CallProteoGrouper op_pro_group = new CallProteoGrouper(mzid2group_op, grouped_mzid_op, true, false, "MS:1001874", true);
            me_pro_group.Use_mzidlib2ProteinGroup();
            xe_pro_group.Use_mzidlib2ProteinGroup();
            oe_pro_group.Use_mzidlib2ProteinGroup();
            mp_pro_group.Use_mzidlib2ProteinGroup();
            xp_pro_group.Use_mzidlib2ProteinGroup();
            op_pro_group.Use_mzidlib2ProteinGroup();
            log.append(new Date() + " DEBUG group proteins prosperity\n\n");
            log.setCaretPosition(log.getDocument().getLength());
        }
    }

    private void Mzid2csv() throws IOException, InterruptedException {

        log.append("Merge results and output...\n");
        String mzid2trans = files_list.getGroup_combined_threshold_file();
        String out_put_csv_psm = files_list.getPsm_combined_threshold_csvfile();
        String out_put_csv_pro = files_list.getProtein_combined_threshold_csvfile();
        CallMzIdentMLToCSV mxo_trans_psm = new CallMzIdentMLToCSV(mzid2trans, out_put_csv_psm, "true", "exportPSMs");
        CallMzIdentMLToCSV mxo_trans_pro = new CallMzIdentMLToCSV(mzid2trans, out_put_csv_pro, "true", "exportProteinGroups");
        mxo_trans_psm.Use_mzidlib2trans();
        mxo_trans_pro.Use_mzidlib2trans();
        log.append("trans to csv  prosperity\n\n");
        log.append("iPeak result file :\ncombined mzid file: " + mzid2trans + "\npsm csv file: " + out_put_csv_psm + "\n"
                + "protein group file: " + out_put_csv_pro + "\n");
        log.setCaretPosition(log.getDocument().getLength());

        if (Boolean.valueOf(debug)) {

            String e_mzid2trans = files_list.getEgroup_combined_threshold_file();
            String e_out_put_csv_psm = files_list.getEpsm_combined_threshold_csvfile();
            String e_out_put_csv_pro = files_list.getEprotein_combined_threshold_csvfile();
            CallMzIdentMLToCSV e_mxo_trans_psm = new CallMzIdentMLToCSV(e_mzid2trans, e_out_put_csv_psm, "true", "exportPSMs");
            CallMzIdentMLToCSV e_mxo_trans_pro = new CallMzIdentMLToCSV(e_mzid2trans, e_out_put_csv_pro, "true", "exportProteinGroups");
            e_mxo_trans_psm.Use_mzidlib2trans();
            e_mxo_trans_pro.Use_mzidlib2trans();
            log.append(new Date() + " DEBUG trans to csv  prosperity\n\n");
            log.append(new Date() + " DEBUG iPeak result file :\ncombined mzid file: " + mzid2trans + "\npsm csv file: " + out_put_csv_psm + "\n"
                    + "protein group file: " + out_put_csv_pro + "\n");
            log.setCaretPosition(log.getDocument().getLength());

            String grouped_mzid_me = files_list.getMsgfPP_efdr_threshold_group();
            String grouped_mzid_xe = files_list.getXtandemPP_efdr_threshold_group();
            String grouped_mzid_oe = files_list.getOmssaPP_efdr_threshold_group();
            String grouped_mzid_mp = files_list.getMsgfPP_pfdr_threshold_group();
            String grouped_mzid_xp = files_list.getXtandemPP_pfdr_threshold_group();
            String grouped_mzid_op = files_list.getOmssaPP_pfdr_threshold_group();
            String csv_psm_me = files_list.getMsgfPP_efdr_threshold_group_psmcsv();
            String csv_psm_xe = files_list.getXtandemPP_efdr_threshold_group_psmcsv();
            String csv_psm_oe = files_list.getOmssaPP_efdr_threshold_group_psmcsv();
            String csv_psm_mp = files_list.getMsgfPP_pfdr_threshold_group_psmcsv();
            String csv_psm_xp = files_list.getXtandemPP_pfdr_threshold_group_psmcsv();
            String csv_psm_op = files_list.getOmssaPP_pfdr_threshold_group_psmcsv();
            String csv_pro_me = files_list.getMsgfPP_efdr_threshold_group_procsv();
            String csv_pro_xe = files_list.getXtandemPP_efdr_threshold_group_procsv();
            String csv_pro_oe = files_list.getOmssaPP_efdr_threshold_group_procsv();
            String csv_pro_mp = files_list.getMsgfPP_pfdr_threshold_group_procsv();
            String csv_pro_xp = files_list.getXtandemPP_pfdr_threshold_group_procsv();
            String csv_pro_op = files_list.getOmssaPP_pfdr_threshold_group_procsv();
            CallMzIdentMLToCSV me_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_me, csv_psm_me, "true", "exportPSMs");
            CallMzIdentMLToCSV xe_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_xe, csv_psm_xe, "true", "exportPSMs");
            CallMzIdentMLToCSV oe_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_oe, csv_psm_oe, "true", "exportPSMs");
            CallMzIdentMLToCSV mp_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_mp, csv_psm_mp, "true", "exportPSMs");
            CallMzIdentMLToCSV xp_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_xp, csv_psm_xp, "true", "exportPSMs");
            CallMzIdentMLToCSV op_trans_psm = new CallMzIdentMLToCSV(grouped_mzid_op, csv_psm_op, "true", "exportPSMs");
            CallMzIdentMLToCSV me_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_me, csv_pro_me, "true", "exportProteinGroups");
            CallMzIdentMLToCSV xe_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_xe, csv_pro_xe, "true", "exportProteinGroups");
            CallMzIdentMLToCSV oe_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_oe, csv_pro_oe, "true", "exportProteinGroups");
            CallMzIdentMLToCSV mp_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_mp, csv_pro_mp, "true", "exportProteinGroups");
            CallMzIdentMLToCSV xp_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_xp, csv_pro_xp, "true", "exportProteinGroups");
            CallMzIdentMLToCSV op_trans_pro = new CallMzIdentMLToCSV(grouped_mzid_op, csv_pro_op, "true", "exportProteinGroups");
            me_trans_psm.Use_mzidlib2trans();
            xe_trans_psm.Use_mzidlib2trans();
            oe_trans_psm.Use_mzidlib2trans();
            mp_trans_psm.Use_mzidlib2trans();
            op_trans_psm.Use_mzidlib2trans();
            xp_trans_psm.Use_mzidlib2trans();
            me_trans_pro.Use_mzidlib2trans();
            xe_trans_pro.Use_mzidlib2trans();
            oe_trans_pro.Use_mzidlib2trans();
            mp_trans_pro.Use_mzidlib2trans();
            xp_trans_pro.Use_mzidlib2trans();
            op_trans_pro.Use_mzidlib2trans();
        }
    }

    private void trans2summary_report() throws IOException {

        log.append("Merge results and summarizing...\n");
        log.setCaretPosition(log.getDocument().getLength());
        WriteProteins2Summary mxo_percolator_pro = new WriteProteins2Summary(files_list.getProtein_combined_threshold_csvfile(), files_list.getProtein_combined_threshold_summaryfile());
        mxo_percolator_pro.transcsv2summary();
        Integer mxo_Percolator_protein_number = mxo_percolator_pro.getProtein_number();
        WritePSM2Summary mxo_percolator_psm = new WritePSM2Summary(files_list.getPsm_combined_threshold_csvfile(), files_list.getPsm_combined_threshold_summaryfile());
        mxo_percolator_psm.transcsv2summary();
        Integer mxo_percolator_pep_number = mxo_percolator_psm.getPeptide_number();
        Integer mxo_percolator_psm_number = mxo_percolator_psm.getPSM_number();

        FileWriter writeFile = new FileWriter(files_list.getStatistic_analysis_file());
        writeFile.write("the number of psm,peptide,protein\n");
        writeFile.write("\tPSM\tPeptide\tProtein\n");
        writeFile.write("iPeak\t" + mxo_percolator_psm_number + "\t" + mxo_percolator_pep_number + "\t" + mxo_Percolator_protein_number + "\n");

        if (Boolean.valueOf(debug)) {

            WriteProteins2Summary mxo_evalue_pro = new WriteProteins2Summary(files_list.getEprotein_combined_threshold_csvfile(), files_list.getEprotein_combined_threshold_summaryfile());
            mxo_evalue_pro.transcsv2summary();
            Integer mxo_Evalue_protein_number = mxo_evalue_pro.getProtein_number();
            WritePSM2Summary mxo_evalue_psm = new WritePSM2Summary(files_list.getEpsm_combined_threshold_csvfile(), files_list.getEpsm_combined_threshold_summaryfile());
            mxo_evalue_psm.transcsv2summary();
            Integer mxo_evalue_pep_number = mxo_evalue_psm.getPeptide_number();
            Integer mxo_evalue_psm_number = mxo_evalue_psm.getPSM_number();

            WriteProteins2Summary xtandem_percolator_pro = new WriteProteins2Summary(files_list.getXtandemPP_pfdr_threshold_group_procsv(), files_list.getXtandemPP_pfdr_threshold_group_prosummary());
            xtandem_percolator_pro.transcsv2summary();
            Integer xtandem_Percolator_protein_number = xtandem_percolator_pro.getProtein_number();
            WriteProteins2Summary xtandem_evalue_pro = new WriteProteins2Summary(files_list.getXtandemPP_efdr_threshold_group_procsv(), files_list.getXtandemPP_efdr_threshold_group_prosummary());
            xtandem_evalue_pro.transcsv2summary();
            Integer xtandem_Evalue_protein_number = xtandem_evalue_pro.getProtein_number();

            WriteProteins2Summary omssa_percolator_pro = new WriteProteins2Summary(files_list.getOmssaPP_pfdr_threshold_group_procsv(), files_list.getOmssaPP_pfdr_threshold_group_prosummary());
            omssa_percolator_pro.transcsv2summary();
            Integer omssa_Percolator_protein_number = omssa_percolator_pro.getProtein_number();
            WriteProteins2Summary omssa_evalue_pro = new WriteProteins2Summary(files_list.getOmssaPP_efdr_threshold_group_procsv(), files_list.getOmssaPP_efdr_threshold_group_prosummary());
            omssa_evalue_pro.transcsv2summary();
            Integer omssa_Evalue_protein_number = omssa_evalue_pro.getProtein_number();

            WriteProteins2Summary msgf_percolator_pro = new WriteProteins2Summary(files_list.getMsgfPP_pfdr_threshold_group_procsv(), files_list.getMsgfPP_pfdr_threshold_group_prosummary());
            msgf_percolator_pro.transcsv2summary();
            Integer msgf_Percolator_protein_number = msgf_percolator_pro.getProtein_number();
            WriteProteins2Summary msgf_evalue_pro = new WriteProteins2Summary(files_list.getMsgfPP_efdr_threshold_group_procsv(), files_list.getMsgfPP_efdr_threshold_group_prosummary());
            msgf_evalue_pro.transcsv2summary();
            Integer msgf_Evalue_protein_number = msgf_evalue_pro.getProtein_number();

            WritePSM2Summary xtandem_percolator_psm = new WritePSM2Summary(files_list.getXtandemPP_pfdr_threshold_group_psmcsv(), files_list.getXtandemPP_pfdr_threshold_group_psmsummary());
            xtandem_percolator_psm.transcsv2summary();
            Integer xtandem_percolator_pep_number = xtandem_percolator_psm.getPeptide_number();
            Integer xtandem_percolator_psm_number = xtandem_percolator_psm.getPSM_number();
            WritePSM2Summary xtandem_evalue_psm = new WritePSM2Summary(files_list.getXtandemPP_efdr_threshold_group_psmcsv(), files_list.getXtandemPP_efdr_threshold_group_psmsummary());
            xtandem_evalue_psm.transcsv2summary();
            Integer xtandem_evalue_pep_number = xtandem_evalue_psm.getPeptide_number();
            Integer xtandem_evalue_psm_number = xtandem_evalue_psm.getPSM_number();

            WritePSM2Summary omssa_percolator_psm = new WritePSM2Summary(files_list.getOmssaPP_pfdr_threshold_group_psmcsv(), files_list.getOmssaPP_pfdr_threshold_group_psmsummary());
            omssa_percolator_psm.transcsv2summary();
            Integer omssa_percolator_pep_number = omssa_percolator_psm.getPeptide_number();
            Integer omssa_percolator_psm_number = omssa_percolator_psm.getPSM_number();
            WritePSM2Summary omssa_evalue_psm = new WritePSM2Summary(files_list.getOmssaPP_efdr_threshold_group_psmcsv(), files_list.getOmssaPP_efdr_threshold_group_psmsummary());
            omssa_evalue_psm.transcsv2summary();
            Integer omssa_evalue_pep_number = omssa_evalue_psm.getPeptide_number();
            Integer omssa_evalue_psm_number = omssa_evalue_psm.getPSM_number();

            WritePSM2Summary msgf_percolator_psm = new WritePSM2Summary(files_list.getMsgfPP_pfdr_threshold_group_psmcsv(), files_list.getMsgfPP_pfdr_threshold_group_psmsummary());
            msgf_percolator_psm.transcsv2summary();
            Integer msgf_percolator_pep_number = msgf_percolator_psm.getPeptide_number();
            Integer msgf_percolator_psm_number = msgf_percolator_psm.getPSM_number();

            WritePSM2Summary msgf_evalue_psm = new WritePSM2Summary(files_list.getMsgfPP_efdr_threshold_group_psmcsv(), files_list.getMsgfPP_efdr_threshold_group_psmsummary());
            msgf_evalue_psm.transcsv2summary();
            Integer msgf_evalue_pep_number = msgf_evalue_psm.getPeptide_number();
            Integer msgf_evalue_psm_number = msgf_evalue_psm.getPSM_number();

            writeFile.write("FDRAnalyse\t" + mxo_evalue_psm_number + "\t" + mxo_evalue_pep_number + "\t" + mxo_Evalue_protein_number + "\n");
            writeFile.write("msgf\t" + msgf_evalue_psm_number + "\t" + msgf_evalue_pep_number + "\t" + msgf_Evalue_protein_number + "\n");
            writeFile.write("MP\t" + msgf_percolator_psm_number + "\t" + msgf_percolator_pep_number + "\t" + msgf_Percolator_protein_number + "\n");
            writeFile.write("omssa\t" + omssa_evalue_psm_number + "\t" + omssa_evalue_pep_number + "\t" + omssa_Evalue_protein_number + "\n");
            writeFile.write("OP\t" + omssa_percolator_psm_number + "\t" + omssa_percolator_pep_number + "\t" + omssa_Percolator_protein_number + "\n");
            writeFile.write("xtandem\t" + xtandem_evalue_psm_number + "\t" + xtandem_evalue_pep_number + "\t" + xtandem_Evalue_protein_number + "\n");
            writeFile.write("XP\t" + xtandem_percolator_psm_number + "\t" + xtandem_percolator_pep_number + "\t" + xtandem_Percolator_protein_number + "\n");

        }
        writeFile.close();
    }

    private Vector<String> get_featureFromSplitResult(String engine_name, File[] search_result_list) throws Exception {
        HashMap<String, Integer> enginename2engineID = new HashMap<String, Integer>();
        enginename2engineID.put("xtandem", 3);
        enginename2engineID.put("omssa", 2);
        enginename2engineID.put("msgf", 1);
        int engineID = enginename2engineID.get(engine_name);
        log.append("\n" + engine_name + "percolator...\n");
        log.setCaretPosition(log.getDocument().getLength());
        String out_dir = "";
        if (engineID == 3) {
            out_dir = files_list.getXtandem_out_dir();
        } else if (engineID == 2) {
            out_dir = files_list.getOmssa_out_dir();
        } else if (engineID == 1) {
            out_dir = files_list.getMsgf_out_dir();
        } else {
            log.append("Cannot get features! Caused by wrong output directory!\n");
            JOptionPane.showMessageDialog(null, "Cannot get features! Caused by wrong output directory!", "ERROR", JOptionPane.ERROR_MESSAGE);
            return null;
        }
        Vector<String> features_files_list = new Vector<String>();
        for (File result_file : search_result_list) {
            String primaryfilename = result_file.getName();
            primaryfilename = primaryfilename.substring(0, primaryfilename.lastIndexOf("."));
            primaryfilename = out_dir + primaryfilename;
            String result_file_path = result_file.getAbsolutePath();
            IPeakPercolator resul2feature;
            if (engineID == 2) {
                resul2feature = new IPeakPercolator(engineID, result_file_path, primaryfilename, database,  umod, decoyregrex);
            } else {
                resul2feature = new IPeakPercolator(engineID, result_file_path, primaryfilename, decoyregrex);
            }
            resul2feature.get_Features();
            features_files_list.add(resul2feature.getFeature_file());
        }
        return features_files_list;
    }

    private void combined_features(Vector<String> feature_file_list, String all_feature_file) throws IOException {

        File test_file = new File(all_feature_file);
        if (test_file.exists()) {
            test_file.delete();
        }
        BufferedReader readtitle = new BufferedReader(new FileReader(feature_file_list.get(0)));
        String title = readtitle.readLine();
        readtitle.close();
        Vector<String> featureinfors = new Vector<String>();
        for (String feature_file : feature_file_list) {
            log.append(feature_file + "\n");
            log.setCaretPosition(log.getDocument().getLength());
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

    private Set<String> split_percolatorResult(String pertxtfile) throws IOException {
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

    private Vector<String> addPercolator2mzid(Set<String> per_name, String mziddir, String out_dir) throws IOException {
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

    @Override
    protected Void doInBackground() throws Exception {
        //System.out.println("Begin...");
        log.setText("Begin at " + new Date() + "\n");
        if (!isCancelled()) {
            searchtype = set_filelist_and_get_analyse_method();
        } else {
            return null;
        }
        if (analyse_step != 0) {
            if (!isCancelled()) {
                single_step_run();
            } else {
                return null;
            }
        } else {
            if (!isCancelled()) {
                basic_analyse();
            } else {
                return null;
            }
        }
        return null;
    }

    @Override
    public void done() {
        try {
            if (!isCancelled()) {
                get();
                JOptionPane.showMessageDialog(null, "Done", "", JOptionPane.INFORMATION_MESSAGE);
            }
        } catch (InterruptedException e) {
            catchInformation(e);
        } catch (ExecutionException e) {
            catchInformation(e);
        }
        IPeakGui.getStartButton().setEnabled(true);
        IPeakGui.getStopButton().setEnabled(true);
    }

    private void catchInformation(Exception e) {
        String sep = "\r\n";
        String s = "";
        s += e.getMessage() + sep;
        StackTraceElement[] trace = e.getStackTrace();
        for (int i = 0; i < trace.length; i++) {
            if (i < 10) {
                s += "    at " + trace[i] + sep;
            } else {
                s += "    more " + (trace.length - 10) + " ..." + sep;
                break;
            }
        }
        Throwable ourCause = e.getCause();
        trace = ourCause.getStackTrace();
        s += "Caused by" + sep;
        for (int i = 0; i < trace.length; i++) {
            if (i < 10) {
                s += "    at " + trace[i] + sep;
            } else {
                s += "    more " + (trace.length - 10) + " ..." + sep;
                break;
            }
        }
        log.append(s);
        log.setCaretPosition(log.getDocument().getLength());
        JOptionPane.showMessageDialog(null, s, "", JOptionPane.ERROR_MESSAGE);
    }
}
