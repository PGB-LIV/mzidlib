package bgi.ipeak.util;

import java.io.File;
import java.util.Set;
import java.util.Vector;

public class IPeakFileList {
	//dir 
	private String combined_out_dir;
	private String omssa_out_dir;
	private String xtandem_out_dir;
	private String msgf_out_dir;
	
	//raw file name and output file regex;
	private String msgf_raw_file;
	private String xtandem_raw_file;
	private String omssa_raw_file;
	private String outfile_regex;
	
	//if the msgf_raw_file xtandem_raw_file omssa_raw_file is a dir,get the file list from the dir;
	//these 6 list used for search type 1;
	private Set<String> xtandem_per_list;
	private Set<String> omssa_per_list;
	private Set<String> msgf_per_list;
	
	private Vector<String> xtandem_addP_list;
	private Vector<String> omssa_addP_list;
	private Vector<String> msgf_addP_list;
	
	private String xtandem_primaryfilename;
	private String omssa_primaryfilename;
	private String msgf_primaryfilename;
	
	//after the trans step,the mzid file; used for search type 0;
	private String msgf_mzid;
	private String xtandem_mzid;
	private String omssa_mzid;
	
	//run percolaotor files 
	private String xtandem_feature_file;
	private String xtandem_percolator_result_file;
	private String xtandem_percolator_target_File;
	private String xtandem_percolator_decoy_file;
	
	private String omssa_feature_file;
	private String omssa_percolator_result_file;
	private String omssa_percolator_target_File;
	private String omssa_percolator_decoy_file;
	
	private String msgf_feature_file;
	private String msgf_percolator_result_file;
	private String msgf_percolator_target_File;
	private String msgf_percolator_decoy_file;
	
	/*the mzid file added percolator score;
		if search type ==0 ,this mzid file come from the msgf_mzid plus msgf_percolator result_file;xtandem,omssa,is the same;
		if search type ==1 ,this mzid file is combined from the addP_list;	
		xtandem_addP_list;
		omssa_addP_list;
		msgf_addP_list;
	 */
	private String xtandem_mzidAddPer;
	private String msgf_mzidAddPer;
	private String omssa_mzidAddPer;
	
	//the mzid file added percolator and proteins sequence file;
	private String xtandem_mzidAddProPer;
	private String msgf_mzidAddProPer;
	private String omssa_mzidAddProPer;
	
	//combined result file ,and fdr analyse result file
	private String combined_resultUsePerScore;		
	private String combined_resultUsePerScore_primaryfilename;
	private String combined_result_debugfile;	

	private String Ecombined_resultUseEvalue;	
	private String Ecombined_resultUseEvalue_primaryfilename;
	private String Ecombined_result_debugfile;	
	
	private String xtandemPP_efdr;
	private String omssaPP_efdr;
	private String msgfPP_efdr;
	private String xtandemPP_pfdr;
	private String omssaPP_pfdr;
	private String msgfPP_pfdr;
	
	//peptide leval fdr analyse result;
	private String pepFDR_Ecombined_resultUseEvalue;
	private String pepFDR_combined_resultUsePerScore;
	private String pepFDR_xtandemPP_efdr;
	private String pepFDR_omssaPP_efdr;
	private String pepFDR_msgfPP_efdr;
	private String pepFDR_xtandemPP_pfdr;
	private String pepFDR_omssaPP_pfdr;
	private String pepFDR_msgfPP_pfdr;
	
	private String statistic_analysis_file;
	
	// note used for ipeak protein group type;	
	private String combined_threshold_file;
	private String group_combined_threshold_file;
	private String psm_combined_threshold_csvfile;
	private String protein_combined_threshold_csvfile;
	private String psm_combined_threshold_summaryfile;
	private String protein_combined_threshold_summaryfile;
		
	private String Ecombined_threshold_file;
	private String Egroup_combined_threshold_file;
	private String Epsm_combined_threshold_csvfile;
	private String Eprotein_combined_threshold_csvfile;
	private String Epsm_combined_threshold_summaryfile;
	private String Eprotein_combined_threshold_summaryfile;
	
	private String xtandemPP_efdr_threshold;
	private String omssaPP_efdr_threshold;
	private String msgfPP_efdr_threshold;
	private String xtandemPP_pfdr_threshold;
	private String omssaPP_pfdr_threshold;
	private String msgfPP_pfdr_threshold;
	
	private String xtandemPP_efdr_threshold_group;
	private String omssaPP_efdr_threshold_group;
	private String msgfPP_efdr_threshold_group;
	private String xtandemPP_pfdr_threshold_group;
	private String omssaPP_pfdr_threshold_group;
	private String msgfPP_pfdr_threshold_group;
	
	private String xtandemPP_efdr_threshold_group_psmcsv;
	private String omssaPP_efdr_threshold_group_psmcsv;
	private String msgfPP_efdr_threshold_group_psmcsv;
	private String xtandemPP_pfdr_threshold_group_psmcsv;
	private String omssaPP_pfdr_threshold_group_psmcsv;
	private String msgfPP_pfdr_threshold_group_psmcsv;
	
	private String xtandemPP_efdr_threshold_group_procsv;
	private String omssaPP_efdr_threshold_group_procsv;
	private String msgfPP_efdr_threshold_group_procsv;
	private String xtandemPP_pfdr_threshold_group_procsv;
	private String omssaPP_pfdr_threshold_group_procsv;
	private String msgfPP_pfdr_threshold_group_procsv;
	
	private String xtandemPP_efdr_threshold_group_psmsummary;
	private String omssaPP_efdr_threshold_group_psmsummary;
	private String msgfPP_efdr_threshold_group_psmsummary;
	private String xtandemPP_pfdr_threshold_group_psmsummary;
	private String omssaPP_pfdr_threshold_group_psmsummary;
	private String msgfPP_pfdr_threshold_group_psmsummary;
	
	private String xtandemPP_efdr_threshold_group_prosummary;
	private String omssaPP_efdr_threshold_group_prosummary;
	private String msgfPP_efdr_threshold_group_prosummary;
	private String xtandemPP_pfdr_threshold_group_prosummary;
	private String omssaPP_pfdr_threshold_group_prosummary;
	private String msgfPP_pfdr_threshold_group_prosummary;
	// note used for ipeak protein group type;	
	
	public IPeakFileList(String msgf_path,String omssa_path,String xtandem_path,String out_dir,String regex_id) {
	
		msgf_raw_file=msgf_path;
		omssa_raw_file=omssa_path;
		xtandem_raw_file=xtandem_path;
		combined_out_dir=out_dir;	
		outfile_regex= regex_id;
		SetFilename();
		Set_MzidentmllibReportFile();
	}
	private void SetFilename() {
		if(combined_out_dir.lastIndexOf("/") !=(combined_out_dir.length()-1)){
			combined_out_dir+="/";
		}
		statistic_analysis_file=combined_out_dir+outfile_regex+"stat.txt";
    	xtandem_out_dir=combined_out_dir+"xtandem/";
    	omssa_out_dir=combined_out_dir+"omssa/";
    	msgf_out_dir=combined_out_dir+"msgf/";
    	
    	build_path();
    	
    	
    	msgf_primaryfilename=msgf_out_dir +outfile_regex+"_msgf";
    	xtandem_primaryfilename=xtandem_out_dir +outfile_regex+"_xtandem";
    	omssa_primaryfilename=omssa_out_dir +outfile_regex+"_omssa";
    	
    
    	String mzid_tag=".mzid";
    	msgf_mzid=msgf_primaryfilename+mzid_tag;
    	xtandem_mzid=xtandem_primaryfilename+mzid_tag;
    	omssa_mzid=omssa_primaryfilename+mzid_tag;
    	
    	
    	String addPmzid_tat="AddP.mzid";
    	xtandem_mzidAddPer=xtandem_primaryfilename+addPmzid_tat;
    	msgf_mzidAddPer=msgf_primaryfilename+addPmzid_tat;
    	omssa_mzidAddPer=omssa_primaryfilename+addPmzid_tat;
    	
    	
    	String addPPmzid_tag="AddPP.mzid";
    	xtandem_mzidAddProPer=xtandem_primaryfilename+addPPmzid_tag;
    	msgf_mzidAddProPer=msgf_primaryfilename+addPPmzid_tag;
    	omssa_mzidAddProPer=omssa_primaryfilename+addPPmzid_tag;

    	String featrue_tag=".features";
    	String percolator_result_tag=".per.txt";
    	String percolator_result_target_tag=".target.txt";
    	String percolator_result_decoy_tag=".decoy.txt";
    	xtandem_feature_file=xtandem_primaryfilename+featrue_tag;
    	xtandem_percolator_result_file=xtandem_primaryfilename+percolator_result_tag;
    	xtandem_percolator_target_File=xtandem_primaryfilename+percolator_result_target_tag;
    	xtandem_percolator_decoy_file=xtandem_primaryfilename+percolator_result_decoy_tag;

    	omssa_feature_file=omssa_primaryfilename+featrue_tag;
    	omssa_percolator_result_file=omssa_primaryfilename+percolator_result_tag;
    	omssa_percolator_target_File=omssa_primaryfilename+percolator_result_target_tag;
    	omssa_percolator_decoy_file=omssa_primaryfilename+percolator_result_decoy_tag;

    	msgf_feature_file=msgf_primaryfilename+featrue_tag;
    	msgf_percolator_result_file=msgf_primaryfilename+percolator_result_tag;
    	msgf_percolator_target_File=msgf_primaryfilename+percolator_result_target_tag;
    	msgf_percolator_decoy_file=msgf_primaryfilename+percolator_result_decoy_tag;
    	
    	String evalue_fdr_tag="_efdr";
    	String percolatorscore_fdr_tag="_pfdr";
    	xtandemPP_efdr=xtandem_primaryfilename+evalue_fdr_tag+mzid_tag;
    	omssaPP_efdr=omssa_primaryfilename+evalue_fdr_tag+mzid_tag;
    	msgfPP_efdr=msgf_primaryfilename+evalue_fdr_tag+mzid_tag;
    	xtandemPP_pfdr=xtandem_primaryfilename+percolatorscore_fdr_tag+mzid_tag;
    	omssaPP_pfdr=omssa_primaryfilename+percolatorscore_fdr_tag+mzid_tag;
    	msgfPP_pfdr=msgf_primaryfilename+percolatorscore_fdr_tag+mzid_tag;
    	
    	String combined_tag = "_combined";
    	String debug_tag ="debug.csv";    	
    	combined_resultUsePerScore=combined_out_dir+outfile_regex+combined_tag+percolatorscore_fdr_tag+mzid_tag;    	
    	combined_resultUsePerScore_primaryfilename=combined_out_dir+outfile_regex+combined_tag+percolatorscore_fdr_tag;
    	combined_result_debugfile=combined_out_dir+outfile_regex+percolatorscore_fdr_tag+debug_tag;

    	Ecombined_resultUseEvalue=combined_out_dir+outfile_regex+combined_tag+evalue_fdr_tag+mzid_tag;	
      	Ecombined_resultUseEvalue_primaryfilename=combined_out_dir+outfile_regex+combined_tag+evalue_fdr_tag;  	
    	Ecombined_result_debugfile=combined_out_dir+outfile_regex+evalue_fdr_tag+debug_tag; 

    	//peptdie leval fdr resutl
    	String peptide_leval_tag= "_pepFDR";
    	pepFDR_xtandemPP_efdr = xtandem_primaryfilename+peptide_leval_tag+evalue_fdr_tag+mzid_tag;
    	pepFDR_omssaPP_efdr = omssa_primaryfilename+peptide_leval_tag+evalue_fdr_tag+mzid_tag;
    	pepFDR_msgfPP_efdr = msgf_primaryfilename+peptide_leval_tag+evalue_fdr_tag+mzid_tag;
    	pepFDR_xtandemPP_pfdr = xtandem_primaryfilename+peptide_leval_tag+percolatorscore_fdr_tag+mzid_tag;
    	pepFDR_omssaPP_pfdr = omssa_primaryfilename+peptide_leval_tag+percolatorscore_fdr_tag+mzid_tag;
    	pepFDR_msgfPP_pfdr = msgf_primaryfilename+peptide_leval_tag+percolatorscore_fdr_tag+mzid_tag;
    	   	
    	pepFDR_combined_resultUsePerScore = combined_out_dir+outfile_regex+peptide_leval_tag+combined_tag+percolatorscore_fdr_tag+mzid_tag;
    	pepFDR_Ecombined_resultUseEvalue = combined_out_dir+outfile_regex+peptide_leval_tag+combined_tag+evalue_fdr_tag+mzid_tag;	

    	
	}
	private void Set_MzidentmllibReportFile(){		
    	Ecombined_threshold_file=combined_out_dir+outfile_regex +"combined_evalue_threshold.mzid";
    	Egroup_combined_threshold_file=combined_out_dir+outfile_regex+"combined_evalue_threshold_group.mzid";
    	Epsm_combined_threshold_csvfile=combined_out_dir+outfile_regex+"combined_evalue_psm.csv";
    	Eprotein_combined_threshold_csvfile=combined_out_dir+outfile_regex+"combined_evalue_pro.csv";
    	Epsm_combined_threshold_summaryfile=combined_out_dir+outfile_regex+"combined_evalue_psmsummary.txt";
    	Eprotein_combined_threshold_summaryfile=combined_out_dir+outfile_regex+"combined_evalue_prosummary.txt";
		
    	combined_threshold_file=combined_out_dir+outfile_regex +"combined_percolator_threshold.mzid";
    	group_combined_threshold_file=combined_out_dir+outfile_regex+"combined_percolator_threshold_group.mzid";
    	psm_combined_threshold_csvfile=combined_out_dir+outfile_regex+"combined_percolator_psm.csv";
    	protein_combined_threshold_csvfile=combined_out_dir+outfile_regex+"combined_percolator_pro.csv";
    	psm_combined_threshold_summaryfile=combined_out_dir+outfile_regex+"combined_percolator_psmsummary.txt";
    	protein_combined_threshold_summaryfile=combined_out_dir+outfile_regex+"combined_percolator_prosummary.txt";
		
		xtandemPP_efdr_threshold=xtandem_primaryfilename+"efdrT.mzid";
    	omssaPP_efdr_threshold=omssa_primaryfilename+"efdrT.mzid";
    	msgfPP_efdr_threshold=msgf_primaryfilename+"efdrT.mzid";
    	xtandemPP_pfdr_threshold=xtandem_primaryfilename+"pfdrT.mzid";
    	omssaPP_pfdr_threshold=omssa_primaryfilename+"pfdrT.mzid";
    	msgfPP_pfdr_threshold=msgf_primaryfilename+"pfdrT.mzid";
    	
    	xtandemPP_efdr_threshold_group=xtandem_primaryfilename+"efdrTG.mzid";
    	omssaPP_efdr_threshold_group=omssa_primaryfilename+"efdrTG.mzid";
    	msgfPP_efdr_threshold_group=msgf_primaryfilename+"efdrTG.mzid";
    	xtandemPP_pfdr_threshold_group=xtandem_primaryfilename+"pfdrTG.mzid";
    	omssaPP_pfdr_threshold_group=omssa_primaryfilename+"pfdrTG.mzid";
    	msgfPP_pfdr_threshold_group=msgf_primaryfilename+"pfdrTG.mzid";
    	
    	xtandemPP_efdr_threshold_group_psmcsv=xtandem_primaryfilename+"efdr_psm.csv";
    	omssaPP_efdr_threshold_group_psmcsv=omssa_primaryfilename+"efdr_psm.csv";
    	msgfPP_efdr_threshold_group_psmcsv=msgf_primaryfilename+"efdr_psm.csv";
    	xtandemPP_pfdr_threshold_group_psmcsv=xtandem_primaryfilename+"pfdr_psm.csv";
    	omssaPP_pfdr_threshold_group_psmcsv=omssa_primaryfilename+"pfdr_psm.csv";
    	msgfPP_pfdr_threshold_group_psmcsv=msgf_primaryfilename+"pfdr_psm.csv";
    	
    	xtandemPP_efdr_threshold_group_procsv=xtandem_primaryfilename+"efdr_pro.csv";
    	omssaPP_efdr_threshold_group_procsv=omssa_primaryfilename+"efdr_pro.csv";
    	msgfPP_efdr_threshold_group_procsv=msgf_primaryfilename+"efdr_pro.csv";
    	xtandemPP_pfdr_threshold_group_procsv=xtandem_primaryfilename+"pfdr_pro.csv";
    	omssaPP_pfdr_threshold_group_procsv=omssa_primaryfilename+"pfdr_pro.csv";
    	msgfPP_pfdr_threshold_group_procsv=msgf_primaryfilename+"pfdr_pro.csv";
    	
       	xtandemPP_efdr_threshold_group_psmsummary=xtandem_primaryfilename+"efdr_psmsummary.txt";
    	omssaPP_efdr_threshold_group_psmsummary=omssa_primaryfilename+"efdr_psmsummary.txt";
    	msgfPP_efdr_threshold_group_psmsummary=msgf_primaryfilename+"efdr_psmsummary.txt";
    	xtandemPP_pfdr_threshold_group_psmsummary=xtandem_primaryfilename+"pfdr_psmsummary.txt";
    	omssaPP_pfdr_threshold_group_psmsummary=omssa_primaryfilename+"pfdr_psmsummary.txt";
    	msgfPP_pfdr_threshold_group_psmsummary=msgf_primaryfilename+"pfdr_psmsummary.txt";
    	
    	xtandemPP_efdr_threshold_group_prosummary=xtandem_primaryfilename+"efdr_prosummary.txt";
    	omssaPP_efdr_threshold_group_prosummary=omssa_primaryfilename+"efdr_prosummary.txt";
    	msgfPP_efdr_threshold_group_prosummary=msgf_primaryfilename+"efdr_prosummary.txt";
    	xtandemPP_pfdr_threshold_group_prosummary=xtandem_primaryfilename+"pfdr_prosummary.txt";
    	omssaPP_pfdr_threshold_group_prosummary=omssa_primaryfilename+"pfdr_prosummary.txt";
    	msgfPP_pfdr_threshold_group_prosummary=msgf_primaryfilename+"pfdr_prosummary.txt";
	}
	
	private void build_path() {
		File combined_file = new File(combined_out_dir);
		File xtandem_file=new File(xtandem_out_dir);
		File omssa_file=new File(omssa_out_dir);
		File msgf_file=new File(msgf_out_dir);
		if(combined_file!= null && !combined_file.exists()){
			combined_file.mkdir();
			xtandem_file.mkdir();
			omssa_file.mkdir();
			msgf_file.mkdir();			
		}
		else{
			if(xtandem_file!= null && !xtandem_file.exists()){
				xtandem_file.mkdir();
			}
			if(omssa_file!= null && !omssa_file.exists()){
				omssa_file.mkdir();
			}
			if(msgf_file!= null && !msgf_file.exists()){
				msgf_file.mkdir();	
			}
		}
	}
	
	
	public String getCombined_out_dir() {
		return combined_out_dir;
	}
	public String getCombined_resultUsePerScore_primaryfilename() {
		return combined_resultUsePerScore_primaryfilename;
	}

	public String getCombined_resultUsePerScore() {
		return combined_resultUsePerScore;
	}
	public String getMsgf_feature_file() {
		return msgf_feature_file;
	}
	public String getMsgf_mzid() {
		return msgf_mzid;
	}
	public String getMsgf_mzidAddPer() {
		return msgf_mzidAddPer;
	}
	public String getMsgf_mzidAddProPer() {
		return msgf_mzidAddProPer;
	}
	public String getMsgf_out_dir() {
		return msgf_out_dir;
	}
	public String getMsgf_percolator_decoy_file() {
		return msgf_percolator_decoy_file;
	}
	public String getMsgf_percolator_result_file() {
		return msgf_percolator_result_file;
	}
	public String getMsgf_percolator_target_File() {
		return msgf_percolator_target_File;
	}
	public String getMsgf_raw_file() {
		return msgf_raw_file;
	}
	public String getOmssa_feature_file() {
		return omssa_feature_file;
	}
	public String getOmssa_mzid() {
		return omssa_mzid;
	}
	public String getOmssa_mzidAddPer() {
		return omssa_mzidAddPer;
	}
	public String getOmssa_mzidAddProPer() {
		return omssa_mzidAddProPer;
	}
	public String getOmssa_out_dir() {
		return omssa_out_dir;
	}
	public String getOmssa_percolator_decoy_file() {
		return omssa_percolator_decoy_file;
	}
	public String getOmssa_percolator_result_file() {
		return omssa_percolator_result_file;
	}
	public String getOmssa_percolator_target_File() {
		return omssa_percolator_target_File;
	}
	public String getOmssa_raw_file() {
		return omssa_raw_file;
	}
	public String getoutfile_regex() {
		return outfile_regex;
	}
	public String getXtandem_feature_file() {
		return xtandem_feature_file;
	}
	public String getXtandem_mzid() {
		return xtandem_mzid;
	}
	public String getXtandem_mzidAddPer() {
		return xtandem_mzidAddPer;
	}
	public String getXtandem_mzidAddProPer() {
		return xtandem_mzidAddProPer;
	}
	public String getXtandem_out_dir() {
		return xtandem_out_dir;
	}
	public String getXtandem_percolator_decoy_file() {
		return xtandem_percolator_decoy_file;
	}
	public String getXtandem_percolator_result_file() {
		return xtandem_percolator_result_file;
	}
	public String getXtandem_percolator_target_File() {
		return xtandem_percolator_target_File;
	}
	public String getXtandem_raw_file() {
		return xtandem_raw_file;
	}
	public String getCombined_result_debugfile() {
		return combined_result_debugfile;
	}
	public String getMsgf_primaryfilename() {
		return msgf_primaryfilename;
	}
	public String getOmssa_primaryfilename() {
		return omssa_primaryfilename;
	}
	public String getXtandem_primaryfilename() {
		return xtandem_primaryfilename;
	}
	public String getCombined_threshold_file() {
		return combined_threshold_file;
	}
	public String getGroup_combined_threshold_file() {
		return group_combined_threshold_file;
	}
	public String getProtein_combined_threshold_csvfile() {
		return protein_combined_threshold_csvfile;
	}
	public String getPsm_combined_threshold_csvfile() {
		return psm_combined_threshold_csvfile;
	}
	public String getEcombined_result_debugfile() {
		return Ecombined_result_debugfile;
	}
	public String getEcombined_resultUseEvalue() {
		return Ecombined_resultUseEvalue;
	}
	public String getEcombined_resultUseEvalue_primaryfilename() {
		return Ecombined_resultUseEvalue_primaryfilename;
	}
	public String getEcombined_threshold_file() {
		return Ecombined_threshold_file;
	}
	public String getEgroup_combined_threshold_file() {
		return Egroup_combined_threshold_file;
	}
	public String getEprotein_combined_threshold_csvfile() {
		return Eprotein_combined_threshold_csvfile;
	}
	public String getEpsm_combined_threshold_csvfile() {
		return Epsm_combined_threshold_csvfile;
	}
	public String getMsgfPP_efdr() {
		return msgfPP_efdr;
	}
	public String getMsgfPP_efdr_threshold() {
		return msgfPP_efdr_threshold;
	}
	public String getMsgfPP_efdr_threshold_group() {
		return msgfPP_efdr_threshold_group;
	}
	public String getMsgfPP_efdr_threshold_group_procsv() {
		return msgfPP_efdr_threshold_group_procsv;
	}
	public String getMsgfPP_efdr_threshold_group_psmcsv() {
		return msgfPP_efdr_threshold_group_psmcsv;
	}
	public String getMsgfPP_pfdr() {
		return msgfPP_pfdr;
	}
	public String getMsgfPP_pfdr_threshold() {
		return msgfPP_pfdr_threshold;
	}
	public String getMsgfPP_pfdr_threshold_group() {
		return msgfPP_pfdr_threshold_group;
	}
	public String getMsgfPP_pfdr_threshold_group_procsv() {
		return msgfPP_pfdr_threshold_group_procsv;
	}
	public String getMsgfPP_pfdr_threshold_group_psmcsv() {
		return msgfPP_pfdr_threshold_group_psmcsv;
	}
	public String getOmssaPP_efdr() {
		return omssaPP_efdr;
	}
	public String getOmssaPP_efdr_threshold() {
		return omssaPP_efdr_threshold;
	}
	public String getOmssaPP_efdr_threshold_group() {
		return omssaPP_efdr_threshold_group;
	}
	public String getOmssaPP_efdr_threshold_group_procsv() {
		return omssaPP_efdr_threshold_group_procsv;
	}
	public String getOmssaPP_efdr_threshold_group_psmcsv() {
		return omssaPP_efdr_threshold_group_psmcsv;
	}
	public String getOmssaPP_pfdr() {
		return omssaPP_pfdr;
	}
	public String getOmssaPP_pfdr_threshold() {
		return omssaPP_pfdr_threshold;
	}
	public String getOmssaPP_pfdr_threshold_group() {
		return omssaPP_pfdr_threshold_group;
	}
	public String getOmssaPP_pfdr_threshold_group_procsv() {
		return omssaPP_pfdr_threshold_group_procsv;
	}
	public String getOmssaPP_pfdr_threshold_group_psmcsv() {
		return omssaPP_pfdr_threshold_group_psmcsv;
	}
	public String getXtandemPP_efdr() {
		return xtandemPP_efdr;
	}
	public String getXtandemPP_efdr_threshold() {
		return xtandemPP_efdr_threshold;
	}
	public String getXtandemPP_efdr_threshold_group() {
		return xtandemPP_efdr_threshold_group;
	}
	public String getXtandemPP_efdr_threshold_group_procsv() {
		return xtandemPP_efdr_threshold_group_procsv;
	}
	public String getXtandemPP_efdr_threshold_group_psmcsv() {
		return xtandemPP_efdr_threshold_group_psmcsv;
	}
	public String getXtandemPP_pfdr() {
		return xtandemPP_pfdr;
	}
	public String getXtandemPP_pfdr_threshold() {
		return xtandemPP_pfdr_threshold;
	}
	public String getXtandemPP_pfdr_threshold_group() {
		return xtandemPP_pfdr_threshold_group;
	}
	public String getXtandemPP_pfdr_threshold_group_procsv() {
		return xtandemPP_pfdr_threshold_group_procsv;
	}
	public String getXtandemPP_pfdr_threshold_group_psmcsv() {
		return xtandemPP_pfdr_threshold_group_psmcsv;
	}
	public String getEprotein_combined_threshold_summaryfile() {
		return Eprotein_combined_threshold_summaryfile;
	}
	public String getEpsm_combined_threshold_summaryfile() {
		return Epsm_combined_threshold_summaryfile;
	}
	public String getMsgfPP_efdr_threshold_group_prosummary() {
		return msgfPP_efdr_threshold_group_prosummary;
	}
	public String getMsgfPP_efdr_threshold_group_psmsummary() {
		return msgfPP_efdr_threshold_group_psmsummary;
	}
	public String getMsgfPP_pfdr_threshold_group_prosummary() {
		return msgfPP_pfdr_threshold_group_prosummary;
	}
	public String getMsgfPP_pfdr_threshold_group_psmsummary() {
		return msgfPP_pfdr_threshold_group_psmsummary;
	}
	public String getOmssaPP_efdr_threshold_group_prosummary() {
		return omssaPP_efdr_threshold_group_prosummary;
	}
	public String getOmssaPP_efdr_threshold_group_psmsummary() {
		return omssaPP_efdr_threshold_group_psmsummary;
	}
	public String getOmssaPP_pfdr_threshold_group_prosummary() {
		return omssaPP_pfdr_threshold_group_prosummary;
	}
	public String getOmssaPP_pfdr_threshold_group_psmsummary() {
		return omssaPP_pfdr_threshold_group_psmsummary;
	}
	public String getProtein_combined_threshold_summaryfile() {
		return protein_combined_threshold_summaryfile;
	}
	public String getPsm_combined_threshold_summaryfile() {
		return psm_combined_threshold_summaryfile;
	}
	public String getXtandemPP_efdr_threshold_group_prosummary() {
		return xtandemPP_efdr_threshold_group_prosummary;
	}
	public String getXtandemPP_efdr_threshold_group_psmsummary() {
		return xtandemPP_efdr_threshold_group_psmsummary;
	}
	public String getXtandemPP_pfdr_threshold_group_prosummary() {
		return xtandemPP_pfdr_threshold_group_prosummary;
	}
	public String getXtandemPP_pfdr_threshold_group_psmsummary() {
		return xtandemPP_pfdr_threshold_group_psmsummary;
	}
	public String getStatistic_analysis_file() {
		return statistic_analysis_file;
	}
	public Set<String> getMsgf_per_list() {
		return msgf_per_list;
	}
	public Set<String> getOmssa_per_list() {
		return omssa_per_list;
	}
	public Set<String> getXtandem_per_list() {
		return xtandem_per_list;
	}	
	public Vector<String> getMsgf_addP_list() {
		return msgf_addP_list;
	}
	public Vector<String> getOmssa_addP_list() {
		return omssa_addP_list;
	}
	public Vector<String> getXtandem_addP_list() {
		return xtandem_addP_list;
	}
	public String getOutfile_regex() {
		return outfile_regex;
	}
	public String getPepFDR_combined_resultUsePerScore() {
		return pepFDR_combined_resultUsePerScore;
	}
	public String getPepFDR_Ecombined_resultUseEvalue() {
		return pepFDR_Ecombined_resultUseEvalue;
	}
	public String getPepFDR_msgfPP_efdr() {
		return pepFDR_msgfPP_efdr;
	}
	public String getPepFDR_msgfPP_pfdr() {
		return pepFDR_msgfPP_pfdr;
	}
	public String getPepFDR_omssaPP_efdr() {
		return pepFDR_omssaPP_efdr;
	}
	public String getPepFDR_omssaPP_pfdr() {
		return pepFDR_omssaPP_pfdr;
	}
	public String getPepFDR_xtandemPP_efdr() {
		return pepFDR_xtandemPP_efdr;
	}
	public String getPepFDR_xtandemPP_pfdr() {
		return pepFDR_xtandemPP_pfdr;
	}	
	public void setXtandem_per_list(Set<String> xtandem_per_list) {
		this.xtandem_per_list = xtandem_per_list;
	}
	public void setOmssa_per_list(Set<String> omssa_per_list) {
		this.omssa_per_list = omssa_per_list;
	}
	public void setMsgf_per_list(Set<String> msgf_per_list) {
		this.msgf_per_list = msgf_per_list;
	}
	public void setOmssa_addP_list(Vector<String> omssa_addP_list) {
		this.omssa_addP_list = omssa_addP_list;
	}
	public void setXtandem_addP_list(Vector<String> xtandem_addP_list) {
		this.xtandem_addP_list = xtandem_addP_list;
	}
	public void setMsgf_addP_list(Vector<String> msgf_addP_list) {
		this.msgf_addP_list = msgf_addP_list;
	}
}
