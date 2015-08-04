package bgi.ipeak;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.MzIdentML;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import bgi.ipeak.io.PSMSummaryInformation;
import bgi.ipeak.io.ProteinSummaryInformation;

public class Mzid2Summary {
	 @Option(name="-f",required=true,usage="(required),input mzid file")
	private String mzid_file;
	 @Option(name="-o",required=true,usage="(required),output file name")
	private String output_file;
	 @Option(name="-d",required=true,usage="(required),fasta format database")
	private String database_file;
	 @Option(name="-t",required=false,usage="(options),threshold score value,default(0.01)")
	private Double threshold_score =0.01;
	 @Option(name="-tN",required=true,usage="(required),threshold score cvParm")
	private String main_score_cv="";
	 @Option(name="-betterScoreAreLower",required=false,usage="(options),batter score are lower,default(true)")
	private Boolean betterScoreAreLower=true;
	 @Option(name="-rank",required=false,usage="(options),use only rank 1 psm,default (true)")
	private Boolean use_only_rank1=true;
	 @Option(name="-decoyRegex",required=false,usage="(options),decoy tag,default(###REV###)")
	private String decoy_tag="###REV###";
	
	
	private Vector<PSMSummaryInformation> psm_infomation=new Vector<PSMSummaryInformation>();
	private Vector<ProteinSummaryInformation> pro_infomation=new Vector<ProteinSummaryInformation>();
	private HashMap<String, String> spectrum_file_list=new HashMap<String, String>();
	private HashMap<String,String> Database_path =new HashMap<String, String>();
	private Integer protein_number=0;
	private Integer peptide_number=0;
	private Integer PSM_number=0;
	
	public static void main(String[] args) throws IOException {

    	
    	Mzid2Summary mzid2summary=new Mzid2Summary();
		CmdLineParser parser = new CmdLineParser(mzid2summary);
		try {
	    	parser.setUsageWidth(100);
	    	parser.parseArgument(args);
	    	System.err.println("\niPeak v1.0(2013-11)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n" +
        				"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    } catch (CmdLineException e) {
	    	System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n" +
						"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    	System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.Mzid2Summary  [options...]\n");
	    	parser.printUsage(System.err);
	    	System.err.println("\nError:");
	    	System.err.println(e.getMessage());
	    	System.err.println("");
	    	return;
	    }    
		mzid2summary.read_inmzid();
		mzid2summary.get_protein_information();
		mzid2summary.write2txt();
	}
	
	public Mzid2Summary(String mzid_file,String output_file,String database_file,Double threshold_score,
			Boolean betterScoreAreLower,String cvParam,String DecoyRegex,Boolean use_only_rank1) {

		this.mzid_file=mzid_file;
		this.output_file=output_file;
		this.database_file=database_file;
		this.threshold_score=threshold_score;
		this.betterScoreAreLower=betterScoreAreLower;
		this.main_score_cv=cvParam;
		this.decoy_tag=DecoyRegex;
		this.use_only_rank1=use_only_rank1;
	}
	public Mzid2Summary() {

		
	}
	public void trans_mzid2txt() throws IOException {
		psm_infomation=new Vector<PSMSummaryInformation>();
		pro_infomation=new Vector<ProteinSummaryInformation>();
		spectrum_file_list=new HashMap<String, String>();
		Database_path =new HashMap<String, String>();
		protein_number=0;
		peptide_number=0;
		PSM_number=0;		
		read_inmzid();
		get_protein_information();
		tatistics_numbers();
		write2txt();		
	}
	private void read_inmzid() {

		MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid_file));
		MzIdentML mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
		

		SequenceCollection aCollection =mzIdentML.getSequenceCollection();
		HashMap<String, DBSequence> id_dbseq = new HashMap<String, DBSequence>();
		for (DBSequence db : aCollection.getDBSequence()) {
			id_dbseq.put(db.getId(), db);
		}
		HashMap<String, Peptide> id_peptide = new HashMap<String, Peptide>();
		for (Peptide pep : aCollection.getPeptide()) {
			id_peptide.put(pep.getId(), pep);
		}
		HashMap<String, PeptideEvidence> id_pepE = new HashMap<String, PeptideEvidence>();
		for (PeptideEvidence db : aCollection.getPeptideEvidence()) {
			id_pepE.put(db.getId(), db);
		}
		

		List<SpectraData> spectrum_data= mzIdentML.getDataCollection().getInputs().getSpectraData();
		for (SpectraData spectraData : spectrum_data) {
			spectrum_file_list.put(spectraData.getId(), spectraData.getLocation());
		}
		List<SearchDatabase> databases= mzIdentML.getDataCollection().getInputs().getSearchDatabase();
		for (SearchDatabase searchDatabase : databases) {
			Database_path.put(searchDatabase.getId(), searchDatabase.getLocation());
		}		 
		List<SpectrumIdentificationList> SpectrumIdentificationList =
				mzIdentML.getDataCollection().getAnalysisData().getSpectrumIdentificationList();
		for (SpectrumIdentificationList spectrumIdentificationList2 : SpectrumIdentificationList) {
			List<SpectrumIdentificationResult> spectrum_result= spectrumIdentificationList2.getSpectrumIdentificationResult();

			for (SpectrumIdentificationResult spectrumIdentificationResult : spectrum_result) {				
				List<SpectrumIdentificationItem> psm=spectrumIdentificationResult.getSpectrumIdentificationItem();

				for (SpectrumIdentificationItem spectrumIdentificationItem : psm) {

					if(use_only_rank1){
						if(spectrumIdentificationItem.getRank()>1){
							continue;
						}
					}
					PSMSummaryInformation psm_infor = new PSMSummaryInformation();
					psm_infor.setSpectrum_index(spectrumIdentificationResult.getSpectrumID());
					psm_infor.setRank(spectrumIdentificationItem.getRank());
					
					String peptide_seq = id_peptide.get(spectrumIdentificationItem.getPeptideRef()).getPeptideSequence();
					Vector<String> mods=new Vector<String>();
					for (uk.ac.ebi.jmzidml.model.mzidml.Modification modification : id_peptide.get(spectrumIdentificationItem.getPeptideRef()).getModification()) {
						Integer location=modification.getLocation();
						mods.add(modification.getMonoisotopicMassDelta()+"@"+location);
					}
					psm_infor.setMods(mods);
					psm_infor.setPeptide_sequence(peptide_seq);
					psm_infor.setCharge(spectrumIdentificationItem.getChargeState());
					psm_infor.setMz(spectrumIdentificationItem.getExperimentalMassToCharge());
					psm_infor.setTheory_mz(spectrumIdentificationItem.getCalculatedMassToCharge());
					psm_infor.setMz_error(spectrumIdentificationItem.getExperimentalMassToCharge()-spectrumIdentificationItem.getCalculatedMassToCharge());
					psm_infor.setPsm_id(spectrumIdentificationItem.getId());
					List<PeptideEvidenceRef> peptide_evidence = spectrumIdentificationItem.getPeptideEvidenceRef();
					
	
					Vector<String> proteins = new Vector<String>();					
					Boolean is_decoy=false;
					for (PeptideEvidenceRef peptideEvidenceRef : peptide_evidence) {
						String protein_ref = id_pepE.get(peptideEvidenceRef.getPeptideEvidenceRef()).getDBSequenceRef();
						String protein_acc=id_dbseq.get(protein_ref).getAccession();
						proteins.add(protein_acc);
						if(protein_acc.contains(decoy_tag)){
							is_decoy=true;
						}
					}
					psm_infor.setProteins(proteins);
					psm_infor.setIs_decoy(is_decoy);
					

					List<CvParam> cvParams=spectrumIdentificationItem.getCvParam();
					HashMap<String, String> scores=new HashMap<String, String>();
					for (CvParam cvParam : cvParams) {
						if(cvParam.getAccession().equals(main_score_cv)){
							HashMap<String, String> main_score = new HashMap<String, String>();
							main_score.put(cvParam.getName(), cvParam.getValue());
							psm_infor.setMain_threshold_score(main_score);
							if(Double.valueOf(cvParam.getValue())<threshold_score){
								psm_infor.setPass_threshold(betterScoreAreLower);
							}
							else{
								psm_infor.setPass_threshold(!betterScoreAreLower);
							}
						}
						else{
							scores.put(cvParam.getName(),cvParam.getValue());
						}
						
					}
					psm_infor.setScores(scores);
					psm_infomation.add(psm_infor);
				}				
			}
		}		
	}


	private void get_protein_information() throws IOException{
		ProteinSummary pro_summary = new ProteinSummary(psm_infomation,database_file);
		pro_infomation= pro_summary.get_proteinSummary();
	}
	
	private void write2txt() throws IOException {
		BufferedWriter writeFile=new BufferedWriter(new FileWriter( output_file+"PSMSummary.txt"));
		psm_infomation.get(0).print_tile(writeFile);
		for (PSMSummaryInformation psm : psm_infomation) {
			psm.print_psm_infor(writeFile);
		}
		writeFile.close();
		
		BufferedWriter writePro=new BufferedWriter(new FileWriter(output_file+"ProteinSummary.txt"));
		writePro.write("GroupID\tAccession\tMass\tPeptideNum\tUniqePepNum\tSpectrumNum\tUniqSpectra\tUniquePep\tRazorPep\tSameSet\tDescription\n");
		writePro.flush();
		for (ProteinSummaryInformation pro : pro_infomation) {
			pro.print_protein(writePro);
		}		
		writePro.close();
	}

	private void tatistics_numbers() {
		PSM_number=0;
		Vector<String> peptides=new Vector<String>();
		for (PSMSummaryInformation psm : psm_infomation) {
			if(!psm.getIs_decoy() && psm.getPass_threshold()){
				PSM_number++;
				if(!peptides.contains(psm.getPeptide_sequence())){
					peptides.add(psm.getPeptide_sequence());
				}
			}
		}
		peptide_number=peptides.size();
		protein_number=pro_infomation.size();
	}
	public Integer getPSM_number() {
		return PSM_number;
	}
	public Integer getProtein_number() {
		return protein_number;
	}
	public Integer getPeptide_number() {
		return peptide_number;
	}
}
