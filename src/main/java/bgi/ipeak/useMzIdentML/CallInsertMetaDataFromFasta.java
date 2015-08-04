package bgi.ipeak.useMzIdentML;

import uk.ac.liv.mzidlib.fasta.InsertMetaDataFromFasta;

public class CallInsertMetaDataFromFasta {
	private String mzid2add;
	private String output_mzid;
	private String fasta_file;
	private String accessionSplitRegex;
	
	public static void main(String[] args) {
		

	}
	public CallInsertMetaDataFromFasta(String mzid2add,String output_mzid,String fasta_file,String accessionSplitregex) {
		
		this.mzid2add=mzid2add;
		this.output_mzid=output_mzid;
		this.fasta_file=fasta_file;
		this.accessionSplitRegex=accessionSplitregex.replaceAll("/", "");
	}
	@SuppressWarnings("unused")
	public void AddFasta_UseMzidlib(){
		InsertMetaDataFromFasta insertMD = new InsertMetaDataFromFasta(mzid2add, output_mzid, fasta_file, accessionSplitRegex);
	}

}
