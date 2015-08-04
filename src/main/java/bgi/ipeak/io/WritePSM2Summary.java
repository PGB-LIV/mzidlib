package bgi.ipeak.io;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import au.com.bytecode.opencsv.CSVReader;

public class WritePSM2Summary {
	@Option(name="-f",required=true,usage="(required) The input psm csv file")
	private String psm_csv;
	@Option(name="-o",required=true,usage="(required) The output psm csv file")
	private String output_psm;
	private Integer peptide_number;
	private Integer PSM_number;

	public static void main(String[] args) throws IOException {
		WritePSM2Summary writepsm=new WritePSM2Summary();
		CmdLineParser parser = new CmdLineParser(writepsm);
		try {
	    	parser.setUsageWidth(100);
	    	parser.parseArgument(args);
	    	System.err.println("\niPeak v1.0(2013-07)\nWritten by dchaoqin in the\n" +
        				"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    } catch (CmdLineException e) {
	    	System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by duchaoqin in the\n" +
						"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    	System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.org.cn.ipeak.io.WritePSM2Summary [options...]\n");
	    	parser.printUsage(System.err);
	    	System.err.println("\nError:");
	    	System.err.println(e.getMessage());
	    	System.err.println("");
	    	return;
	    }
	    writepsm.transcsv2summary();
	}
	public WritePSM2Summary() {
		
	}
	
	
	public WritePSM2Summary(String psm_csv,String output_psm ) {
		
		this.psm_csv=psm_csv;
		this.output_psm=output_psm;
	}
	public void transcsv2summary() throws IOException {
		CSVReader csvReader=new CSVReader(new FileReader(psm_csv));
		FileWriter writeFile=new FileWriter(output_psm);
		String[] nextLine;
		int lineCounter =0;
		HashMap<String, Integer> headerToColumnMap = new HashMap<String, Integer>();		 
        HashMap<String, Integer> pep_number = new HashMap<String,Integer>();
        
        PSM_number=0;
        writeFile.write("spectrum_index\tSequence\tModification\tCharge\tCalc_mz\tExp_mz\tFdrScore\tProteins\tRank\tIsDecoy\n");
		while ((nextLine = csvReader.readNext()) != null){
		     if (lineCounter == 0) {		    	 
                 for (int i = 0; i < nextLine.length; i++) {
//                	 System.out.println(nextLine[i].trim());
                     headerToColumnMap.put(nextLine[i].trim(), i);
                 }
                 lineCounter++;
             } 
		     else { 
		    	 String spectrumID = nextLine[headerToColumnMap.get("Spectrum ID")];
                 String rank=nextLine[headerToColumnMap.get("rank")];
                 Boolean pass_threshold=Boolean.valueOf( nextLine[headerToColumnMap.get("Pass Threshold")]);
                 String calc_mz= nextLine[headerToColumnMap.get("Calc m/z")];
                 String exp_mz= nextLine[headerToColumnMap.get("Exp m/z")];
                 String charge= nextLine[headerToColumnMap.get("Charge")];
                 String sequence= nextLine[headerToColumnMap.get("Sequence")];
                 String modifications= nextLine[headerToColumnMap.get("Modifications")];
                 String combined_fdrScore;
                 if(headerToColumnMap.containsKey("combined FDRScore")){
                	 combined_fdrScore= nextLine[headerToColumnMap.get("combined FDRScore")];
                 }
                 else{                	 
                	 combined_fdrScore= nextLine[headerToColumnMap.get("FDRScore")];
                 }
                
                 String protein =nextLine[headerToColumnMap.get("proteinacc_start_stop_pre_post_;")];
                 Boolean is_decoy=Boolean.valueOf(nextLine[headerToColumnMap.get("Is decoy")]);    
                 String[] protein_list=protein.split(";");
                 Vector<String> new_protein_list=new Vector<String>();
                 for (int i = 0; i < protein_list.length; i++) {
                	 String the_protein=protein_list[i].replace("\"", "");
                	 int j=the_protein.lastIndexOf("_");
                	 j=the_protein.lastIndexOf("_",j-1);
                	 j=the_protein.lastIndexOf("_",j-1);
                	 j=the_protein.lastIndexOf("_",j-1);
                	 the_protein=the_protein.substring(0,j);
                	 new_protein_list.add(the_protein);                	 
                 }
                 if(is_decoy || !pass_threshold){
                	 continue;
                 }
                 if(!is_decoy){
                	 PSM_number++;
                	 if(pep_number.containsKey(sequence)){
                		pep_number.put(sequence, pep_number.get(sequence)+1); 
                	 }
                	 else{
                		 pep_number.put(sequence, 1);
                	 }
                	 writeFile.write(spectrumID+"\t"+sequence+"\t"+modifications+"\t"+charge+"\t"+calc_mz+"\t"+exp_mz+"\t"
                 +combined_fdrScore+"\t"+new_protein_list+"\t"+rank+"\t"+is_decoy+"\n");  
                	 
                 }
		     }
		}
		csvReader.close();
		writeFile.close();
		peptide_number=pep_number.size();
	}
	public Integer getPSM_number() {
		return PSM_number;
	}
	public Integer getPeptide_number() {
		return peptide_number;
	}
}
