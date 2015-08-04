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

public class WriteProteins2Summary {
	@Option(name="-f",required=true,usage="(required) The input psm csv file")
	private String pro_csv;
	@Option(name="-o",required=true,usage="(required) The output psm summary file")
	private String output_pro;
	private Integer protein_number;
	
	public static void main(String[] args) throws IOException {
		WriteProteins2Summary writeprotein=new WriteProteins2Summary();
		CmdLineParser parser = new CmdLineParser(writeprotein);
		try {
	    	parser.setUsageWidth(100);
	    	parser.parseArgument(args);
	    	System.err.println("\niPeak v1.0(2013-07)\nWritten by duchaoqin (duchaoqin@genomics.cn) in the\n" +
        				"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    } catch (CmdLineException e) {
	    	System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by duchaoiqn (duchaoqin@genomics.cn) in the\n" +
						"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    	System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.org.cn.ipeak.io.WriteProteins2Summary [options...]\n");
	    	parser.printUsage(System.err);
	    	System.err.println("\nError:");
	    	System.err.println(e.getMessage());
	    	System.err.println("");
	    	return;
	    }
	    writeprotein.transcsv2summary();
	}
	public WriteProteins2Summary() {
		// TODO �Զ���ɵĹ��캯����
	}
	public WriteProteins2Summary(String pro_csv,String output_pro) {
		// TODO �Զ���ɵĹ��캯����
		this.pro_csv=pro_csv;
		this.output_pro=output_pro;
	}
	public void transcsv2summary() throws IOException {
		CSVReader csvReader=new CSVReader(new FileReader(pro_csv));
		FileWriter writeFile=new FileWriter(output_pro);
		String[] nextLine;
		int lineCounter =0;
		HashMap<String, Integer> headerToColumnMap = new HashMap<String, Integer>();		 
        HashMap<Integer, String> columnToHeaderMap = new HashMap<Integer, String>();
        
        String pag_id_this="";
        protein_number=0;
        Vector<String> same_set_protein = new Vector<String>();
        
        String protein_accession="";
        String description="";
        String peptide_number="";
        String proteoGrouper_PDHscore="";
        String razor_peptides = "";
        String uniq_peptides = "";
        writeFile.write("PAG_ID\tProtein_Accession\tProteoGrouper_PDHscore\tPeptideNumber\tUniq_Peptides\tRazorPeptides\tAnnotation\tSame_set_protein\n");
		while ((nextLine = csvReader.readNext()) != null){
		
		     if (lineCounter == 0) {		    	 
                 for (int i = 0; i < nextLine.length; i++) {
//                	 System.out.println(nextLine[i].trim());
                     headerToColumnMap.put(nextLine[i].trim(), i);
                     columnToHeaderMap.put(i, nextLine[i]);
                 }
                 lineCounter++;
             }
		     else {
		    	
		    	 Vector<String> inforStrings=new Vector<String>();
		    	 for (Integer i =0;i<nextLine.length;i++) {
					if(nextLine[i].startsWith("\"")){
						Integer j ;
						if(nextLine[i].lastIndexOf("\"")==nextLine[i].length()-1){
							inforStrings.add(nextLine[i]);
						}
						else{
							for (j=i+1;j<nextLine.length;j++){								
								if(nextLine[i].lastIndexOf("\"")==nextLine[i].length()-1){
									nextLine[i]+=nextLine[j];
									inforStrings.add(nextLine[i]);
									continue;
								}
								else{
									nextLine[i]+=nextLine[j];
								}
							}
							i=j;
						}
					}
					else{
						inforStrings.add(nextLine[i]);
					}					
				}
		    	 
		    	 String pag_id=inforStrings.get(headerToColumnMap.get("PAG ID")); 
		    	 
		   
		    	 if(!pag_id.equals(pag_id_this)){
//		    		 System.out.println(pag_id+"\t"+pag_id_this);
		    		 if(!pag_id_this.equals("")){
			    		 writeFile.write(pag_id_this+"\t"+protein_accession+"\t"+proteoGrouper_PDHscore+"\t"+peptide_number+"\t"
		    				 	+uniq_peptides+"\t"+razor_peptides+"\t"+description+"\t"+same_set_protein+"\n");
			    		 protein_number++;
		    		 }	    		 
		    		 same_set_protein.clear();
		    		 pag_id_this=pag_id;		    		 
		    	 }
		    	 
		    	 String group_membership=inforStrings.get(headerToColumnMap.get("group membership"));
		    	 Boolean is_decoy=Boolean.valueOf(inforStrings.get(headerToColumnMap.get("Is decoy")));
//		    	 Boolean pass_threshold=Boolean.valueOf(inforStrings.get(headerToColumnMap.get("protein accession")]);
		    	 
		    
		    	 if(is_decoy){
		    		 continue;
		    	 }
   	 
		    	 if(group_membership.contains("anchor protein")){		    		 
		    		 protein_accession=inforStrings.get(headerToColumnMap.get("protein accession"));		    		 
		    		 description=inforStrings.get(headerToColumnMap.get("description"));
		    		 peptide_number=inforStrings.get(headerToColumnMap.get("distinct peptide sequences"));
		    		 proteoGrouper_PDHscore=inforStrings.get(headerToColumnMap.get("ProGrouper:PDH score"));
		    		 uniq_peptides=inforStrings.get(headerToColumnMap.get("unique peptides"));
		    		 razor_peptides=inforStrings.get(headerToColumnMap.get("razor peptides"));
		    		 if(uniq_peptides.equals("")){
		    			 uniq_peptides="null";
		    		 }
		    		 if(razor_peptides.equals("")){
		    			 razor_peptides="null";
		    		 }
//		    		 System.out.println("anch\t"+protein_accession+"\t"+razor_peptides+"\n");
		    	 }
		    	 else if(group_membership.contains("same-set")){
		    		 String the_protein_accession=inforStrings.get(headerToColumnMap.get("protein accession"));
		    		 same_set_protein.add(the_protein_accession);
		    	 }
		    	 else if(group_membership.contains("subsumable")){
//		    		 String the_protein_accession=inforStrings.get(headerToColumnMap.get("protein accession")];
		    		 continue;
		    	 }
		    	 else if(group_membership.contains("sub-set")){
//		    		 String the_protein_accession=inforStrings.get(headerToColumnMap.get("protein accession")];
		    		 continue;
		    	 }
		     }
		}
		csvReader.close();
		writeFile.close();
	}
	public Integer getProtein_number() {
		return protein_number;
	}
}


