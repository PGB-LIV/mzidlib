package bgi.ipeak.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.percolator.PlusPropertiesToMzid;

public class SplitPercolatorResult {
	@Option(name="-perfile",required=true,usage="(required) features file")
	private String perfile="";
	@Option(name="-outdir",required=false,usage="(optional) features file")
	private String outdir="";
	@Option(name="-mziddir",required=true,usage="(required) features file")
	private String mziddir="";
	@Option(name="-s",required=true,usage="(required) The search engine code.[1:MSGF+(default),2:OMSSA,3:X!Tandem]")
	private int searchEngine=1;
	
	public static void main(String[] args) throws IOException {
	
		SplitPercolatorResult spp=new SplitPercolatorResult();
		CmdLineParser parser = new CmdLineParser(spp);
		try {
	    	parser.setUsageWidth(100);
	    	parser.parseArgument(args);
	    	System.err.println("\niPeak v1.0(2013-11)\nWritten by  in the\n" +
        				"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    } catch (CmdLineException e) {
	    	System.err.println("\niPeakPercolator v1.0(2013-11)\nWritten by  in the\n" +
						"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    	System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.org.cn.ipeak.core.IPeakPercolator [options...]\n");
	    	parser.printUsage(System.err);
	    	System.err.println("\nError:");
	    	System.err.println(e.getMessage());
	    	System.err.println("");
	    	return;
	    }
		spp.split_percolatorResult();
	}

	private void split_percolatorResult() throws IOException {
		File per_file = new File(perfile);
		if(outdir.equals("")){
			outdir=per_file.getParent()+"/";
		}
		if(outdir.lastIndexOf("/") != outdir.length()-1){
			outdir+="/";
		}
		BufferedReader read =new BufferedReader(new FileReader(perfile));
		String line=read.readLine();
		String title = line;
		HashMap<String,Vector<String> > per_file_list=new HashMap<String, Vector<String>>();
		while((line=read.readLine())!=null){
			String[] items=line.split(":");
			String id= outdir+items[0]+".per.txt";
			if(per_file_list.containsKey(id)){
				Vector<String> infor =per_file_list.get(id);
				infor.add(line);
				per_file_list.put(id, infor);
			}
			else{
				Vector<String> infor =new Vector<String>();
				infor.add(line);
				per_file_list.put(id, infor);
			}						
		}
		read.close();
		Set<String> out_names= per_file_list.keySet();
		for (String out : out_names) {
			BufferedWriter writer=new BufferedWriter(new FileWriter(out));
			writer.append(title);
			writer.newLine();
			for (String line_infor : per_file_list.get(out)) {
				writer.append(line_infor);
				writer.newLine();
			}
			writer.close();
		}
		System.out.println("split *.per.txt file ok!");
		addPercolator2mzidP(out_names);
		
	}
	private void addPercolator2mzidP(Set<String> per_name) throws IOException {
		System.out.println("Add percolator score to mzid file");
		if(mziddir.lastIndexOf("/")!=mziddir.length()-1){
			mziddir=mziddir+"/";
		}
		String tag = "";
		String tag_out="";
		if(searchEngine==1){
			tag="_msgf.mzid";
			tag_out="_msgfAddP.mzid";
		}
		else if(searchEngine==2){
			tag="_omssa.mzid";
			tag_out="_omssaAddP.mzid";
		}
		else if(searchEngine==3){
			tag="_xtandem.mzid";
			tag_out="_xtandemAddP.mzid";
		}
		
		Vector<String> addp_mzdi_list=new Vector<String>();
		for (String file_name : per_name) {
			File per_file =new File(file_name);
			String name=per_file.getName().substring(0,per_file.getName().lastIndexOf(".per.txt"));
			String mzid_path=mziddir+name+tag;
			String mzidP_path=outdir+name+tag_out;
			File mz_file=new File(mzid_path);
			System.out.println("plus "+file_name +" to mzid file "+mzid_path);
			if(mz_file.exists()){				
				PlusPropertiesToMzid pmp=new PlusPropertiesToMzid(mzid_path,file_name ,mzidP_path);
				pmp.export();
			}
			else{
				System.err.println("can not found mzid file: "+mzid_path);
			}
			addp_mzdi_list.add(mzidP_path);
		}
//		Combined_mzid(addp_mzdi_list);
	}
}
