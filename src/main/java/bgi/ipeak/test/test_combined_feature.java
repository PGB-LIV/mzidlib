package bgi.ipeak.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;

public class test_combined_feature {
	
	private Vector<String> feature_files;
	private String out_file;
	public static void main(String[] args) throws IOException {
		test_combined_feature spp=new test_combined_feature();
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
		spp.combine_feature();
	}
	private void combine_feature() throws IOException {
		out_file="/ifs1/ST_PRO/USER/duchaoqin/project/mulpitle_search/last_paper/sub_px0071_1/test_com/ipeak_analyse/xtandem/TITLE_201B7_xtandem.features";
		feature_files=new Vector<String>();
		feature_files.add("est_com/ipeak_analyse/xtandem/TITLE_201B7_11.features");
		feature_files.add("test_com/ipeak_analyse/xtandem/TITLE_201B7_9.features");
		combined_features(feature_files, out_file);
	}
	private void combined_features(Vector<String > feature_file_list, String all_feature_file) throws IOException {
	
		File test_file=new File(all_feature_file);
		if(test_file.exists()){
			test_file.delete();
		}		
		BufferedReader readtitle=new BufferedReader(new FileReader(feature_file_list.get(0)));
		String title=readtitle.readLine();
		readtitle.close();
		Vector<String> featureinfors=new Vector<String>();
		for (String feature_file:feature_file_list){
			System.out.println(feature_file);
			BufferedReader read=new BufferedReader(new FileReader(feature_file));
			String line=read.readLine();
			while((line=read.readLine())!=null){
				featureinfors.add(line);
			}
			read.close();
		}
		
		BufferedWriter writeFile=new BufferedWriter(new FileWriter(all_feature_file));
		writeFile.write(title);
		for (String infor : featureinfors) {
			writeFile.append(infor);
			writeFile.newLine();
		}
		writeFile.close();
	}
}
