package bgi.ipeak.test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.MzIdentML;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class ThresholdPSM {

	@Option(name="-f",required=true,usage="(required) input mzid file path")
	String mzid_path;
	@Option(name="-o",required=true,usage="(required) output mzid file path")
	String out_path;
	private String threshold_score_acc="MS:1001491";
    private Double threshold_score=1.0;
	private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;


	public static void main(String[] args) throws FileNotFoundException {
	
		ThresholdPSM threshold=new ThresholdPSM();
		CmdLineParser parser = new CmdLineParser(threshold);
		try {
	    	parser.setUsageWidth(100);
	    	parser.parseArgument(args);
	    	System.err.println("\niPeak v1.0(2013-11)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n" +
        				"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    } catch (CmdLineException e) {
	    	System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n" +
						"Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
	    	System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.ThresholdPSM [options...]\n");
	    	parser.printUsage(System.err);
	    	System.err.println("\nError:");
	    	System.err.println(e.getMessage());
	    	System.err.println("");
	    	return;
	    } 
		threshold.ThresholdAndPrint();
	}
	public ThresholdPSM(String mzid_file_path,String out_path,Double threshold_score,String threshold_score_acc) {
	
		this.mzid_path=mzid_file_path;
		this.out_path=out_path;
		this.threshold_score_acc=threshold_score_acc;
		this.threshold_score=threshold_score;
	}
	public ThresholdPSM() {
	
	}

	public void ThresholdAndPrint() throws FileNotFoundException {
		System.out.println("Threshold the result: "+mzid_path+"\n"+"by score: "+threshold_score_acc + "value <"+threshold_score+"\n"+
				"output file: "+out_path+"\n");
		mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(mzid_path));
		mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
		DataCollection data=mzIdentML.getDataCollection();
		for (SpectrumIdentificationList spectrumIdentificationList : data.getAnalysisData().getSpectrumIdentificationList()) {
			List<SpectrumIdentificationResult> spe_new_result =new ArrayList<SpectrumIdentificationResult>() ;
			for( SpectrumIdentificationResult spectrumIdentificationResult : spectrumIdentificationList.getSpectrumIdentificationResult()) {
				List<SpectrumIdentificationItem> spe_item =new ArrayList<SpectrumIdentificationItem>();
				for (SpectrumIdentificationItem spectrumIdentificationItem : spectrumIdentificationResult.getSpectrumIdentificationItem()) {
					Boolean pass_item=false;
					for (CvParam cvParam : spectrumIdentificationItem.getCvParam()) {
						if(cvParam.getAccession().equals(threshold_score_acc)){
							if(Double.valueOf(cvParam.getValue())>threshold_score){
								pass_item=true;
							}
						}
					}
					if(!pass_item){
						spe_item.add(spectrumIdentificationItem);
					}
				}
				if(!spe_item.isEmpty()){
					spectrumIdentificationResult.getSpectrumIdentificationItem().clear();
					spectrumIdentificationResult.getSpectrumIdentificationItem().addAll(spe_item);
					spe_new_result.add(spectrumIdentificationResult);
				}				
			}
			spectrumIdentificationList.getSpectrumIdentificationResult().clear();
			spectrumIdentificationList.getSpectrumIdentificationResult().addAll(spe_new_result);
		}
        MzIdentMLMarshaller m = new MzIdentMLMarshaller();        
        m.marshal(mzIdentML, new FileOutputStream(out_path));
	}
	

}
