package uk.ac.liv.mzidlib.multiplesearch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FdrAndMzIdentInformationContainer {
	
	String xmlFileName;
	String searchEngine;
	Map<String, List<List<String>>> peptideModifications;
	Map<String, String> peptideSequence;
	Map<String, List<List<Object>>> spectrumInfo;
	List<String> sorted_spectrumID;
	List<String> sorted_peptideID;
	List<Double> sorted_evalues; 
	List<Double> sorted_scores;
	List<String> sorted_decoyornot;
	List<Double> sortedSimpleFDR;
	List<Double> sorted_qValues; 
	List<Double> sorted_estimatedFDR;
	
	FdrAndMzIdentInformationContainer(){
		xmlFileName  = new String("");
		searchEngine = new String("");
		peptideModifications = new HashMap <String, List<List<String>>>(1);
		peptideSequence = new HashMap <String, String>(1) ;
		spectrumInfo = new HashMap <String, List<List<Object>>>(1) ;
		sorted_spectrumID = new ArrayList <String>(1);
		sorted_peptideID = new ArrayList <String>(1);
		sorted_evalues = new ArrayList <Double> (1); 
		sorted_scores = new ArrayList <Double>(1);
		sorted_decoyornot = new ArrayList <String> (1);
		sortedSimpleFDR = new ArrayList <Double>(1);
		sorted_qValues = new ArrayList <Double>(1); 
		sorted_estimatedFDR = new ArrayList <Double>(1);
	}
	
	

	public void populateData(String xmlToRead,String searchEngine,
			Map <String, List<List<String>>> pepMod,
			Map <String, String> pepSeq,
			Map <String, List<List<Object>>> specInfo,
			List <String> sorted_spec,List <String> sorted_pepID,
			List <Double> sorted_eval, List <Double> sorted_score,
			List <String> sorted_decoy,List <Double> sorted_FDR,
			List <Double> sorted_qValues, List <Double> sorted_estFDR){
		
		this.xmlFileName 		= xmlToRead;
		this.searchEngine 		= searchEngine;
		this.peptideModifications 	= new HashMap<String, List<List<String>>>(pepMod);
		this.peptideSequence 		= new HashMap <String, String>(pepSeq);
		this.spectrumInfo 			= new HashMap <String, List<List<Object>>>(specInfo) ;
		this.sorted_spectrumID 		= new ArrayList <String>(sorted_spec);
		this.sorted_peptideID 		= new ArrayList <String>(sorted_pepID);
		this.sorted_evalues 		= new ArrayList <Double>(sorted_eval); 
		this.sorted_scores 			= new ArrayList <Double>(sorted_score);
		this.sorted_decoyornot 		= new ArrayList <String>(sorted_decoy);
		this.sortedSimpleFDR 		= new ArrayList <Double>(sorted_FDR);
		this.sorted_qValues         = new ArrayList <Double>(sorted_qValues); 
		this.sorted_estimatedFDR    = new ArrayList <Double>(sorted_estFDR);
	}
	
}
