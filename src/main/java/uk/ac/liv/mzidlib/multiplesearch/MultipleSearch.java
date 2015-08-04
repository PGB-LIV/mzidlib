/**
 * @author Ritesh Krishna , University of Liverpool
 * @date  May 04,2010
 */

package uk.ac.liv.mzidlib.multiplesearch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class MultipleSearch {

	public int decoyRatio;
	
	
	/* The data read from MzIdentML file after parsing */
	public Map<String, List<List<String>>> peptideModificationHash;
	public Map<String, String> peptideIdAndSequenceHash;
	public Map<String, List<List<Object>>> spectrumInformationHash;
	
	 /* For each Spectrum Result id  we need to store :: - peptides associated- scores - evalue	- decoy = true/false  */
	public List<String> spectrumResult = new ArrayList<String>();
	public List<String> peptideNames = new ArrayList<String>();
	public List<Double> evalues = new ArrayList<Double>();
	public List<Double> scores = new ArrayList<Double>();
	public List<String> decoyOrNot = new ArrayList<String>();
	
	/* Once the evalues/scores list is sorted, we need to remember the original indices of the entries (pre-sorted order)
	 so that we can map back that which value belongs to which spectrumResult,peptideNames etc.*/
	public Integer[] sortOrderForEvalues; 

	 /* the above information in sorted order */
	public List<String> sorted_spectrumResult = new ArrayList<String>();
	public List<String> sorted_peptideNames = new ArrayList<String>();
	public List<Double> sorted_evalues = new ArrayList<Double>();
	public List<Double> sorted_scores = new ArrayList<Double>();
	public List<String> sorted_decoyOrNot = new ArrayList<String>();

	/* Store the estimated FDR, q-Value, FDR Score here */
	public List<Double> estimated_simpleFDR = new ArrayList<Double>();
	public List<Double> estimated_qvalue = new ArrayList<Double>();	
	public List<Double> estimated_fdrscore = new ArrayList<Double>();
	
	
	/**
	 * Constructor takes decoy ratio 
	 * @param decoyratio
	 */
	
	MultipleSearch(int decoyratio){
		
		if(decoyratio <= 0)
			decoyratio = 1;
		
		this.decoyRatio = decoyratio;
	}
	/**
	 * Read the mzIdentML file and set the data structures for further processing
	 * @param mzIdFile
	 * @param searchEngine
	 */
	public void readMzIdentMLData(String mzIdFile, String searchEngine){
		
		MzIdentMLReader mzRead = new MzIdentMLReader(mzIdFile,searchEngine); 
		mzRead.parseDocument();
		
		peptideModificationHash = mzRead.getPeptideModification();
		peptideIdAndSequenceHash = mzRead.getPeptideSequences();
		spectrumInformationHash = mzRead.getSpectrumInfo();
		
	}
	
	/**
	 *Overloaded version for Mascot
	 */
	public void readMzIdentMLData(String mzIdFile, String searchEngine,String mascotDecoyIdentifier){
		
		MzIdentMLReader mzRead = new MzIdentMLReader(mzIdFile,searchEngine, mascotDecoyIdentifier); 
		mzRead.parseDocument();
		
		peptideModificationHash = mzRead.getPeptideModification();
		peptideIdAndSequenceHash = mzRead.getPeptideSequences();
		spectrumInformationHash = mzRead.getSpectrumInfo();
		
	}
	
	
	
	/**
	 * Make simple arrayLists from complicated data structures returned after 
	 * reading mzIdentML file. These simple arrayLists will be used for further
	 * calculations of various scores.  
	 */
	private void getEvalueSortedPeptideList(){		
		
		Set<String> spectraID = spectrumInformationHash.keySet();
		Iterator<String> spectraIdIterator = spectraID.iterator();
		
		try{
		while(spectraIdIterator.hasNext()){
			
			String specRes = spectraIdIterator.next();
			List<List<Object>> spectrumRes = spectrumInformationHash.get(specRes); 
			
			for(int i = 0 ; i < spectrumRes.size(); i++){
				
				List<Object> spectrumItem = spectrumRes.get(i);
					
				spectrumResult.add(specRes);
				String pepName = spectrumItem.get(1).toString();
				peptideNames.add(pepName);
				evalues.add(Double.parseDouble(spectrumItem.get(6).toString()));
				scores.add(Double.parseDouble(spectrumItem.get(7).toString()));
				
				String checkForDecoy = new String("false");
				
				if(spectrumItem.size() > 8){
					List<List<Object>> DBRefs = (List<List<Object>>) spectrumItem.get(8);
					
					String decoyCheckForFirstDBRef = new String();
					
					int flag = 0;
					for (int p = 0 ; p < DBRefs.size(); p++){
						List<Object> singleDBRefs = DBRefs.get(p);
						checkForDecoy = singleDBRefs.get(2).toString();
						
						if(p == 0)
							decoyCheckForFirstDBRef = checkForDecoy;
						
						if(!checkForDecoy.equalsIgnoreCase(decoyCheckForFirstDBRef)){
							System.out.println("Different DB REferences for same Spectrum Item have peptide " + specRes + " --" + pepName + "labelled as Decoy and Non-decoy !");
							System.out.println("Assigning the peptide to be a decoy in this case. Press enter to continue....");
							//System.in.read();
							flag = 1;
						}
						
						if (flag == 1){
							checkForDecoy = "false";
							break;
						}
					} 
				} // end of - if( spectrumItem.size() > 8){
				
				decoyOrNot.add(checkForDecoy);
			}
		} // End of while loop
		
		
		 // Call the sorting routine to find the indices of sorted evalues
		TreeSortForIndices sortClass = new TreeSortForIndices();
		sortOrderForEvalues = sortClass.sortTheValueColumn (evalues.toArray(new Double [0]), true); 
		
		// Arrange the values in the order determined by the sorting operation, such that, each index an list can be map to
		// the entries in other lists for the "same" index
		for(int i=0; i <sortOrderForEvalues.length ; i++){
			int index = sortOrderForEvalues[i];
			sorted_spectrumResult.add(spectrumResult.get(index));
			sorted_peptideNames.add(peptideNames.get(index));
			sorted_evalues.add(evalues.get(index));
			sorted_scores.add(scores.get(index));
			sorted_decoyOrNot.add(decoyOrNot.get(index));
		}
		
		// Clear the memory for the items no more needed
		spectrumResult.clear();
		peptideNames.clear();
		evalues.clear();
		scores.clear();
		decoyOrNot.clear();
		
		}catch(Exception e){
			e.printStackTrace();
		}
				
	}
	
	/**
	 * RK 15-08-12
	 * Insert a fake decoy
	 */
	void insertFakeDecoyInTheEnd(){
		
		Double epsilon = 0.1d;
		sorted_spectrumResult.add("Fake_spec");
		sorted_peptideNames.add("Fake Pep");
		sorted_evalues.add(sorted_evalues.get(sorted_evalues.size()-1) + epsilon);
		sorted_scores.add(sorted_scores.get(sorted_scores.size()-1) + epsilon);
		sorted_decoyOrNot.add("true");
		
	}
	/**
	 * RK 15-08-12
	 * Remove the Fake Decoy
	 */
	void removeFakeDecoyFromTheEnd(){
		int index = estimated_fdrscore.size() - 1;
		sorted_spectrumResult.remove(index);
		sorted_peptideNames.remove(index);
		sorted_evalues.remove(index);
		sorted_scores.remove(index);
		sorted_decoyOrNot.remove(index);
		estimated_fdrscore.remove(index);
		estimated_simpleFDR.remove(index);
		estimated_qvalue.remove(index);
	}
	
	// Write the sorted data into a file
	public void writeTheSortedDataToFile(String fileName)throws Exception{
		
		// Basic check to see that each peptide has a evalue
		try{
			if (sorted_peptideNames.size() != sorted_evalues.size()){
				throw (new myDataMissingException("Number of entries = " + sorted_peptideNames.size() +
						"in sorted_peptideNames don't match with the number of entries = "+ sorted_evalues.size() + "in sorted_evalues"));
			}
		}catch(myDataMissingException e){
			System.out.println(e.toString());
			e.printStackTrace();
		}
		
	//	System.out.println("Peptide = " + sorted_peptideNames.size() + " Evalue = " + sorted_evalues.size()
	//			+ " DecoyorNot = " + sorted_decoyOrNot.size() + " Score = " + sorted_scores.size()
	//			+ " Simple FDR = " + estimated_simpleFDR.size() + " Est QVal = " + estimated_qvalue.size()
	//			+ " Est FDR = " + estimated_fdrscore.size());
		//System.in.read();
		
		Writer out = new BufferedWriter(new FileWriter(fileName));
		
		for (int i = 0 ; i < sorted_peptideNames.size(); i++){
			
			String outStr = sorted_spectrumResult.get(i) + "\t" + 
			                sorted_peptideNames.get(i) + "\t" + sorted_decoyOrNot.get(i) + "\t" + 
			                sorted_evalues.get(i).toString() + "\t" + sorted_scores.get(i).toString() + "\t" +
			                estimated_simpleFDR.get(i) + "\t" + estimated_qvalue.get(i) + "\t" +
			                // "\n";
			                estimated_fdrscore.get(i) + "\n";
			
			out.write(outStr);
		}
		out.close();
	}
	
	
	/**
	 *  Compute FDR score and q value etc using the method described in
	 *  Jones et al. Proteomics, 2009,9, 1220-1229 
	 */
	
	public void computeFDRusingJonesMethod() throws Exception{

		// Sort the data using evalue
		getEvalueSortedPeptideList();
		
		insertFakeDecoyInTheEnd(); // RK 15-08-12
		
		// Force the first sorted e-value to be zero.
		//sorted_evalues.set(0, 0f); // RK 21-10-10
		
		computeSimpleFDR(); // Step 1
		computeQValues();   // Step 2
		computeFDRScore();  // Step 3
		
		removeFakeDecoyFromTheEnd(); // RK 15-08-12
	}
	
	/**
	 * Compute simple FDR
	 * 
	 */
	private void computeSimpleFDR(){
		
		int falsePositiveCount = 0;
		 
		for (int i = 0 ; i < sorted_peptideNames.size(); i++){			
			if(sorted_decoyOrNot.get(i).equals("true"))
				falsePositiveCount++;
		
			double falsePositiveDivRatio = (double)falsePositiveCount / (double)decoyRatio;
			
			estimated_simpleFDR.add(falsePositiveDivRatio / (double)(i+1));
			estimated_qvalue.add((double)0);
		}
	}
	
	/**
	 * 	Compute q-value
	 *  
	 */
	private void computeQValues(){
	
		double immediateMinFdr = estimated_simpleFDR.get(estimated_simpleFDR.size() - 1);
		estimated_qvalue.set(estimated_qvalue.size()-1, immediateMinFdr);
		for (int i = estimated_simpleFDR.size() - 1 ; i > 0 ; i--){
			double currentFDR = estimated_simpleFDR.get(i - 1);
			
			if( currentFDR < immediateMinFdr)
				immediateMinFdr = currentFDR;
			
			estimated_qvalue.set(i-1, immediateMinFdr);
		}
	}
	
	/**
	 * Compute FDR Score
	 *  
	 */
	private void computeFDRScore() throws Exception{
		
		// Initialize the estimated_fdrscore by adding elements containing 0. This needs to be
		// done because of back tracking involved in the algorithm, so simple .add() wouldn't work
		for (int i = 0 ; i < sorted_peptideNames.size(); i++){
			estimated_fdrscore.add((double)0);
		}
		
		//Double prev_evalue = sorted_evalues.get(0); 			// RK 21-10-10
		//Double prev_qvalue = estimated_qvalue.get(0);
		//Double prev_prev_evalue =  sorted_evalues.get(0); // previous to previous evalue in case of straight vertical rise in q value without any change in evalue
		
		double prev_evalue = 0d;	    // RK 21-10-10
		double prev_qvalue = 0d;     // RK 21-10-10
		double prev_prev_evalue = 0d;// RK 21-10-10 
		
		int counter_backwardStep = 0; // RK 21-10-10
		//int counter_backwardStep = 1; // First time, the counter is 1 to account for the origin
		
		//int i = 1;// RK 21-10-10
		int i = 0; // RK 21-10-10
		for (; i < sorted_peptideNames.size(); i++){
			double current_evalue = sorted_evalues.get(i);
			double current_qvalue = estimated_qvalue.get(i);
			
			if (current_qvalue > prev_qvalue){
				
				double slope;
				double intercept;
				
				int id = 0;
				if(current_evalue != prev_evalue){
					slope = (current_qvalue - prev_qvalue)/(current_evalue - prev_evalue);
					id  = 1; //RK 22-10-10
				}else{
					slope = (current_qvalue - prev_qvalue)/(current_evalue - prev_prev_evalue);
					id = 2; //RK 22-10-10
				}
				
				//intercept = current_qvalue - slope * current_evalue; //RK 22-10-10
				//RK 22-10-10 - Tells us which co-ordinates to use for calculating intercepts.
				if(id == 1) 
					intercept = prev_qvalue - slope * prev_evalue;
				else
					intercept = prev_qvalue - slope * prev_prev_evalue;
				
				if (counter_backwardStep > 0){ // compute the FDR score for flat q-value region
					for (int k = 0 ;k <= counter_backwardStep ; k++){
						int index = i - counter_backwardStep + k;
						//System.out.println("i = " + i + "Count = " + counter_backwardStep +  " k = " + k + " index = " + index);
						double fdrScore = slope * sorted_evalues.get(index) + intercept;
						estimated_fdrscore.set(index, fdrScore);
						//System.out.println("i = " + i + "Count = " + counter_backwardStep +  " k = " + k + " index::1 = " + index + " slope = " + slope + " intercept = " + intercept + " e-val = " + sorted_evalues.get(index) + " FDR = " + fdrScore);
					}
				}else{ 							// In case an immediate increment in q value is found
					double fdrScore = slope * current_evalue + intercept;
					estimated_fdrscore.set(i, fdrScore);
					//System.out.println("index::2 i = " + i + " FDR = " + fdrScore);
				}
				
				// Re-initialise the variables
				counter_backwardStep = 0;
				if(current_evalue > prev_evalue){ // the previous e-value will change only if the current e-value is different
					prev_prev_evalue = prev_evalue;
					prev_evalue = current_evalue;
				}
				prev_qvalue = current_qvalue;
				
			}else{
				counter_backwardStep++;
			}
		}
		
		// In case if we miss to update the very last values in estimated_fdrscore because 
		// current_qvalue == prev_qvalue
		if(estimated_fdrscore.get(i-1) == 0){
			double lastFdrValue = 0;
			i = i - 1;
			while(estimated_fdrscore.get(i) == 0 && i > 0){
				i--;
			}
			if(i == 0){
				System.out.println("\n Can't compute FDR. Likely that all the PSM hits are corrects, so nothing to do.");
			}else{
				lastFdrValue = estimated_fdrscore.get(i);
				while(i < estimated_fdrscore.size()){
					estimated_fdrscore.set(i, lastFdrValue);
					i++;
				}
			}
		}
		
	}
	
	/*************************************************************************/
	/**
	 * Extract all this data from XML 
	 */
	public Map<String, List<List<String>>> getFromXMLPeptideModificationHash(){
		return peptideModificationHash;
	}
	public Map<String, String> getFromXMLPeptideSequenceHash(){
		return peptideIdAndSequenceHash;
	}
	public Map<String, List<List<Object>>> getFromXMLSpectrumInfoHash(){
		return spectrumInformationHash;
	}

	/**
	 * Return the computed result structure. The elements in all the following 
	 * lists are arranged to be in correspondence with each other according
	 * to their indices. Eg - sorted_spectrumResult[i] has sorted_evalues[i] 
	 * and estimated_simpleFDR[i] etc. etc. We can pull whatever information
	 *  we want individually and create a matrix like structure when needed, for 
	 *  further processing.     
	 */
	public List<String> getSorted_spectrumResult(){
			return sorted_spectrumResult;
	}
	
	public List<String> getSorted_peptideNames(){ 
			return sorted_peptideNames;
	}
	
	public List<Double> getSorted_evalues(){
		return sorted_evalues;
	}
	
	public List<Double> getSorted_scores(){
			return sorted_scores;
	}
	
	public List<String> getSorted_decoyOrNot(){ 
			return sorted_decoyOrNot;
	}
	
	public List<Double> getSorted_simpleFDR(){
			return estimated_simpleFDR;
	}
	
	public List<Double> getSorted_qValues(){ 
			return estimated_qvalue;	
	}
	
	public List<Double> getSorted_estimatedFDR(){
		return estimated_fdrscore ;
	}
	
	
	/**
	 * Clear all the data structures in the class to release all the memory
	 */
	public void clearAllData(){
		peptideModificationHash.clear();
		peptideIdAndSequenceHash.clear();
		spectrumInformationHash.clear();
		sorted_spectrumResult.clear();
		sorted_peptideNames.clear();
		sorted_evalues.clear();
		sorted_scores.clear();
		sorted_decoyOrNot.clear();
		estimated_simpleFDR.clear();
		estimated_qvalue.clear();	
		estimated_fdrscore.clear();
	}
	/*************************************************************************/
	
	/**
	 * Test the functions implemented in the class....
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//String xmlToRead 	= args[0];
		//String searchEngine = args[1];
		//String outFileName 	= args[2];
		//int decoyRatio = Integer.parseInt(args[3]);

		
		//String xmlToRead = "exampleMzIDFiles/125merge_omssa.mzid";
		//String searchEngine = "omssa";
		//String outFileName = "output/sortedEvalues_125_omssa.txt";
		
		//String xmlToRead = "exampleMzIDFiles/Mascot_SS_exp13.mzid";
		//String searchEngine = "mascot";
		//String outFileName = "output/sortedEvalues_mascot.txt";
		String mascotDecoyTag  = "Rev";
		
		String xmlToRead = "exampleMzIDFiles/From_tandem_448.mzid";
		String searchEngine = "X!Tandem";
		String outFileName = "output/sortedEvaluesFrom_tandem_448.txt";
		
		///////////////////////////////////////////////////////////////
		//String xmlToRead = "exampleMzIDFiles/Toxo1D_Slice10_2_omssa.mzid";
		//String searchEngine = "omssa";
		//String outFileName = "output/sortedEvalues_ToxoSlice10_omssa.txt";
		
		//String xmlToRead = "exampleMzIDFiles/Toxo1D_Slice10_2_mascot.mzid";
		//String searchEngine = "mascot";
		//String outFileName = "output/sortedEvalues_ToxoSlice10_mascot.txt";
		//String mascotDecoyTag  = "REV";
		
		//String xmlToRead = "exampleMzIDFiles/Toxo1D_Slice10_2_tandem.mzid";
		//String searchEngine = "X!Tandem";
		//String outFileName = "output/sortedEvalues_ToxoSlice10_tandem.txt";
		
		int decoyRatio = 1;
		/////////////////////////////////////////////////////////////////

		// Read the MzIdentML file and collect relevant information 	
		MultipleSearch m = new MultipleSearch(decoyRatio);
		
		if(searchEngine.equals("mascot"))
			m.readMzIdentMLData(xmlToRead, searchEngine,mascotDecoyTag);
		else
			m.readMzIdentMLData(xmlToRead, searchEngine);
				
		m.computeFDRusingJonesMethod();
		
		// Write the sorted data to a file
		m.writeTheSortedDataToFile(outFileName);
		
		
		
	}
	

}
