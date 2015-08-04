/**
 * 
 */
package bgi.ipeak.percolator.paser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import de.proteinms.xtandemparser.interfaces.Modification;
import de.proteinms.xtandemparser.xtandem.Domain;
import de.proteinms.xtandemparser.xtandem.FragmentIon;
import de.proteinms.xtandemparser.xtandem.ModificationMap;
import de.proteinms.xtandemparser.xtandem.Peptide;
import de.proteinms.xtandemparser.xtandem.PeptideMap;
import de.proteinms.xtandemparser.xtandem.Spectrum;
import de.proteinms.xtandemparser.xtandem.SpectrumPeak;
import de.proteinms.xtandemparser.xtandem.SupportData;
import de.proteinms.xtandemparser.xtandem.XTandemFile;
import bgi.ipeak.percolator.util.PSMImpl;
import bgi.ipeak.percolator.util.Peak;

/**
 * @author Administrator
 *
 */
public class XTandemParser {
	/**
	 * 
	 */
	private String decoyRegular="###";	//decoy regular expression
	private String xtandem_xml;		//
	private String output_feature_file;		//feature output file path;
	private Vector<XTandemFeatures> featureList=new Vector<XTandemFeatures>();	//feature list vector
	private Boolean use_ms2=false;

	public XTandemParser(String xtandem_xml,String output_feature_file,String decoyRegex){
		this.xtandem_xml=xtandem_xml;
		this.output_feature_file=output_feature_file;
		this.decoyRegular=decoyRegex;		
	}
	public XTandemParser(){
		
	}
	public void get_feature_file() throws Exception {
		parseXTandemXml();
		writeFeatures();
		clear();
	}
	

	private void parseXTandemXml() throws Exception {
		System.out.println("Parse X!Tandem xml search result "+xtandem_xml+"...");
		XTandemFile xf=new XTandemFile(xtandem_xml);
		de.proteinms.xtandemparser.parser.XTandemParser xp=xf.getXTandemParser();
		boolean byion=xf.getInputParameters().isScoring_bIons() && xf.getInputParameters().isScoring_yIons();
		boolean czion=xf.getInputParameters().isScoring_cIons() && xf.getInputParameters().isScoring_zIons();
		HashMap<String,String> rawPepMap=xp.getRawPeptideMap();
		ModificationMap modMap=xf.getModificationMap();
		String[] varmod=xf.getInputParameters().getResiduePotModMass().split(",");
		HashMap<String, Double> varModifiableSites=new HashMap<String, Double>();
		for(String modif:varmod){
			double mass=Double.parseDouble(modif.split("@")[0]);
			modif=modif.split("@")[1];
			Pattern pattern=Pattern.compile("\\[(\\w+)\\]");
			Matcher matcher=pattern.matcher(modif);
			if(matcher.find()){
				String[] multiSites=matcher.group(1).split("");
				for(String ms:multiSites)
					if(!varModifiableSites.containsKey(ms) && !ms.isEmpty())
						varModifiableSites.put(ms,mass);
			}else{
				if(!varModifiableSites.containsKey(modif))
					varModifiableSites.put(modif,mass);
			}
		}
		ArrayList<Spectrum> specList=xf.getSpectraList();
		PeptideMap pepMap=xf.getPeptideMap();
		for(int i=0;i<specList.size();i++){
			Spectrum specTemp=specList.get(i);
			int specNum=specTemp.getSpectrumNumber();
			XTandemFeatures aFeature=new XTandemFeatures();
			bgi.ipeak.percolator.util.PeptideInfor aPep=new bgi.ipeak.percolator.util.PeptideInfor();
			aFeature.setCharge(Byte.parseByte(String.valueOf(specTemp.getPrecursorCharge())));
			aPep.setCharge(Byte.parseByte(String.valueOf(specTemp.getPrecursorCharge())));
			aFeature.setMaxMatchedIonIntensity(specTemp.getMaxFragIonIntensity());
			aFeature.setRetentionTime(Double.parseDouble(specTemp.getPrecursorRetentionTime().isEmpty()?"0":specTemp.getPrecursorRetentionTime()));
			aFeature.setIndex(specTemp.getSpectrumId()-1);
			
			Peptide bestPeptide=pepMap.getPeptideByIndex(specNum, 1);
			if(bestPeptide==null)
				continue;
			Domain firstDomainBestPep=bestPeptide.getDomains().get(0);
//			System.out.println("parse the "+i+"th spectrum");
			ArrayList<String> proteins=new ArrayList<String>();
			boolean targetDecoy=true;
			Pattern proPattern=Pattern.compile(this.decoyRegular);
			for(int j=0;j<bestPeptide.getDomains().size();j++){
				Domain dm=bestPeptide.getDomains().get(j);
				String pro=dm.getProteinKey().split(" ")[0];
				Matcher match=proPattern.matcher(pro);
				if(match.find() && targetDecoy)
					targetDecoy=false;
				proteins.add(pro);
			}
			aPep.setProteins(proteins);
						
			aFeature.setnLog10Evalue(-Math.log(firstDomainBestPep.getDomainExpect())/Math.log(10));
			aFeature.setScore(firstDomainBestPep.getDomainHyperScore());
			aFeature.setNextScore(firstDomainBestPep.getDomainNextScore());
			aFeature.setDeltaScore(firstDomainBestPep.getDomainHyperScore()-firstDomainBestPep.getDomainNextScore());
			aFeature.setDeltaMass(firstDomainBestPep.getDomainDeltaMh());
			aFeature.setPrecursorMass(specTemp.getPrecursorMh()-firstDomainBestPep.getDomainDeltaMh()-1.007276);
			aFeature.setDeltaMassPPM(firstDomainBestPep.getDomainDeltaMh()/aFeature.getPrecursorMass()*1000000);
			aFeature.setAbsDeltaMass(Math.abs(aFeature.getDeltaMass()));
			aFeature.setAbsDeltaMassPPM(Math.abs(aFeature.getDeltaMassPPM()));
			aPep.setPeptideSequence(firstDomainBestPep.getDomainSequence());
			aPep.setPeptideLength(firstDomainBestPep.getDomainSequence().length());
			double isoDm_0=Math.abs(aFeature.getDeltaMass());
			double isoDm_1=Math.abs(aFeature.getDeltaMass()-1.007825);
			double isoDm_2=Math.abs(aFeature.getDeltaMass()-2*1.007825);
			double isoDm=isoDm_0;
			if(isoDm_1<isoDm)
				isoDm=isoDm_1;
			if(isoDm_2<isoDm)
				isoDm=isoDm_2;
			aFeature.setIsoDeltaMass(isoDm);
			aFeature.setIsoDeltaMassPPM(isoDm/aFeature.getPrecursorMass()*1000000);
			
			@SuppressWarnings("unchecked")
			Vector<FragmentIon[]> fragIonList=xf.getFragmentIonsForPeptide(bestPeptide, 
					firstDomainBestPep, xf.getInputParameters().getSpectrumMonoIsoMassError());
			boolean[][] ionSeries=byion?new boolean[6][firstDomainBestPep.getDomainSequence().length()-1]
			                            :new boolean[2][firstDomainBestPep.getDomainSequence().length()-1];
			double[] ionSeriesInt=byion?new double[12]:new double[4];
			int[] actualIonSeries=byion?new int[12]:new int[4];
			ArrayList<Double> massErrorList=new ArrayList<Double>();
			ArrayList<bgi.ipeak.percolator.util.FragmentIon> fragList=new ArrayList<bgi.ipeak.percolator.util.FragmentIon>();
			double massError=0,fragIonCount=0,totalMatchedIonInt=0;
			for(int j=0;j<fragIonList.size();j++)
	        	for(int k=0;k<fragIonList.get(j).length;k++){
	        		int ionType=fragIonList.get(j)[k].getType();
	        		massError+=fragIonList.get(j)[k].getTheoreticalExperimentalMassError();
	        		massErrorList.add(fragIonList.get(j)[k].getTheoreticalExperimentalMassError());
	        		fragIonCount++;
	        		totalMatchedIonInt+=fragIonList.get(j)[k].getIntensity();
	        		bgi.ipeak.percolator.util.FragmentIon tempFragmentIon=new bgi.ipeak.percolator.util.FragmentIon();
	        		tempFragmentIon.setIntensity(fragIonList.get(j)[k].getIntensity());
	        		tempFragmentIon.setMz(fragIonList.get(j)[k].getMZ());
	        		tempFragmentIon.setID(ionType);
	        		tempFragmentIon.setNumber(fragIonList.get(j)[k].getNumber());
	        		tempFragmentIon.setMassError(fragIonList.get(j)[k].getTheoreticalExperimentalMassError());
	        		tempFragmentIon.setDoubleCharged(fragIonList.get(j)[k].getCharge()==1.0?false:true);
	        		fragList.add(tempFragmentIon);
	        		if(byion){
	        		if(ionType>=3 && ionType<=5)
	        			if(fragIonList.get(j)[k].getCharge()==1.0){
	        				actualIonSeries[ionType-3]++;
	        				ionSeriesInt[ionType-3]+=fragIonList.get(j)[k].getIntensity();
	        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
	        					ionSeries[ionType-3][fragIonList.get(j)[k].getNumber()-1]=true;
	        			}else{
	        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
	        					ionSeries[ionType][fragIonList.get(j)[k].getNumber()-1]=true;
	        				actualIonSeries[ionType]++;
	        				ionSeriesInt[ionType]+=fragIonList.get(j)[k].getIntensity();
	        			}
	        		else if(ionType>=8 && ionType<=10)
	        			if(fragIonList.get(j)[k].getCharge()==1.0){
	        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
	        					ionSeries[ionType-8][aPep.getPeptideLength()-fragIonList.get(j)[k].getNumber()-1]=true;
	        				actualIonSeries[ionType-2]++;
	        				ionSeriesInt[ionType-2]+=fragIonList.get(j)[k].getIntensity();
	        			}else{
	        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
	        					ionSeries[ionType-5][aPep.getPeptideLength()-fragIonList.get(j)[k].getNumber()-1]=true;
	        				actualIonSeries[ionType+1]++;
	        				ionSeriesInt[ionType+1]+=fragIonList.get(j)[k].getIntensity();
	        			}
	        		}else if(czion){
	        			if(ionType==6)
	        				if(fragIonList.get(j)[k].getCharge()==1.0){
		        				actualIonSeries[ionType-6]++;
		        				ionSeriesInt[ionType-6]+=fragIonList.get(j)[k].getIntensity();
		        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
		        					ionSeries[ionType-6][fragIonList.get(j)[k].getNumber()-1]=true;
		        			}else{
		        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
		        					ionSeries[ionType-5][fragIonList.get(j)[k].getNumber()-1]=true;
		        				actualIonSeries[ionType-5]++;
		        				ionSeriesInt[ionType-5]+=fragIonList.get(j)[k].getIntensity();
		        			}
		        		else if(ionType==11)
		        			if(fragIonList.get(j)[k].getCharge()==1.0){
		        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
		        					ionSeries[ionType-11][aPep.getPeptideLength()-fragIonList.get(j)[k].getNumber()-1]=true;
		        				actualIonSeries[ionType-9]++;
		        				ionSeriesInt[ionType-9]+=fragIonList.get(j)[k].getIntensity();
		        			}else{
		        				if(fragIonList.get(j)[k].getNumber()>0 && fragIonList.get(j)[k].getNumber()<aPep.getPeptideLength())
		        					ionSeries[ionType-10][aPep.getPeptideLength()-fragIonList.get(j)[k].getNumber()-1]=true;
		        				actualIonSeries[ionType-8]++;
		        				ionSeriesInt[ionType-8]+=fragIonList.get(j)[k].getIntensity();
		        			}
	        		}
	        	}
			for(int j=0;j<ionSeries.length;j++){
				int longest=1;
				for(int k=0;k<ionSeries[j].length;k++){
					longest=ionSeries[j][k]==false?1:longest+1;
					aFeature.setLongest(aFeature.getLongest()<longest?longest:aFeature.getLongest());
				}
			}
			
			aFeature.setFragDmError(fragIonCount==0?0:massError/fragIonCount);
			/*@SuppressWarnings("static-access")
			InSilicoDigester digester = new InSilicoDigester(bestPeptide, firstDomainBestPep,
					xf.getModificationMap(), xf.getMassesMap(), specTemp.getPrecursorCharge(),
					xf.getInputParameters().getSpectrumMonoIsoMassError());*/
			SupportData supData = xf.getSupportData(bestPeptide.getSpectrumNumber());
	        ArrayList<Double> mzList = supData.getXValuesFragIonMass2Charge();
	        ArrayList<Double> intList = supData.getYValuesFragIonMass2Charge();
	        SpectrumPeak[] peaks = new SpectrumPeak[mzList.size()];
	        Peak[] pList=new Peak[mzList.size()];
	        for (int j = 0; j < mzList.size(); j++) {
	            peaks[j] = new SpectrumPeak();
	            peaks[j].setMz(mzList.get(j));
	            peaks[j].setIntensity(intList.get(j));
	            Peak apeak=new Peak();
	            apeak.setIndex(j);
	            apeak.setIntensity(intList.get(j));
	            apeak.setMass(mzList.get(j));
	            pList[j]=apeak;
	        }
	        aFeature.setPeakArray(pList);
	        aFeature.setLnTotalIntensity(aFeature.calLogTotInt());
	        aFeature.setMatchedIonTotalInte(totalMatchedIonInt==0?0:Math.log(totalMatchedIonInt));
	        aFeature.setRelMatchedIonTotalInte(totalMatchedIonInt/Math.pow(Math.E, aFeature.getLnTotalIntensity()));
	       
	        int[] fracIonSeries=new int[ionSeriesInt.length];
			int[] relIntIonSeries=new int[ionSeriesInt.length];
			for(int k=0;k<ionSeriesInt.length;k++){
				//System.out.println(ionSeriesInt[k]);
				fracIonSeries[k]=(int)(Math.rint(100*actualIonSeries[k]/(double)aPep.getPeptideLength()));
				relIntIonSeries[k]=(int)(aFeature.getLnTotalIntensity()==-1?0:Math.rint(100*ionSeriesInt[k]/Math.pow(Math.E,aFeature.getLnTotalIntensity())));
				//System.out.println(relIntIonSeries[k]);
				//System.exit(0);
			}
			aFeature.setFracIonSeries(fracIonSeries);
			aFeature.setInteFracIonSeries(relIntIonSeries);
			
			double[] massErrorArray=new double[massErrorList.size()==0?2:massErrorList.size()];
			for(int j=0;j<massErrorList.size();j++)
				massErrorArray[j]=massErrorList.get(j);
			Percentile percent=new Percentile();
			aFeature.setFragDmMedian(percent.evaluate(massErrorArray, 50));
			aFeature.setFragDmMedianPPM(percent.evaluate(massErrorArray, 50)/aFeature.getPrecursorMass()*1000000);
			aFeature.setFragDmIqr(percent.evaluate(massErrorArray,75)-percent.evaluate(massErrorArray,25));
			aFeature.setFragDmIqrPPM(aFeature.getFragDmIqr()/aFeature.getPrecursorMass()*1000000);
			
	        /*int matchedBIons=0,theoreticalBIons=0;
	        matchedBIons=digester.getMatchedIons(Ion.B_ION, peaks).size()+digester.getMatchedIons(Ion.BH2O_ION, peaks).size()+
	        			digester.getMatchedIons(Ion.BNH3_ION, peaks).size();
	        theoreticalBIons=digester.getTheoreticIons(Ion.B_ION).length+digester.getTheoreticIons(Ion.BH2O_ION).length+
	        			digester.getTheoreticIons(Ion.BNH3_ION).length;
	        spf.setbIonFraction(matchedBIons/(double)theoreticalBIons);
	        int matchedYIons=0,theoreticalYIons=0;
	        matchedYIons=digester.getMatchedIons(Ion.Y_ION, peaks).size()+digester.getMatchedIons(Ion.YH2O_ION, peaks).size()+
						digester.getMatchedIons(Ion.YNH3_ION, peaks).size();
	        theoreticalYIons=digester.getTheoreticIons(Ion.Y_ION).length+digester.getTheoreticIons(Ion.YH2O_ION).length+
						digester.getTheoreticIons(Ion.YNH3_ION).length;
	        spf.setyIonFraction(matchedYIons/(double)theoreticalYIons);*/
			
			String domainKey=firstDomainBestPep.getDomainKey();
			aFeature.setbIons(rawPepMap.containsKey("b_ions_"+domainKey)?Integer.parseInt(rawPepMap.get("b_ions_"+domainKey)):0);
			aFeature.setyIons(rawPepMap.containsKey("y_ions_"+domainKey)?Integer.parseInt(rawPepMap.get("y_ions_"+domainKey)):0);
			aFeature.setbScore(rawPepMap.containsKey("b_score_"+domainKey)?Double.parseDouble(rawPepMap.get("b_score_"+domainKey)):0);
			aFeature.setyScore(rawPepMap.containsKey("y_score_"+domainKey)?Double.parseDouble(rawPepMap.get("y_score_"+domainKey)):0);
			
			String seq=firstDomainBestPep.getDomainSequence();
			if(seq.startsWith("R") || seq.startsWith("K"))
				aFeature.setEntryN(true);
			if(seq.endsWith("R") || seq.endsWith("K")){
				aFeature.setEntryN(true);
				aFeature.setEntryC(true);
			}
			
			int varModSites=0,modifiable=0;
			ArrayList<Modification> mdl=modMap.getVariableModifications(domainKey);
			if(varModifiableSites.containsKey("["))
				modifiable++;
			if(varModifiableSites.containsKey("]"))
				modifiable++;
			String[] acidsResidue=seq.split("");
			for(String acids:acidsResidue)
				if(varModifiableSites.containsKey(acids))
					modifiable++;
			ArrayList<String> modString=new ArrayList<String>();
			@SuppressWarnings("unused")
			String modSites="";
			for(int j=0;j<mdl.size();j++){
				modSites+=mdl.get(j).getName()+";";
				String residue=mdl.get(j).getName().split("@")[1];
				if(varModifiableSites.containsKey(residue)){
					varModSites++;
					if(!modString.contains(residue))
						modString.add(residue);
				}
			}
			mdl=modMap.getFixedModifications(domainKey);
			for(int j=0;j<mdl.size();j++)
				modSites+=mdl.get(j).getName()+";";
			aFeature.setVarMods(modifiable==0?0:varModSites/(double)modifiable);
			aFeature.setTitle(supData.getFragIonSpectrumDescription());
			
			PSMImpl psm=new PSMImpl(aPep, fragList.toArray(new bgi.ipeak.percolator.util.FragmentIon[0]));
			psm.setMissedCleavages(firstDomainBestPep.getMissedCleavages());
			psm.setNumberOfIonsMatched(fragList.size());
			psm.setPeptide(aPep);
			psm.setPSMId(aFeature.getIndex());
			psm.setRank(1);
			psm.setScore(aFeature.getScore());
			psm.setTarget(targetDecoy);
			aFeature.setiPSM(psm);
			String fileName=new File(xtandem_xml).getName();
			aFeature.setFileName(fileName.substring(0,fileName.lastIndexOf(".")));
			aFeature.setTitle(specTemp.getLabel());
			featureList.add(aFeature);
		}
	}
	
	private void writeFeatures() throws IOException {
		BufferedWriter out=new BufferedWriter(new FileWriter(output_feature_file));
		String titles="PSMid\tTargetDecoy\tHyperScore\tNextScore\tDeltaScore\tMass\tCharge\tDeltaMass" +
				"\tDeltaMassPPM\tabsDM\tabsDMppm\tIsoDM\tIsoDMppm" +
				"\tMissedCleavages\tVarModRatio\tTotalIntensity\tMatchedIonInt\trelTotMatchedIonInt\tMaxMatchedIonInt\t" +
				"BIons\tBIonScore\tYIons\tYIonScore\t";
		if(use_ms2){
			titles+="FragError\tFragDeltaMed\tFragDeltaMedPPM\tFragDeltaIqr\tFragDeltaIqrPPM\t";
		}
		titles+="Longest\tPepLen\tLog10Evalue";
		for(int i=0;i<featureList.get(0).getFracIonSeries().length;i++){
			titles+="\tfracIonSeries_"+i+"\trelIonSeriesInt_"+i;
		}
		titles+="\tPeptide\tProteins";
		out.write(titles);
		out.newLine();
		out.flush();
		for(int i=0;i<featureList.size();i++){
			XTandemFeatures xf=featureList.get(i);
			PSMImpl psm=xf.getiPSM();
			String psmid=xf.getFileName()+":"+xf.getIndex()+"."+psm.getRank();
			out.write(psmid+"\t"+(psm.isTarget()?1:-1)+"\t"+xf.getScore()+"\t"+xf.getNextScore()+"\t"+
					xf.getDeltaScore()+"\t"+xf.getPrecursorMass()+"\t"+xf.getCharge()+"\t"+xf.getDeltaMass()+"\t"+
					xf.getDeltaMassPPM()+"\t"+xf.getAbsDeltaMass()+"\t"+xf.getAbsDeltaMassPPM()+"\t"+
					xf.getIsoDeltaMass()+"\t"+xf.getIsoDeltaMassPPM()+"\t"+psm.getMissedCleavages()+"\t"+
					xf.getVarMods()+"\t"+xf.getLnTotalIntensity()+"\t"+xf.getMatchedIonTotalInte()+"\t"+
					xf.getRelMatchedIonTotalInte()+"\t"+xf.getMaxMatchedIonIntensity()+"\t"+xf.getbIons()+"\t"+
					xf.getbScore()+"\t"+xf.getyIons()+"\t"+xf.getyScore()+"\t");
			if(use_ms2){
				out.write( xf.getFragDmError()+"\t"+xf.getFragDmMedian()+"\t"+xf.getFragDmMedianPPM()+"\t"
						+xf.getFragDmIqr()+"\t"+xf.getFragDmIqrPPM()+"\t");
			}
			out.write(xf.getLongest()+"\t"+psm.getPeptide().getPeptideLength()+"\t"+xf.getnLog10Evalue());
			for(int j=0;j<xf.getFracIonSeries().length;j++)
				out.write("\t"+xf.getFracIonSeries()[j]+"\t"+xf.getInteFracIonSeries()[j]);
			out.write("\t"+psm.getPeptide().getPeptideSequence());
			String pros="";
			for(int j=0;j<psm.getPeptide().getProteins().size();j++)
				pros+="\t"+psm.getPeptide().getProteins().get(j);
			out.write(pros);
			out.newLine();
			out.flush();
		}
		out.close();
	}
	private void clear() {
		featureList.clear();
	}
}
