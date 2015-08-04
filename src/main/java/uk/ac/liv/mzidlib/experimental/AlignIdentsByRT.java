package uk.ac.liv.mzidlib.experimental;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
/**
 *
 * @author jonesar
 */
public class AlignIdentsByRT {
    
    private static String inputMasterIdents =  "resources/CPTACIdents/mam_042408o_CPTAC_study6_6B011.mzid";
    private static String inputAlignIdents = "resources/CPTACIdents/mam_042408o_CPTAC_study6_6E004.mzid";
    
    //Progenesis master = mam_042408o_CPTAC_study6_6E004
    //Versus progenesis mam_042408o_CPTAC_study6_6B011


    
    private MzIdentMLUnmarshaller unmarshaller;

    private Map<String, Peptide> peptideIdHashMap = new HashMap<>();   
    
    private Map<SpectrumIdentificationItem,Double> masterSIIToRT = new HashMap<>();
    private Map<SpectrumIdentificationItem,Double> alignSIIToRT = new HashMap<>();
    
    private Map<SpectrumIdentificationItem,Double> alignSIIToMappedRT = new HashMap<>();
    
    private Map<Double,Integer> masterRTSortIndex = null;
    private Map<Double,Integer> alignRTSortIndex = null;
    
    private Map<Double,Double> alignRTToDelta = new HashMap<>();           //Map of RTs from the align run, with the known deltas from matched peptides
    //private HashMap<Double,Integer> alignRTsWithDeltasSortIndex =new HashMap();     //Sort order of the RTs that have known deltas mapped to sort index
       
    private Map<String, List<SpectrumIdentificationItem>> masterUniquePepToSII = new HashMap<>();
    private Map<String, List<SpectrumIdentificationItem>> alignUniquePepToSII = new HashMap<>();
    
    private Map<String, SpectrumIdentificationItem> masterUniquePepToMedianSII = null;
    private Map<String, SpectrumIdentificationItem> alignUniquePepToMedianSII = null;
    
    private List<Double> masterRTs = null;
    private List<Double> alignRTs = null;
   
    
    private boolean verbose = true;
    
    private MzidLibUtils utils;
    
    private static final String rtAccession = "MS:1001114";  //     SIR: <cvParam accession="MS:1001114" name="retention time(s)"  cvRef="PSI-MS" value="3102.8816" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
    
    /*
     * For testing only
     */
    public static void main(String args[]){
        
        AlignIdentsByRT align = new AlignIdentsByRT();
        
    }

    public AlignIdentsByRT() {
        
        utils = new MzidLibUtils();
        
        readIdents(true,inputMasterIdents);
        readIdents(false,inputAlignIdents);
        findMedianRTForEachPeptide(true);
        findMedianRTForEachPeptide(false);
        processAlignIdents();
    }
    
     /*
    * Code logic is as follows:
    * 
    * For master idents:
    * 
    * - The concept of an mzid peptide_charge is unique
    * - May be multiple PSMs for the same peptide_charge combination
    * - We probably want to determine the median RT of the all PSMs matched to same peptide_charge entity
    * - Create a HashMap{peptide_charge} ==> arraylist of PSMs (for tracking)
    * - Make arraylist of all median RT values
    * - Sort arrayList
    * - Make HashMap{RT} => sorted index
    * 
    * 
    * 
    * Repeat for align idents:
    * 
    * - For each peptide_charge, calculate median RT
    * - If matches to master, get exact delta and create new RT values for all PSMs, with the delta
    *        
    * - save all with no match for later in an arraylist
    * - Go through no matches, calling routine estimate_rt_delta
    * 
    * estimate_rt_delta: 
    * - Find closest RT in the matched list
    * - Get the position of this RT in the sorted list
    * - Retrieve the deltas from the surrounding PSMs
    * - Calculate median delta, and assign to this peptide and other PSMs
    * 
    */
    
    
    
    private void readIdents(boolean isMaster, String inputFile) {

        try {
            //URL xmlFileURL = JmzIdentMLParser.class.getClassLoader().getResource("Mascot_MSMS_example.mzid");

            if (inputFile != null) {
                
                //HashMap<String, PeptideEvidence> peptideEvidenceIdHashMap = new HashMap();                              //HashMap with PE IDs and PE objects
                                           //HashMap with Peptide IDs and Peptide objects
                Map<String, List<SpectrumIdentificationItem>> peptideEvidenceIDToSIIIDsHashMap = new HashMap<>();     //Pointers from PE ID to SII IDs

                unmarshaller = new MzIdentMLUnmarshaller(new File(inputFile));

                Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
                while (iterPeptide.hasNext()) {
                    Peptide peptide = iterPeptide.next();
                    peptideIdHashMap.put(peptide.getId(), peptide);
                }
                
                Iterator<SpectrumIdentificationResult> iterSIR = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
                while (iterSIR.hasNext()) {
                    SpectrumIdentificationResult sir  = iterSIR.next();
                    double rt = getRTFromSIR(sir);
                    
                    //for(SpectrumIdentificationItem sii : sir.getSpectrumIdentificationItem()){
                        
                        //We only want one SII per Spectrum, since some peptides can cause problems if mapped in multiple similar forms due to I / L differences etc
                        SpectrumIdentificationItem sii = sir.getSpectrumIdentificationItem().get(0);
                        if(isMaster){
                            masterSIIToRT.put(sii, rt);
                        }
                        else{
                            alignSIIToRT.put(sii, rt);
                        }

                        if(sii.getRank() != 1){
                            System.out.println("Error - SII in first position is not top ranked" + sii.getId());
                        }
                        
                        if (sii.isPassThreshold()) {
                            List<PeptideEvidenceRef> peptideEvidenceRefList = sii.getPeptideEvidenceRef();
                            for (int i = 0; i < peptideEvidenceRefList.size(); i++) {

                                Peptide peptide = peptideIdHashMap.get(sii.getPeptideRef());
                                String uniquePep = peptide.getId() + "_" + sii.getChargeState();

                                List<SpectrumIdentificationItem> allSIIsForPep = null;

                                if(isMaster){
                                    if(masterUniquePepToSII.containsKey(uniquePep)){
                                        allSIIsForPep=masterUniquePepToSII.get(uniquePep);
                                    }
                                    else{
                                        allSIIsForPep = new ArrayList<>();
                                        masterUniquePepToSII.put(uniquePep, allSIIsForPep);
                                    }
                                }
                                else{
                                    if(alignUniquePepToSII.containsKey(uniquePep)){
                                        allSIIsForPep=alignUniquePepToSII.get(uniquePep);
                                    }
                                    else{
                                        allSIIsForPep = new ArrayList<>();
                                        alignUniquePepToSII.put(uniquePep, allSIIsForPep);
                                    }
                                }
                                allSIIsForPep.add(sii);
                            }
                        }
                   // }
                }
            } else {
                System.out.println("FILE NOT FOUND:" + inputFile);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /*
     * Helper method to find the median RT value from all PSMs matched to one unique peptide (peptideID_charge)
     * 
     */
    private void findMedianRTForEachPeptide(boolean isMaster){
        Map<String, List<SpectrumIdentificationItem>> uniquePepToSII = null;
        Map<String, SpectrumIdentificationItem> uniquePepToMedianSII = new HashMap<>();
        Map<SpectrumIdentificationItem, Double> siiToRT = null;
        List<Double> allRTs = new ArrayList<>();
        Map<Double,Integer> rtSortIndex = new HashMap<>();
        
        if(isMaster){
           uniquePepToSII  = masterUniquePepToSII;
           siiToRT = masterSIIToRT;
        }
        else{
           uniquePepToSII  = alignUniquePepToSII;
           siiToRT = alignSIIToRT;
        }

        
        for(String pep : uniquePepToSII.keySet()){
            List<SpectrumIdentificationItem> allSII = uniquePepToSII.get(pep);
            SpectrumIdentificationItem medianSII = null;
            
            List<Double> allRT = new ArrayList<>();
            for(SpectrumIdentificationItem sii : allSII){
               allRT.add(siiToRT.get(sii));                
            }
            double rtMedian = CalcMedianWithNearestValue(allRT);
            
            if(allSII.isEmpty()){
                System.out.println("FATAL ERROR, no SIIs for " + pep);
            }
            
            for(SpectrumIdentificationItem sii : allSII){
               //System.out.println(""+siiToRT.get(sii)+ "=" + rtMedian + " (" + pep + ")");
               if(siiToRT.get(sii).equals(rtMedian)){
                   medianSII = sii;
                   break;
               }
            }    
            
            
            if(medianSII!=null){
                uniquePepToMedianSII.put(pep, medianSII);
                if(!allRTs.contains(siiToRT.get(medianSII))){   //TODO - don't want to duplicate RT values, but this could be really slow, optimisation possible
                    allRTs.add(rtMedian);
                }
            }
            else{
                System.out.println("FATAL ERROR, no median SII set for pep" + pep);
            }
        }
        
        //Sort RT values and put sort index into a HashMap
        Collections.sort(allRTs);
        int i = 0;
        for(double rt : allRTs){
            rtSortIndex.put(rt, i);
            i++;
        }
        
        if(isMaster){
            masterUniquePepToMedianSII = uniquePepToMedianSII;
            masterRTs = allRTs;
            masterRTSortIndex = rtSortIndex;
        }
        else{
            alignUniquePepToMedianSII = uniquePepToMedianSII;
            alignRTs = allRTs;
            alignRTSortIndex = rtSortIndex;
        }
    }
    
    
    private void processAlignIdents(){
        
        System.out.println("TempAlign: " + alignUniquePepToMedianSII.size());
        System.out.println("TempMaster: " + masterUniquePepToMedianSII.size());
        
        int siiCounter = 0;
        
                
        List<SpectrumIdentificationItem> unmappedSIIs = new ArrayList<>();
        
        List<Double> allSortedMappedRTs = new ArrayList<>();
        Map<Double, List<SpectrumIdentificationItem>> mappedRTToSII = new HashMap<>();

        
        for(String pep : alignUniquePepToMedianSII.keySet()){
            SpectrumIdentificationItem alignMedianSII = alignUniquePepToMedianSII.get(pep);
            double alignRT = alignSIIToRT.get(alignMedianSII);
            
            //System.out.println("Processing pep: " + pep);
            
            if(masterUniquePepToMedianSII.containsKey(pep)){
                SpectrumIdentificationItem masterMedianSII = masterUniquePepToMedianSII.get(pep);            
                double masterRT = masterSIIToRT.get(masterMedianSII);                
                double delta = masterRT - alignRT;                
                //System.out.println("\tDelta for pep" + pep + " = " + delta);
                alignRTToDelta.put(alignRT, delta);
                alignSIIToMappedRT.put(alignMedianSII, masterRT);           //Give master RT value to this SII
                
                alignMedianSII.getUserParam().add(utils.makeUserParam("mapped RT", ""+masterRT));
                alignMedianSII.getUserParam().add(utils.makeUserParam("mapped RT delta", ""+delta));
                alignMedianSII.getUserParam().add(utils.makeUserParam("RT exact map", "true"));
                
                List<SpectrumIdentificationItem> allSIIsWithThisMappedRT = null;
                if(!mappedRTToSII.containsKey(masterRT)){
                    allSIIsWithThisMappedRT = new ArrayList<SpectrumIdentificationItem>();
                    mappedRTToSII.put(masterRT, allSIIsWithThisMappedRT);
                }
                else{
                    allSIIsWithThisMappedRT = mappedRTToSII.get(masterRT);
                }
                allSIIsWithThisMappedRT.add(alignMedianSII);
                
                
                allSortedMappedRTs.add(masterRT);
                
            }
            else{
                //System.out.println("\tNo match, will need to estimate delta");
                unmappedSIIs.add(alignMedianSII);
                //Double delta = estimateRTDelta()
                
            }
            siiCounter++;
        }
        
        //Have to collect all mapped deltas, before estimating deltas for unmapped SIIs
        for(SpectrumIdentificationItem sii : unmappedSIIs){
             double alignRT = alignSIIToRT.get(sii);
             double estDelta = estimateRTDelta(alignRT);
             double mappedRT = alignRT + estDelta;
             alignSIIToMappedRT.put(sii, mappedRT);
             
             sii.getUserParam().add(utils.makeUserParam("mapped RT", ""+mappedRT));
             sii.getUserParam().add(utils.makeUserParam("mapped RT delta", ""+estDelta));
             sii.getUserParam().add(utils.makeUserParam("RT exact map", "false"));
             allSortedMappedRTs.add(mappedRT);
             
            List<SpectrumIdentificationItem> allSIIsWithThisMappedRT = null;
            if(!mappedRTToSII.containsKey(mappedRT)){
                allSIIsWithThisMappedRT = new ArrayList<SpectrumIdentificationItem>();
                mappedRTToSII.put(mappedRT, allSIIsWithThisMappedRT);
            }
            else{
                allSIIsWithThisMappedRT = mappedRTToSII.get(mappedRT);
            }
            allSIIsWithThisMappedRT.add(sii);
        }

        Collections.sort(allSortedMappedRTs);
        
        
        
        //mean = StatUtils.mean(values, 0, 3);

        //double[] allDeltaValues = new double[alignUniquePepToMedianSII.size()];
        double[] allDeltaValues = new double[10000];            //TODO - this shouldn't be hard-code, see below
        int siiCount=0;
        List<SpectrumIdentificationItem> orderedSIIs = new ArrayList<>();
        
        for(double orderedRT : allSortedMappedRTs){        
            for(SpectrumIdentificationItem sii : mappedRTToSII.get(orderedRT)){
                orderedSIIs.add(sii);
                for(UserParam userParam: sii.getUserParam()){
                    if(userParam.getName().equals("mapped RT delta")){
                        allDeltaValues[siiCount] = Double.parseDouble(userParam.getValue());
                        siiCount++;
                    }
                }
            }
        }
        
        System.out.println("Num mapped SIIs: " + siiCounter);
        System.out.println("Num found SIIs: " + siiCount);      //TODO - fix why these counts are different - unclear
        
        int j = 0;
        
        for(SpectrumIdentificationItem sii : orderedSIIs){      
            Peptide peptide = peptideIdHashMap.get(sii.getPeptideRef());
            String uniquePep = peptide.getId() + "_" + sii.getChargeState();
            
            DescriptiveStatistics stats = new DescriptiveStatistics();
            
            double runningMean = 0.0;
            double runningMedian = 0.0;
            if(j < 5){
                //runningMean = StatUtils.mean(allDeltaValues, 0, 5);
                for( int k = 0; k < 11; k++) {
                   stats.addValue(allDeltaValues[k]);
                }

            }
            else if(j >= siiCount-5){
                //runningMean = StatUtils.mean(allDeltaValues, alignUniquePepToMedianSII.size()-10, alignUniquePepToMedianSII.size()-1);
                for( int k = siiCount-10; k < siiCount; k++) {
                   stats.addValue(allDeltaValues[k]);
                }
            }
            else{
                //runningMean = StatUtils.mean(allDeltaValues, j-5, j+5);                   
                for( int k = j-5; k < j+5; k++) {
                   stats.addValue(allDeltaValues[k]);
                }

            }
            
            runningMean = stats.getMean();
            runningMedian = stats.getPercentile(50);
            j++;
            
            
                        // Add the data from the array


            
            String mappedRT = "";
            String mappedDelta = "";
            String exactMap = "";

            for(UserParam userParam: sii.getUserParam()){
                if(userParam.getName().equals("mapped RT")){
                    mappedRT = userParam.getValue();
                }
                if(userParam.getName().equals("mapped RT delta")){
                    mappedDelta = userParam.getValue();
                }
                if(userParam.getName().equals("RT exact map")){
                    exactMap = userParam.getValue();
                }              

            }

            System.out.println(uniquePep + "\t" + mappedRT + "\t" + mappedDelta + "\t" + exactMap + "\t" + runningMean + "\t" + runningMedian);

        }

        
        
        
        /*
       ArrayList<Double> sortedRTsWithDeltas = new ArrayList(alignRTToDelta.keySet());
       Collections.sort(sortedRTsWithDeltas);
        
        int i = 0;
        for(Double rt : sortedRTsWithDeltas){
            alignRTsWithDeltasSortIndex.put(rt, i);
        }
        */
    }


    /*
     * Helper method to retrieve the RT from an SIR
     * <cvParam accession="MS:1001114" name="retention time(s)"  cvRef="PSI-MS" value="3102.8816" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
     */
    private double getRTFromSIR(SpectrumIdentificationResult sir){
        double rt = -999.0;

        for(CvParam cvParam : sir.getCvParam()){
            if(cvParam.getAccession().equals(rtAccession)){
                rt = Double.parseDouble(cvParam.getValue());
                if(cvParam.getUnitName().contains("minute")){
                    rt /= 60;        //Convert minutes to seconds
                }
            }
        }
        
        if(rt == -999.0){
            System.out.println("FATAL ERROR: No RT for Result: " + sir.getId() + " exiting...");
        }
        
        return rt;
    }
    
    public double estimateRTDelta(double rt){
        
        double estDelta = 0.0;
        
        // TODO: Check if this could return null - AC
        int deltaRTPos = alignRTSortIndex.get(rt);
        
        //Now get nearby deltas        
        int i = deltaRTPos;
        int foundUpwards = 0;
        int foundDownwards = 0;
        
        List<Double> foundDeltas = new ArrayList<>();
        
        //Look upwards for nearby deltas
        while(foundUpwards < 3 && i != alignRTs.size()){
            //System.out.print("i: " + i + " ");
            double nearbyRT = alignRTs.get(i);   
            if(alignRTToDelta.get(nearbyRT) != null){
                foundDeltas.add(alignRTToDelta.get(nearbyRT));
                foundUpwards++;
            }
            i++;
        }
        
        //Reset i to starting point and look downwards for nearby deltas
        i = deltaRTPos;

        while(foundDownwards < 3 && i != -1){
            //System.out.print("i: " + i + " ");
            double nearbyRT = alignRTs.get(i);   
            if(alignRTToDelta.get(nearbyRT) != null){
                foundDeltas.add(alignRTToDelta.get(nearbyRT));
                foundDownwards++;
            }
            i--;
        }
        
        estDelta = CalcMedian(foundDeltas);
        //System.out.println("\tDeltas:" + foundDeltas.toString() + " delta: " + estDelta);
        
        return estDelta;
    }
    
    
    /*
     * Calculate median but for an even set, select the exact value below the median point, rather than taking the mean of tied values
     */
    public static double CalcMedianWithNearestValue(List<Double> values) {
        Collections.sort(values);

        if (values.size() % 2 == 1) //i.e. odd length, so there is a true median
            return values.get((values.size()+1)/2-1);
        else
        {
            double lower = values.get(values.size()/2-1);
            return lower;
            /* //not calculating true median, since we want one value
            double upper = values.get(values.size()/2);

            return (lower + upper) / 2.0;
            * 
            */
        }	
    }
    
    /*
     * Calculate real median
     */
    public static double CalcMedian(List<Double> values) {
        Collections.sort(values);

        if (values.size() % 2 == 1) //i.e. odd length, so there is a true median
            return values.get((values.size()+1)/2-1);
        else
        {
            double lower = values.get(values.size()/2-1);
            double upper = values.get(values.size()/2);
            return (lower + upper) / 2.0;

        }	
    }
   
}
