package uk.ac.liv.mzidlib.proteogrouper;


import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;


/**
 * Class to be used by the Protein Inference Method
 * @author jonesar
 */
public class PDHSetMember {
    
   // private Boolean isMasterProtein = false;
    //private PDHSetMember subsetOfMember = null;

    private boolean isSameSet = false;
    private boolean isSubset = false;
    private boolean isMultiplySubsumed = false;
   
    
    //private PDHSetMember inSameSetAs = null;                        //Pointer up to Master (if this is not a Master)
    //private PDHSet containingSet = null;                            //Pointer up to containing set (should only be set if this is a Master)
    private ProteinDetectionHypothesis pdh = null;                  //PDH that this set member refers to
    //private ArrayList<String> peptides;
    //private ArrayList<String> uniquePeptides;
    
    private List<PDHSetMember> allContainedMembers = null;     //subset, sameset and multiply subsumed
    private String clusterID = null;                                //Clusters link multiple PAGs by shared peptides
    
    public PDHSetMember(ProteinDetectionHypothesis proteinDetectionHypothesis){
        pdh = proteinDetectionHypothesis;
        //peptides = new ArrayList();
        //uniquePeptides = new ArrayList();
        allContainedMembers = new ArrayList<PDHSetMember>();

    }

    public String getClusterID() {
        return clusterID;
    }

    public void setClusterID(String clusterID) {
        this.clusterID = clusterID;
    }
    
    
    

    public Boolean getIsMultiplySubsumed() {
        return isMultiplySubsumed;
    }

    public void setIsMultiplySubsumed(boolean isMultiplySubsumed) {
        this.isMultiplySubsumed = isMultiplySubsumed;
    }
    
    

    public List<PDHSetMember> getAllContainedMembers() {
        return allContainedMembers;
    }

    public void setAllContainedMembers(List<PDHSetMember> allContainedMembers) {
        this.allContainedMembers = allContainedMembers;
    }
    
    

    public Boolean getIsSameSet() {
        return isSameSet;
    }

    public void setIsSameSet(boolean isSameSet) {
        this.isSameSet = isSameSet;
    }

    public Boolean getIsSubset() {
        return isSubset;
    }

    public void setIsSubset(boolean isSubset) {
        this.isSubset = isSubset;
    }

    
    
    /*
    public PDHSetMember getSubsetOfMember() {
        return subsetOfMember;
    }

    public void setSubsetOfMember(PDHSetMember subsetOfMember) {
        this.subsetOfMember = subsetOfMember;
    }

    public ArrayList<PDHSetMember> getAllSameSetPeptides() {
        return allSameSetPeptides;
    }

    public void setAllSameSetPeptides(ArrayList<PDHSetMember> allSameSetPeptides) {
        this.allSameSetPeptides = allSameSetPeptides;
    }

    public ArrayList<PDHSetMember> getAllSubsetPeptides() {
        return allSubsetPeptides;
    }

    public void setAllSubsetPeptides(ArrayList<PDHSetMember> allSubsetPeptides) {
        this.allSubsetPeptides = allSubsetPeptides;
    }

    public ArrayList<String> getPeptides() {
        return peptides;
    }

    public void setPeptides(ArrayList<String> peptides) {
        this.peptides = peptides;
    }
*/
    /*
    public ArrayList<String> getUniquePeptides() {
        return uniquePeptides;
    }

    public void setUniquePeptides(ArrayList<String> uniquePeptides) {
        this.uniquePeptides = uniquePeptides;
    }
    */
       
    
    

    public ProteinDetectionHypothesis getPdh() {
        return pdh;
    }

    public void setPdh(ProteinDetectionHypothesis pdh) {
        this.pdh = pdh;
    }
    
    /*
    public void setContainingSet(PDHSet belongsToset){
        containingSet = belongsToset;
    }
    
    public PDHSet getContainingSet(){        
        return containingSet;
    }
    */
    /*
    public Boolean getIsMasterProtein() {
        return isMasterProtein;
    }

    public void setIsMasterProtein(Boolean isMasterProtein) {
        this.isMasterProtein = isMasterProtein;
    }

*/
    
    /*
    public PDHSetMember getInSameSetAs() {
        return inSameSetAs;
    }

    public void setInSameSetAs(PDHSetMember inSameSetAs) {
        this.inSameSetAs = inSameSetAs;
    }
*/
    

    
    
}
