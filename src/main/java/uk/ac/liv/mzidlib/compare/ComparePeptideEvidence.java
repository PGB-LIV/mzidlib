package uk.ac.liv.mzidlib.compare;

import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;

/**
 *
 * @author Fawaz Ghali
 */
public class ComparePeptideEvidence {
    
    public static boolean comparePeptideEvidence(PeptideEvidence pe1, PeptideEvidence pe2) {
        
        if(pe1.getDBSequenceRef()!=null&&pe2.getDBSequenceRef()!=null&&!pe1.getDBSequenceRef().equals(pe2.getDBSequenceRef()))
            return false;
        if(pe1.getPeptideRef()!=null&&pe2.getPeptideRef()!=null&&!pe1.getPeptideRef().equals(pe2.getPeptideRef()))
            return false;
        if(pe1.isIsDecoy()!=pe2.isIsDecoy())
            return false;
        if(pe1.getStart()!=null&&pe2.getStart()!=null&&pe1.getStart().intValue()!=pe2.getStart().intValue())
            return false;
        if(pe1.getEnd()!=null&&pe2.getEnd()!=null&&pe1.getEnd().intValue()!=pe2.getEnd().intValue())
            return false;
      
        
        
        return true;
    }
    
}
