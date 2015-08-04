package uk.ac.liv.mzidlib.compare;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;

/**
 * 
 * @author Fawaz Ghali 27-Aug-2014 
 */
public class MzidFile {
    
    //metadata
    private AnalysisSoftwareList analysisSoftwareList = new AnalysisSoftwareList();
    private AuditCollection auditCollection = new AuditCollection();
    private Provider provider = new Provider();
    private AnalysisProtocolCollection analysisProtocolCollection = new AnalysisProtocolCollection();
    private CvList cvList = new CvList();
    private AnalysisCollection combinedAnalysisCollection = new AnalysisCollection();
    
    public MzidFile (){
        
    }

}
