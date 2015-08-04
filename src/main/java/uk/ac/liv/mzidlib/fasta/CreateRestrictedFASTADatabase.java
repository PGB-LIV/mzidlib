package uk.ac.liv.mzidlib.fasta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author jonesar
 * Class to read all PDHs with passthreshold=true and create a new FASTA file from these, assuming that InsertMetaDataFromFasta has already been run to insert
 * sequences and descriptions into the file
 */
public class CreateRestrictedFASTADatabase {
    
    
    //private String inputMzid = "resources/iprg2012_msgf_with_proteins_and_seqs.mzid";
    private String inputMzid = "C:/Work/MassSpec/iprg_2012_archive/MSGF/iprg2012_msgf_no_mods_nd_fdr_threshold_groups_pdh_threshold_seqs.mzid";
    private String outputFasta = "test_outputs/iprg2012_msgf_nd_restricted.fasta";   //default for testing
    
    private Map<String,DBSequence> dbSeqMap = new HashMap<>();
    
    
    public CreateRestrictedFASTADatabase(){
        
    }
    
    public static void main(String args[]){
         CreateRestrictedFASTADatabase restrictedFasta = new CreateRestrictedFASTADatabase();
        if (args.length==2){
            restrictedFasta.inputMzid = args[0];
            restrictedFasta.outputFasta=args[1];
        }
       
        restrictedFasta.init();
    }
    public void use(String inputMzid, String outputFasta){
        
        readProteins(inputMzid,outputFasta);
    }
    
    private void init(){
        
        readProteins(inputMzid,outputFasta);
    }
    
    
    private void readProteins(String inputMzid, String outputFile){
        Writer out = null;
        try {
            out = new BufferedWriter(new FileWriter(outputFile));
            MzIdentMLUnmarshaller mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));
            Iterator<DBSequence> iterDBSeq = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
            while(iterDBSeq.hasNext()){
                DBSequence dbSeq = iterDBSeq.next();
                dbSeqMap.put(dbSeq.getId(), dbSeq);
            }
            Iterator<ProteinDetectionHypothesis> iterPDH = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinDetectionHypothesis);
            Map<String,String> accToDesc = new HashMap<>();
            Map<String,String> accToSeq = new HashMap<>();
            while(iterPDH.hasNext()){
                ProteinDetectionHypothesis pdh = iterPDH.next();
                if(pdh.isPassThreshold()){
                    DBSequence dbSeq = dbSeqMap.get(pdh.getDBSequenceRef());
                    if(dbSeq.getSeq()==null ||dbSeq.getSeq().isEmpty()){
                        System.out.print("Error pdh " + pdh.getId() + " DBSequenceRef has either empty DBSequence or null");
                    }
                    accToSeq.put(dbSeq.getAccession(), dbSeq.getSeq());

                    String desc = "";
                    for(CvParam cvParam : dbSeq.getCvParam()){
                        if(cvParam.getAccession().equals("MS:1001088")){
                            desc = cvParam.getValue();
                        }
                    }
                    accToDesc.put(dbSeq.getAccession(), desc);        

                }
            }
            for(String acc : accToSeq.keySet()){
                
                String seq = accToSeq.get(acc);
                String desc = accToDesc.get(acc);
                if(seq!=null && !seq.equals("")&& desc!=null && !desc.equals("")){
	                out.write(">"+acc+" " + desc+"\n");
	                
	                int i=0;
	                while(i<seq.length()-60){
	                    out.write(seq.substring(i, i+60)+"\n");
	                    i=i+60;
	                }
	                if(i<seq.length()){ //catch last bit of sequence
	                    out.write(seq.substring(i, seq.length())+"\n");           
	                }
                }
            }
            out.close();
        } catch (IOException ex) {
              String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
        } finally {
            try {
                out.close();
            } catch (IOException ex) {
                  String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
            }
        }
            
       
        
        
    }
    
    
}
