package uk.ac.liv.mzidlib.fasta;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import uk.ac.ebi.jmzidml.model.mzidml.MzIdentML;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 *
 * @author jonesar
 */
public class InsertMetaDataFromFasta {
    
    
    private Map<String, String> accToSeq = new HashMap<>();
    private Map<String, String> accToDefline = new HashMap<>();
    
    private String inputMzid = "example_files/Toxo-1D_OF.mzid";
    private String inputFasta = "example_files/TgondiiME49_ToxoDB-6_2.fasta";   //default for testing
    private String outputFile = "edited_file_from_fasta.mzid";   //default for testing
    private String accessionRegex = " ";                        //Separate accessions from defline by a space by default
    
    private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    
    private Cv psiCV;
    private Cv unitCV;
    private MzidLibUtils mzidLibUtils;
    
    public static void main(String args[]){
        InsertMetaDataFromFasta insertMD = new InsertMetaDataFromFasta();
    }
    
    /*
     * Default constructor for testing with defaults only
     */
    private InsertMetaDataFromFasta(){
        init();
    }
    
    public InsertMetaDataFromFasta(String mzidIn, String mzidOut, String fastaIn, String accRegex){
        this.inputMzid = mzidIn;
        this.outputFile = mzidOut;
        this.inputFasta = fastaIn;
        this.accessionRegex = accRegex;
        
        init();
    }
    
    private void init(){
        mzidLibUtils = new MzidLibUtils();
        readFasta();
        insertDetailsIntoMzid();
        writeToMzIdentMLFile();
    }
    
    private void readFasta(){
        InputStream fstream = null;
        try {
            fstream = new FileInputStream(inputFasta);
            // Get the object of DataInputStream
            InputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String line;
            String currSequence = "";
            String currProtAcc = null;
            String currDefline = null;
            int recordCounter = 0;
            while ((line = br.readLine()) != null)   {          
                line = line.replaceAll("\n","");                
                      
                if(line.contains(">")){                
                    //Insert previous into hash and reset
                    if(recordCounter != 0){              
                        currSequence = currSequence.replaceAll(" ","");
                        accToSeq.put(currProtAcc, currSequence);
                        //System.out.println("Inserting:" + currProtAcc + "_" + currSequence);
                        accToDefline.put(currProtAcc, currDefline);
                        //System.out.println("Inserting2:" + currProtAcc + "_" + currDefline);
                        
                        currSequence = "";
                    }   
                    
                    line = line.replaceAll(">", "");
                    
                    int splitPos = line.indexOf(accessionRegex);
                    if(splitPos!=-1){
                        currProtAcc = line.substring(0,splitPos);
                        currDefline = line.substring(splitPos+1);
                    }
                    else{
                        System.out.println("Regular expression not found for split: " + line + " regex:" + accessionRegex);
                    }
                    
                    recordCounter++;
                }
                else{
                    currSequence += line;
                }
            }
            //handle last
            accToSeq.put(currProtAcc, currSequence.replaceAll(" ",""));
            accToDefline.put(currProtAcc, currDefline);
            //Close the input stream
            in.close();
        } catch (FileNotFoundException ex) {
              String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 01 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
        } catch (IOException ex) {
              String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
        } finally {
            try {
                fstream.close();
            } catch (IOException ex) {
                  String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
            }
        }
       
        
    }
    
    private void insertDetailsIntoMzid(){
         try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));
            mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
            
            Iterator<Cv> iterCv = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()){
                Cv cv = iterCv.next();
                if(cv.getUri().toLowerCase().contains("psi")){
                    psiCV = cv;
                }
                else if(cv.getUri().toLowerCase().contains("unit")){
                    unitCV = cv;
                }
            }
            
            for(DBSequence dbSeq : mzIdentML.getSequenceCollection().getDBSequence()){
                
                //System.out.println("Lookup: " + dbSeq.getAccession());
                String seq = accToSeq.get(dbSeq.getAccession());
                //System.out.println("seq: " + seq);
                if(seq!=null){
                    dbSeq.setSeq(seq);
                }
                
                String defLine = accToDefline.get(dbSeq.getAccession());
                if(defLine!=null){
                    boolean alreadyPresent = false;
                    for(CvParam cvParam : dbSeq.getCvParam()){
                        if(cvParam.getAccession().equals("MS:1001088")){
                            alreadyPresent = true;
                            cvParam.setValue(defLine);
                        }
                    }
                    if(!alreadyPresent){
                        dbSeq.getCvParam().add(mzidLibUtils.makeCvParam("MS:1001088","protein description",psiCV,defLine));
                    }
                }
                
            }
            
        } catch (OutOfMemoryError error) {
          
              String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ error.getMessage()+"."+
                "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
     
        
        }

    }
    
    // Write the new data into a file
    private void writeToMzIdentMLFile() {
        try {
            MzIdentMLMarshaller m = new MzIdentMLMarshaller();
            
            m.marshal(mzIdentML, new FileOutputStream(outputFile));
            

        } catch (FileNotFoundException ex) {
             String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 01 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
        }

    }
    
   
}
