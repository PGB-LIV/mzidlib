package uk.ac.liv.mzidlib.performance;

import java.io.File;
import java.net.URL;
import java.util.Iterator;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import static uk.ac.ebi.jmzidml.MzIdentMLElement.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author Fawaz Ghali 03-Feb-2015
 */
public class JMzIdentML {

    private static MzIdentMLUnmarshaller unmarshaller;
    private long lStartTime;
    private long lEndTime;
    
    private long difference;
    
    private static String file= "C:\\Users\\fghali\\Desktop\\HV_map_B4_omssa.mzid";


    public JMzIdentML(String file) {
        lStartTime = System.currentTimeMillis();
        System.out.println ("--------------------------------------------------------------");
        unmarshaller = new MzIdentMLUnmarshaller(new File (file));
        lEndTime = System.currentTimeMillis();
       	difference = lEndTime - lStartTime;
	System.out.println("MzIdentMLUnmarshaller init Elapsed milliseconds: " + difference);
        System.out.println ("--------------------------------------------------------------");
        SpectrumIdentificationItem sii;
        Iterator<SpectrumIdentificationItem> siiIter;
        lStartTime = System.currentTimeMillis();
        siiIter  =unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem); 
        lEndTime = System.currentTimeMillis();
       	difference = lEndTime - lStartTime;
	System.out.println("unmarshalCollectionFromXpath SII Elapsed milliseconds: " + difference);
        System.out.println ("--------------------------------------------------------------");
        lStartTime = System.currentTimeMillis();
        sii  =unmarshaller.unmarshal(MzIdentMLElement.SpectrumIdentificationItem); 
        lEndTime = System.currentTimeMillis();
       	difference = lEndTime - lStartTime;
	System.out.println("unmarshal SII Elapsed milliseconds: " + difference);
        System.out.println ("--------------------------------------------------------------");
        
        
 
    }
    
   
    
    public static void main(String args[]) {
        if (args !=null && args.length==1)
        {
            if ( args[0]!=null&& !args[0].equals(""))
                file  = args[0];
            
        }
        
        JMzIdentML jMzIdentML = new JMzIdentML(file);
        

        
    }

}
