package uk.ac.liv.mzidlib.experimental;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.xml.bind.JAXBException;

import org.apache.commons.lang.StringEscapeUtils;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author fghali
 */
public class CountMods {

    private MzIdentMLUnmarshaller unmarshaller;
 
    public CountMods(String mzid, String out) {
        // C:\Users\fghali\Desktop\CPTAC_RT_ALIGNMENTS\mam_042408o_CPTAC_study6_6C008_s2_fdr_threshold.mzid
        unmarshaller = new MzIdentMLUnmarshaller(new File(mzid));
        
        read();
    }

    public static void main(String[] args) {
        CountMods cm = new CountMods(args[0],"");
        
    }
    private void read() {
        List<String> list = new ArrayList<String>();
        try {
            // loop SIR
            Iterator<SpectrumIdentificationResult> iterSIR = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            List<SpectrumIdentificationResult> sirList = new ArrayList<>();
            while (iterSIR.hasNext()) {
                SpectrumIdentificationResult sir = iterSIR.next();
                sirList.add(sir);
                List<SpectrumIdentificationItem> listSII = sir.getSpectrumIdentificationItem();
                // loop SII
                for (SpectrumIdentificationItem sii : listSII) {
                   String peptideRef = sii.getPeptideRef();
                   Peptide peptide = unmarshaller.unmarshal(Peptide.class, StringEscapeUtils.escapeXml(peptideRef));
                   if(peptide!=null && sii.isPassThreshold()){
                       int peptideLength = peptide.getPeptideSequence().length();
                       // loop modifications
                       List<Modification> modificationList = peptide.getModification();
                       for (int i = 0; i < modificationList.size(); i++) {
                           Modification modification = modificationList.get(i);
                           int location = modification.getLocation().intValue();
                           String unimod = modification.getCvParam().get(0).getName();
                           // If the mod is at position 0, then it is an N-terminal mod 
                           // If the mod is at position peptide length + 1, it is an N-terminal mod.
                           if ((location==0) || (location == peptideLength + 1))
                           {
                             list.add(unimod + " N="); 
                           }
                           else{
                               // the mod is at any other position, it is simply a standard variable modification
                             list.add(unimod + " ="); 
                           }
                       }
                   }
                }
            }
                       
            // print out unique values and their counts.
            Set<String> uSet = new HashSet<String>(list);
            for (String entry : uSet) {
                System.out.println(entry + " " + Collections.frequency(list, entry));
            }


        } catch (JAXBException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 04 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }
}
