package uk.ac.liv.mzidlib.converters;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import uk.ac.liv.unimod.CompositionT;
import uk.ac.liv.unimod.ModT;
import uk.ac.liv.unimod.ModificationsT;
import uk.ac.liv.unimod.SpecificityT;
import uk.ac.liv.unimod.UnimodT;


/**
 *
 * @author jonesar
 */
public class ReadUnimod {

    private static String inputUnimod =  "unimod.xml";
    private List<ModT> modList;



    public ReadUnimod(){
        try {
            //Use the getResourceAsStream trick to read the unimod.xml file as  
            //a classpath resource. This enables us to also distribute the unimod.xml file
            //inside the .jar file which simplifies usage of the solution as no extra 
            //classpath or path configurations are needed to let the code below find
            //the unimod.xml file: 
        InputStream stream = ClassLoader.getSystemClassLoader().getResourceAsStream(inputUnimod);
        UnimodT unimod = unmarshal(UnimodT.class,stream);

        ModificationsT mods = unimod.getModifications();

        modList = mods.getMod();

        /*
        for (ModT mod : modList) {

            
                Long id = mod.getRecordId();
                String modName  = mod.getTitle();

                CompositionT comp = mod.getDelta();
                double mass = comp.getMonoMass();


            System.out.println(id + " " + modName + " " + mass);

        }
        */
        } catch (JAXBException ex) {
               String methodName =Thread.currentThread().getStackTrace()[1].getMethodName();
             String className = this.getClass().getName();
             String message= "The task \""+methodName +  "\" in the class \""+ className + "\" was not completed because of "+ ex.getMessage()+"."+
                "\nPlease see the reference guide at 04 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
             System.out.println (message);
        }

        
        

    }

    
    public ModT getModByMass(double testMass, double massError, boolean isMono, char res){

    	List<String> residues = new ArrayList<String>();
    	residues.add(""+res);
    	return getModByMass(testMass, massError, isMono, residues);

    }

    /**
     * This method will look for unimod modification entries that are specific for the given
     * residues and have a mass that falls within the massError window around the given testMass.
     * It returns, from the entries found, the one with the smallest mass difference.
     * 
     * @param testMass
     * @param massError
     * @param isMono
     * @param residues
     * @return
     */
    public ModT getModByMass(double testMass, double massError, boolean isMono, List<String> residues){

        ModT foundMod = null;
        boolean isFound = false;
        
   
        	double diffFound = 1000000.0;   //Choose smallest mass difference

            for (ModT mod : modList) {

                CompositionT comp = mod.getDelta();
                double mass;

                if(isMono){
                    mass = comp.getMonoMass();
                }
                else{
                    mass = comp.getAvgeMass();
                }

                boolean siteMatch = false;
                
                
                //check if the modification is a modification that can 
                //occur on the given site/residue
                for(SpecificityT spec : mod.getSpecificity()){
                    
                    for(String residue : residues){
                        
                        if(residue.equals("[")){
                            residue = "N-term";
                        }
                        else if(residue.equals("]")){
                            residue = "C-term";
                        }
                        
                        String site = spec.getSite();
                        if(site.equals(residue)){
                            siteMatch=true;
                            break;
                        }
                    }
                    if (siteMatch)
                    	break;
                }

                if(mass < testMass + massError && mass > testMass - massError && siteMatch){                    
                	//Choose smallest mass difference
                	if(Math.abs(mass - testMass) < diffFound){
                        //System.out.println("Error: Multiple mods found with same mass, choosi: " + testMass);
                        foundMod = mod; 
                        diffFound = Math.abs(mass - testMass);
                    }
                    isFound = true;
                }               
            }
            
            if(!isFound){
                //System.out.println("No mod found in Unimod with mass:" + testMass + " error: " + massError + " residue: " + residues.toString());
            }

        

        return foundMod;

    }


    public <T> T unmarshal( Class<T> docClass, InputStream inputStream )
        throws JAXBException {
        String packageName = docClass.getPackage().getName();
        JAXBContext jc = JAXBContext.newInstance( packageName );
        Unmarshaller u = jc.createUnmarshaller();
        @SuppressWarnings("unchecked")
		JAXBElement<T> doc = (JAXBElement<T>)u.unmarshal( inputStream );
        return doc.getValue();
    }



}
