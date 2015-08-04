package uk.ac.liv.mzidlib.experimental;
    
import java.io.InputStream;
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
public class WriteUnimodToMSGFPlus {

    
    private static String inputUnimod =  "/resources/unimod.xml";
    private List<ModT> modList;


    public static void main(String args[]){
        
        WriteUnimodToMSGFPlus writeConfig = new WriteUnimodToMSGFPlus();
        writeConfig.printAllMods();
        
    }
    

    public WriteUnimodToMSGFPlus(){
        try{
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

        }
        catch(Exception e){
            e.printStackTrace();

        }
        

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
    public void printAllMods(){

        ModT foundMod = null;
        boolean isFound = false;
        
        try{

            
            for (ModT mod : modList) {
                String modLine = "";
                CompositionT comp = mod.getDelta();
                double mass = comp.getMonoMass();

                
          
                String psiMs = mod.getTitle();
                
                //check if the modification is a modification that can 
                //occur on the given site/residue
                
                /*
                 * spec.getPosition().PROTEIN_N_TERM;
                            spec.getPosition().PROTEIN_C_TERM;
                            spec.getPosition().ANY_N_TERM;
                            spec.getPosition().ANY_C_TERM;
                 */
                for(SpecificityT spec : mod.getSpecificity()){
                        String site = spec.getSite();
                        String position = spec.getPosition().value();
                        
                        if(site.equals("N-term") || site.equals("C-term")){
                            site = "*";
                        }
                        
                        if(position.equals("Anywhere")){
                            position = "any";
                        }
                        else if(position.equals("Protein N-term")){
                            position = "Prot-N-term";                            
                        }
                        else if(position.equals("Protein C-term")){
                            position = "Prot-N-term";                            
                        }
                        else if(position.equals("Any N-term")){
                            position = "N-term";                            
                        }
                        else if(position.equals("Any C-term")){
                            position = "C-term";                            
                        }
                        else{
                            System.out.println("Position not recognized:" + position);
                        }
                        
                        if(spec.getClassification().equals("Post-translational") || spec.getClassification().equals("multiple")||spec.getClassification().equals("Artefact")){
                        
                            modLine=""+mass+","+site+",opt,"+position+","+psiMs;
                            System.out.println(modLine);
                        }
                }
            
            }
            


        }
        catch(Exception e){

            e.printStackTrace();
        }



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
