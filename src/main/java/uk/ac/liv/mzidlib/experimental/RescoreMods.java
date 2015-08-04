package uk.ac.liv.mzidlib.experimental;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 *
 * @author jonesar
 *
 * Tester class to rescore PTMs, based on co-occurrence with unmodifed pair and
 * rarity of AA
 *
 * Notes:
 *
 * Need to get MSGF+ to export top 10 hits per PSM For every modified PTM,
 * generate a new mod liklihood score based on if peptide pair is present,
 * re-score by dividing e-value by 100 say? Use AA frequencies from PTMCuration
 *
 *
 * good example here:
 *
 * C:\Work\MassSpec\iprg_2012_archive\MSGF\..\iPRG2012.mgf	index	16374	16375
 * sample=1 period=1 cycle=3116 experiment=5	-1	SII_16375_1	1	TRUE	985.1687622
 * 985.1687622	3	DGKAPEAEQVLSAAATFPIAQPATDVEAR	DGKAPEAEQVJSAAATFPJAQPATDVEAR
 * DGKAPEAEQVJSAAATFPJAQPATDVEAR	1	110	147	1.57E-18	5.54E-13	0	0
 * sp|P40212|RL13B_YEAST_129_157_R_A	FALSE
 * C:\Work\MassSpec\iprg_2012_archive\MSGF\..\iPRG2012.mgf	index	16374	16375
 * sample=1 period=1 cycle=3116 experiment=5	-1	SII_16375_2	2	TRUE	985.1687622
 * 985.1687622	3	NGKAPEAEQVLSAAATFPIAQPATDVEAR	NGKAPEAEQVJSAAATFPJAQPATDVEAR
 * Deamidated:1	DGKAPEAEQVJSAAATFPJAQPATDVEAR	0	110	147	1.57E-18	5.54E-13	0	0
 * sp|Q12690|RL13A_YEAST_129_157_R_A	FALSE
 *
 * Mod and normal score the same, therefore PTMScoring must down-weight all PTMs
 * compared to normal hits...
 *
 *
 *
 */
public class RescoreMods {

    //private String inputMzid = "resources/iprg2012_msgf_restrict_threshold.mzid";
    //private String outputMzid = "test_outputs/iprg2012_msgf_restrict_rescored_threshold_03April2013.mzid";   //default for testing
    private String inputMzid = "C:/Work/MassSpec/iprg_2012_archive/MSGF/iprg2012_msgf_restrict_nd.mzid";
    private String outputMzid = "C:/Work/MassSpec/iprg_2012_archive/MSGF/iprg2012_msgf_restrict_nd_rescored_threshold_03April2013.mzid";   //default for testing
    //private String inputIPRG2012Answer = "resources/iprg_answer_key.txt";
    //private String inputIPRG2012AnswerSpiked = "resources/iprg_answer_key_spiked.txt";
    private String inputIPRG2012Answer = "/resources/iprg_answer_key_nd.txt";

    private String inputIPRG2012AnswerSpiked = "/resources/iprg_answer_key_spiked_nd.txt";
    private static String inputUnimod =  "/resources/unimod.xml";
    private Map<String, String> iprgAnswerMap = new HashMap<>();
    private Map<String, String> iprgAnswerSpikedMap = new HashMap<>();
    private Map<String, Peptide> peptideIDMap = new HashMap<>();
    private Map<String, Peptide> modPeptideIDMap = new HashMap<>();
    private Map<String, SpectrumIdentificationItem> siiIDMap = new HashMap<>();
    private Map<String, List<SpectrumIdentificationItem>> unmodifiedPepSeqToSIIMap = new HashMap<>();
    private Map<String, List<SpectrumIdentificationItem>> modifiedPepSeqToSIIMap = new HashMap<>();
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private MzidLibUtils utils = new MzidLibUtils();        //Utils class containing makeCvParam methods
    //private int[] commonUnimodIDs = new int[]{35,27,28,385}; //Ox on M, pyro-glu E/Q, ammonia loss on N-terminal C
    //private int[] mediumCommonBioModIDs = new int[]{21,}; 
    private Map<String, String> commonUnimodIDs = new HashMap<>();
    private Map<String, String> mediumCommonUnimodIDs = new HashMap<>();
    //private Double commonModWeight = 0.2;   //Now set all to be same weight, as per PeaksPTM
    //private Double mediumModWeight = 0.5;
    //private Double rareModWeight = 1.0;
    private double commonModWeight = 1.0;   //Now set all to be same weight, as per PeaksPTM - we're not searching for any genuinely rare PTMs, these cancel out the general Mod weight
    private double mediumModWeight = 1.0;
    private double rareModWeight = 1.0;
    private double generalModWeight = 0.1;          //Slightly down-weight var mod over unmod alternative
    private double multipleVarModWeight = 0.01;      //Don't downscore multiple mods
    private double pairedModAndUnmodWeight = 5.0;
    //Thus a common mod without a pair would be downweighted 0.5 (0.1 *5). 
    //A rare mod with a pair would be up-weighted 10 fold, a medium fold 20 fold
    //These numbers are just guesses
    private String cvAccForBaseScore = "MS:1002053";        //For testing we will use MSGF:evalue
    //private Boolean logTransformScore = true;
    private String newCvAcc = "MS:100PTMSCORE";                 //To be replaced with real CV term
    private String newCvName = "PTM score";
    private Map<String, DBSequence> dbSequenceHashMap = new HashMap<>();
    private Map<String, Peptide> peptideHashMap = new HashMap<>();
    private Map<String, PeptideEvidence> peptideEvidenceHashMap = new HashMap<>();
    private Cv psiCV;

    public static void main(String args[]) {

        RescoreMods rescoreMods = new RescoreMods();
        rescoreMods.init();
    }

    public RescoreMods(){
        
    }
    
    public RescoreMods(String inputFileName, String outputFileName, String cvAccForScoreToAdapt,  double commonModificationWeight, double mediumModificationWeight, double rareModificationWeight, double generalModificationWeight, double pairedModificationAndUnmodWeight, double multipleVariableModWeight){
        
        inputMzid = inputFileName;
        outputMzid = outputFileName;
        
        commonModWeight = commonModificationWeight;
        mediumModWeight = mediumModificationWeight;
        rareModWeight = rareModificationWeight;
        generalModWeight = generalModificationWeight;        
        pairedModAndUnmodWeight = pairedModificationAndUnmodWeight;
        cvAccForBaseScore = cvAccForScoreToAdapt;
        //logTransformScore = logTransformPSMScore;
        multipleVarModWeight = multipleVariableModWeight;
        
        System.out.println("We will score with accession ID: " + cvAccForScoreToAdapt);
        System.out.println("Code assumes we are working off an e-value type score, ordered low to high, log transform has not been implemented");
        System.out.println(newCvName + " for all PSMs with common mods will be multiplied by " + commonModWeight);
        System.out.println(newCvName + " for all PSMs with medium common mods will be multiplied by " + mediumModWeight);
        System.out.println(newCvName + " for all PSMs with rare mods will be multiplied by " + rareModWeight);
        System.out.println(newCvName + " for all PSMs with any type of variable mod will be multiplied by " + generalModWeight);
        System.out.println(newCvName + " for all PSMs with with a mod for which unmodified paired peptides has also been identified will be multiplied by " + pairedModAndUnmodWeight);
        
        init();
    }
    
    private void init() {

        System.out.println("Note: code assumes threshold has been run to set passThreshold=true sensibly e.g. FDR < 0.01. "
                + "This is needed so that we can profile high-quality unmodifed PSMs");

        populateModMaps();
        
        System.out.println("currently hard-coded to read IPRG_2012_ND answer key and insert into results ");
        readIPRGAnswers(inputIPRG2012Answer, inputIPRG2012AnswerSpiked);
        readPSMs(inputMzid);
        writeToMzIdentMLFile(outputMzid);
    }

    private void readIPRGAnswers(String inputAnswerFile, String inputAnswerSpikedFile) {

        try {

            InputStream stream = ClassLoader.getSystemClassLoader().getResourceAsStream(inputAnswerFile);
            
            // Get the object of DataInputStream
            DataInputStream in = new DataInputStream(stream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String line;

            while ((line = br.readLine()) != null) {
                line = line.replaceAll("\n", "");
                String[] temp = line.split("\t");
                iprgAnswerMap.put(temp[0], temp[1]);
            }
            br.close();

            InputStream stream2 = ClassLoader.getSystemClassLoader().getResourceAsStream(inputAnswerSpikedFile);
            // Get the object of DataInputStream
            in = new DataInputStream(stream2);
            br = new BufferedReader(new InputStreamReader(in));
            line = null;

            while ((line = br.readLine()) != null) {
                line = line.replaceAll("\n", "");
                String[] temp = line.split("\t");

                iprgAnswerSpikedMap.put(temp[0], temp[1]);
            }
            br.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private void populateModMaps() {

        commonUnimodIDs.put("UNIMOD:35@M", "Anywhere");
        commonUnimodIDs.put("UNIMOD:27@E", "Pep N-term");
        commonUnimodIDs.put("UNIMOD:28@Q", "Pep N-term");
        commonUnimodIDs.put("UNIMOD:7@Q", "Anywhere");
        commonUnimodIDs.put("UNIMOD:7@N", "Anywhere");
        commonUnimodIDs.put("UNIMOD:385@C", "Pep N-term");
        
        System.out.println("Assuming common variable mods:");
        for(String unimodID: commonUnimodIDs.keySet()){
            System.out.println(unimodID + " " + commonUnimodIDs.get(unimodID));
        }
                

        mediumCommonUnimodIDs.put("UNIMOD:21@S", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:21@T", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:21@Y", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:1@K", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:1@A", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:1@M", "Anywhere");
        mediumCommonUnimodIDs.put("UNIMOD:1@S", "Anywhere");

        System.out.println("Assuming \"medium\" mods:");
        for(String unimodID: mediumCommonUnimodIDs.keySet()){
            System.out.println(unimodID + " " + mediumCommonUnimodIDs.get(unimodID));
        }
    }

    private void readPSMs(String inputMzid) {

        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));

        Iterator<SearchModification> iterSearchMod = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SearchModification);

        List<SearchModification> allSearchMods = new ArrayList<>();
        while (iterSearchMod.hasNext()) {
            SearchModification searchMod = iterSearchMod.next();
            allSearchMods.add(searchMod);
        }



        Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);

        int countVarMods = 0;
        int countFixMods = 0;

        while (iterPeptide.hasNext()) {
            Peptide pep = iterPeptide.next();
            peptideIDMap.put(pep.getId(), pep);

            boolean hasVarMods = false;



            if (pep.getModification() != null && !pep.getModification().isEmpty()) {
                for (Modification mod : pep.getModification()) {
                    if (!testIfFoundModWasFixed(allSearchMods, pep, mod)) {
                        modPeptideIDMap.put(pep.getId(), pep);
                        countVarMods++;
                    } else {
                        countFixMods++;
                    }
                }
            }
        }

        System.out.println("Var mods:" + countVarMods + " fixed mods: " + countFixMods);



        Iterator<SpectrumIdentificationItem> iterSII = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);

        while (iterSII.hasNext()) {
            SpectrumIdentificationItem sii = iterSII.next();
            siiIDMap.put(sii.getId(), sii);

            Peptide pep = peptideIDMap.get(sii.getPeptideRef());
            String pepSeq = pep.getPeptideSequence();
            if (modPeptideIDMap.containsKey(sii.getPeptideRef())) {
                List<SpectrumIdentificationItem> siiList = null;
                if (modifiedPepSeqToSIIMap.containsKey(pepSeq)) {
                    siiList = modifiedPepSeqToSIIMap.get(pepSeq);
                } else {
                    siiList = new ArrayList<>();
                }
                siiList.add(sii);
                modifiedPepSeqToSIIMap.put(pepSeq, siiList);
            } else {

                if (sii.isPassThreshold()) {
                    List<SpectrumIdentificationItem> siiList = null;
                    if (unmodifiedPepSeqToSIIMap.containsKey(pepSeq)) {
                        siiList = unmodifiedPepSeqToSIIMap.get(pepSeq);
                    } else {
                        siiList = new ArrayList<>();
                    }
                    siiList.add(sii);
                    unmodifiedPepSeqToSIIMap.put(pepSeq, siiList);
                }
                //i.e. we are discarding unmodified peptides with FDR < 0.01
            }

        }


        for (String id : siiIDMap.keySet()) {
            SpectrumIdentificationItem sii = siiIDMap.get(id);
            Peptide pep = peptideIDMap.get(sii.getPeptideRef());
            String pepSeq = pep.getPeptideSequence();

            boolean scoreFound = false;
            for (CvParam cvParam : sii.getCvParam()) {
                if (cvParam.getAccession().equals(cvAccForBaseScore)) {
                    scoreFound = true;
                    double evalue = Double.parseDouble(cvParam.getValue());
                    double newEvalue = evalue;

                    if (modPeptideIDMap.containsKey(sii.getPeptideRef())) {       //This is a modified peptide
                        newEvalue *= generalModWeight;
                        if (unmodifiedPepSeqToSIIMap.containsKey(pepSeq)) {       //A high-quality unmodified version exists
                            newEvalue *= pairedModAndUnmodWeight;
                            sii.getCvParam().add(utils.makeCvParam("MS:1VXYZ", "Paired Mod classification", psiCV, "Paired mod"));
                        }

                        Modification varMod = null;
                        int varModCounter = 0;
                        for (Modification mod : pep.getModification()) {
                            if (!testIfFoundModWasFixed(allSearchMods, pep, mod)) {
                                varMod = mod;
                                varModCounter++;
                            }
                        }

                        if (varModCounter > 1) {                                  //This peptide has 2 or more var mods
                            newEvalue *= multipleVarModWeight;
                            sii.getCvParam().add(utils.makeCvParam("MS:1XYZ", "Modification classification", psiCV, "Multiple"));
                        } else {
                            if (varMod != null) {                                   //This peptide has one var mod
                                
                                String modRes = "";
                                if(varMod.getLocation()>0){
                                    modRes = ""+pepSeq.charAt(varMod.getLocation()-1);
                                }
                                else{
                                    modRes = ""+pepSeq.charAt(0);
                                }
                                
                                boolean isCommonMod = false;
                                boolean isMediumMod = false;

                                for (CvParam param : varMod.getCvParam()) {                                  
                                    
                                    String modKey = param.getAccession() + "@" + modRes;
                                    if (commonUnimodIDs.get(modKey) != null) {    //present in common mods                                        
                                        String res = commonUnimodIDs.get(param.getAccession());
                                        isCommonMod = true;
                                    } else if (mediumCommonUnimodIDs.get(modKey) != null) {
                                        isMediumMod = true; //TODO not implemented medium N-terminal mod code                                            
                                    }
                                }

                                if (isCommonMod) {
                                    newEvalue *= commonModWeight;
                                    sii.getCvParam().add(utils.makeCvParam("MS:1XYZ", "Modification classification", psiCV, "Common"));
                                } else if (isMediumMod) {
                                    newEvalue *= mediumModWeight;
                                    sii.getCvParam().add(utils.makeCvParam("MS:1XYZ", "Modification classification", psiCV, "Medium"));
                                } else {
                                    newEvalue *= rareModWeight;
                                    sii.getCvParam().add(utils.makeCvParam("MS:1XYZ", "Modification classification", psiCV, "Rare"));
                                }

                            }
                        }
                    }
                    sii.getCvParam().add(utils.makeCvParam(newCvAcc, newCvName, psiCV, "" + newEvalue));  //Add new score to SII                        
                }
            }


            if (!scoreFound) {
                System.out.println("Error - SII without required score found, accession:" + cvAccForBaseScore + " on " + sii.getId());
            }
        }
    }


    /*
     * Helper class to determine if a modification that has been identified was
     * a fixed modification Assumes Unimod IDs have been used for mods
     */
    private boolean testIfFoundModWasFixed(List<SearchModification> searchMods, Peptide pep, Modification mod) {
        boolean isFixed = false;

        SearchModification matchedMod = null;

        String modCvParamAcc = null;
        for (CvParam cvParam : mod.getCvParam()) {
            if (cvParam.getAccession().contains("UNIMOD:")) {
                modCvParamAcc = cvParam.getAccession();
            }
        }

        String pepSeq = pep.getPeptideSequence();
        String modRes = "";
        if (mod.getLocation() == 0) {
            modRes = "." + pepSeq.substring(0, 1);       //Any char plus first residue to cope with both types
        } else if (mod.getLocation() == pepSeq.length() + 1) {
            modRes = "." + pepSeq.substring(mod.getLocation() - 1, mod.getLocation()); //Any char plus last residue
        } else {
            modRes = pepSeq.substring(mod.getLocation() - 1, mod.getLocation());
        }
        //System.out.print("r:" + modRes);

        for (SearchModification searchMod : searchMods) {
            if (matchedMod == null) {
                String searchModCvParamAcc = null;
                for (CvParam cvParam : searchMod.getCvParam()) {
                    if (cvParam.getAccession().contains("UNIMOD:")) {
                        searchModCvParamAcc = cvParam.getAccession();
                    }
                }

                if (modCvParamAcc != null && searchModCvParamAcc != null) {
                    //System.out.println("\tHere");
                    if (modCvParamAcc.equals(searchModCvParamAcc)) {
                        for (String searchModRes : searchMod.getResidues()) {
                            if (modRes.contains(searchModRes)) {
                                matchedMod = searchMod;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (matchedMod != null) {
            isFixed = matchedMod.isFixedMod();
        } else {
            System.out.println("No matching search mod found for:" + mod.getCvParam().get(0).getAccession() + "@" + mod.getLocation());
        }

        return isFixed;
    }

    // Write the new data into a file
    private void writeToMzIdentMLFile(String outputFile) {
        try {
            String outFile = outputFile;
//            if (!outFile.endsWith(".mzid")) {
//                outFile = outFile + ".mzid";
//            }
            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            //cvList = mzIdentML.getCvList();
            CvList cvList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
            Iterator<Cv> iterCv = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                }
            }
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            AnalysisSoftwareList analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            AuditCollection auditCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            Provider provider = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            AnalysisProtocolCollection analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            AnalysisCollection analysisCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);

            writer.write(marshaller.createXmlHeader() + "\n");

            String mzID = mzIdentMLUnmarshaller.getMzIdentMLId();
            if (mzID != null) {
                writer.write(marshaller.createMzIdentMLStartTag(mzID) + "\n");
            } else {
                writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            }



            if (cvList != null) {
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");
            if (analysisSoftwareList != null) {
                AnalysisSoftware analysisSoftware = new AnalysisSoftware();
                analysisSoftware.setName(this.getClass().getSimpleName());
            Date date = new Date() ;
            SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss") ;
            analysisSoftware.setName(this.getClass().getSimpleName()+"_"+dateFormat.format(date)); analysisSoftware.setId(this.getClass().getSimpleName()+"_"+dateFormat.format(date));
                analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);

                marshaller.marshal(analysisSoftwareList, writer);
            }
            writer.write("\n");

            if (provider != null) {
                marshaller.marshal(provider, writer);
            }
            writer.write("\n");

            if (auditCollection != null) {
                marshaller.marshal(auditCollection, writer);
            }
            writer.write("\n");

            SequenceCollection sequenceCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);

            marshaller.marshal(sequenceCollection, writer);

            writer.write("\n");


            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }


            Iterator<SpectrumIdentificationResult> sirIter = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (sirIter.hasNext()) {

                SpectrumIdentificationResult sr = sirIter.next();
                List<SpectrumIdentificationItem> siiList = sr.getSpectrumIdentificationItem();

                String specTitle = null;
                for (CvParam param : sr.getCvParam()) {
                    if (param.getAccession().equals("MS:1000796")) {
                        specTitle = param.getValue();
                    }
                }

                //"sample=1 period=1 cycle=1657 experiment=4" to 1657.4
                String newSpecTitle = null;
                String iprgInsert = null;
                String iprgSpiked = "not spiked";
                if (specTitle != null) {
                    String[] temp = specTitle.split(" ");
                    try {
                        //String cycleValue = temp[2].split("=")[1];
                        String cycleValue = temp[5].split("=")[1];
                        cycleValue = cycleValue.substring(0, cycleValue.length() - 1);
                        //String expValue = temp[3].split("=")[1];
                        String expValue = temp[6].split("=")[1];
                        expValue = expValue.substring(0, expValue.length() - 1);      //for _nd file, spectrum title structure has changed
                        newSpecTitle = cycleValue + "." + expValue;

                        //System.out.println(newSpecTitle);

                        if (iprgAnswerMap.containsKey(newSpecTitle)) {
                            iprgInsert = iprgAnswerMap.get(newSpecTitle);
                        } else {
                            iprgInsert = "no_spec_match";
                        }

                        if (iprgAnswerSpikedMap.containsKey(newSpecTitle)) {
                            iprgSpiked = iprgAnswerSpikedMap.get(newSpecTitle);
                        }
                    } catch (Exception e) {
                        System.out.println("e:" + e.getMessage());
                        //17219: Scan 72138 (rt=3903.62, p=0, c=3435, e=2)
                        //sample=1 period=1 cycle=3276 experiment=3
                    }
                }

                List<SpectrumIdentificationItem> newSiiList = new ArrayList<>();
                for (SpectrumIdentificationItem sii : siiList) {
                    SpectrumIdentificationItem newSII = siiIDMap.get(sii.getId());
                    String pepSeq = peptideIDMap.get(sii.getPeptideRef()).getPeptideSequence();
                    String pepSeqIL = pepSeq.replace('I', 'L');

                    if (iprgInsert != null) {
                        newSII.getCvParam().add(utils.makeCvParam("MS:123456", "iprg2012 answer", psiCV, iprgInsert));
                        if (pepSeqIL.equals(iprgInsert)) {
                            newSII.getCvParam().add(utils.makeCvParam("MS:IPRG_MATCH", "IPRG_MATCH", psiCV, "1"));
                        } else {
                            newSII.getCvParam().add(utils.makeCvParam("MS:IPRG_MATCH", "IPRG_MATCH", psiCV, "0"));
                        }

                    }
                    if (iprgSpiked != null) {
                        newSII.getCvParam().add(utils.makeCvParam("MS:123457", "iprg2012 spiked", psiCV, iprgSpiked));
                        if (pepSeqIL.equals(iprgSpiked)) {
                            newSII.getCvParam().add(utils.makeCvParam("MS:IPRG_SPIKE", "IPRG_SPIKE", psiCV, "1"));
                        } else {
                            newSII.getCvParam().add(utils.makeCvParam("MS:IPRG_SPIKE", "IPRG_SPIKE", psiCV, "0"));
                        }
                    }

                    newSiiList.add(newSII);

                }
                sr.getSpectrumIdentificationItem().clear();
                for (SpectrumIdentificationItem sii : newSiiList) {
                    sr.getSpectrumIdentificationItem().add(sii);
                }

                siList.getSpectrumIdentificationResult().add(sr);



            }

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
                marshaller.marshal(analysisProtocolCollection, writer);
            }
            writer.write("\n");


            writer.write(marshaller.createDataCollectionStartTag() + "\n");

            writer.write("\n");

            if (inputs != null) {
                marshaller.marshal(inputs, writer);
            }
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");


            marshaller.marshal(siList, writer);
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
