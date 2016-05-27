package uk.ac.liv.mzidlib;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;
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
 */
public class AddEmpaiToMzid {

    private Map<String, String> accToSeq = new HashMap<>();
    private Map<String, String> accToDefline = new HashMap<>();
    private Map<String, String[]> accToPeptides = new HashMap<>();
    private Map<String, Integer> accToObservable = new HashMap<>();
    private Map<String, DBSequence> idToDBSeq = new HashMap<>();
    private Map<String, SpectrumIdentificationItem> idToSIIPassThresholdMap = new HashMap<>();

    private String inputMzid = "example_files/Toxo-1D_OF.mzid";
    private String inputFasta = "example_files/TgondiiME49_ToxoDB-6_2.fasta";   //default for testing
    //private String inputFasta = "example_files/SerumAlbumin.fasta";
    private String outputFile = "edited_file_from_fasta.mzid";   //default for testing
    private String accessionRegex = " ";                        //Separate accessions from defline by a space by default
    private String enzymeRegex = "(?<=[KR])(?!P)";                        //Separate accessions from defline by a space by default
    private int missedCleavages = 1;

    private final Double HMASS = 1.007825;

    double lowMr = 1000000.0;
    double highMr = 0.0;
    private int lowCharge = 1;
    private int highCharge = 3;

    private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;

    private Cv psiCV;
    private Cv unitCV;

    private boolean verbose = false;

    private Map<String, Double> aaMap = new HashMap<>();

    private MzidLibUtils mzidLibUtils;

    /*
     * Main method for testing only
     * 
     */
    public static void main(String args[]) {
        AddEmpaiToMzid addEMPAI = new AddEmpaiToMzid();
    }

    private AddEmpaiToMzid() {
        init();
    }

    public AddEmpaiToMzid(String mzidIn, String mzidOut, String fastaIn, String accRegex, String enzRegex, boolean verboseMode) {
        this.inputMzid = mzidIn;
        this.outputFile = mzidOut;
        this.inputFasta = fastaIn;
        this.accessionRegex = accRegex;
        this.verbose = verboseMode;

        if (enzRegex != null && !enzRegex.equals("")) { //defaults to trypsin
            this.enzymeRegex = enzRegex;
        }

        //this.missedCleavages = missCleaves;
        System.out.println("Assuming charge range: " + lowCharge + " to " + highCharge);
        System.out.println("Accession regex /" + accRegex + "/");
        System.out.println("Verbose mode:" + verbose);
        init();
    }

    private void init() {

        mzidLibUtils = new MzidLibUtils();

        fillAAMap();
        //getPepMZ("PEPTIDE",1);

        System.out.println("Preprocess file\n");
        preProcessMZIDFile();
        System.out.println("Read FASTA\n");
        readFasta();
        System.out.println("Get observables\n");
        getObservableCount();
        System.out.println("Insert emPAI\n");
        insertEmPAIIntoMzid();
        System.out.println("Write MZID\n");
        writeToMzIdentMLFile();

        //insertDetailsIntoMzid();
        //writeToMzIdentMLFile();
    }

    private void readFasta() {
        BufferedReader in = null;
        try {
            in = new BufferedReader(new FileReader(inputFasta));
            /*
             while ((line = in.readLine()) != null) {
             lineNum++;
             string = string.trim();
             if (lineNum==refLine) {
             break;
             }
             }
             */
            //FileInputStream fstream = new FileInputStream(inputFasta);
            // Get the object of DataInputStream
            //DataInputStream in = new DataInputStream(fstream);
            //BufferedReader br = new BufferedReader(new InputStreamReader(in));
            String line;
            String currSequence = "";
            String currProtAcc = null;
            String currDefline = null;
            int recordCounter = 0;
            while ((line = in.readLine()) != null) {
                line = line.replaceAll("\n", "");
                line = line.replaceAll("\r", "");

                if (line.contains(">")) {
                    //Insert previous into hash and reset
                    if (recordCounter != 0) {
                        currSequence = currSequence.replaceAll(" ", "");
                        accToPeptides.put(currProtAcc, currSequence.split(enzymeRegex));
                        accToSeq.put(currProtAcc, currSequence);
                        //System.out.println("Inserting:" + currProtAcc + "_" + currSequence);
                        accToDefline.put(currProtAcc, currDefline);
                        //System.out.println("Inserting2:" + currProtAcc + "_" + currDefline);

                        currSequence = "";
                    }

                    line = line.replaceAll(">", "");

                    int splitPos = line.indexOf(accessionRegex);
                    if (splitPos != -1) {
                        currProtAcc = line.substring(0, splitPos);
                        currDefline = line.substring(splitPos + 1);
                    } else {
                        System.out.println("Regular expression not found for split: " + line);

                    }

                    recordCounter++;
                } else {
                    currSequence += line;
                }
            }
            //handle last
            currSequence = currSequence.replaceAll(" ", "");
            accToSeq.put(currProtAcc, currSequence);
            accToDefline.put(currProtAcc, currDefline);
            accToPeptides.put(currProtAcc, currSequence.split(enzymeRegex));
            //System.out.println("Put last:" + currProtAcc + " seq: " + currSequence);
            //Close the input stream
            in.close();
        } catch (FileNotFoundException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 01 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);

        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);

        } finally {
            try {
                in.close();
            } catch (IOException ex) {
                String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
                String className = this.getClass().getName();
                String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                        + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
                System.out.println(message);

            }
        }

    }

    private void getObservableCount() {

        //lowMr=700.0;
        //highMr=2800.0;  -  for testing SerumAlbumin example only
        for (String acc : accToPeptides.keySet()) {

            String[] peptides = accToPeptides.get(acc);

            if (verbose) {
                System.out.println("\nProcessing prot:" + acc + " total peptides:" + peptides.length);
            }

            int observable = 0;
            for (String pep : peptides) {
                double mr = getPepMr(pep);
                if (mr != -999) {
                    if (mr >= lowMr && mr <= highMr) {
                        if (verbose) {
                            System.out.println("\tobserverable pep:" + pep + " mr:" + mr);
                        }
                        observable++;
                    } else {
                        if (verbose) {
                            System.out.println("\tnon observerable pep:" + pep + " mr:" + mr);
                        }
                    }
                } else {

                    if (verbose) {
                        System.out.println("\tError code for pep:" + pep);
                    }
                }

            }
            accToObservable.put(acc, observable);
            if (verbose) {
                System.out.println("obsCount:" + observable);
            }
        }
    }

    private void fillAAMap() {

        aaMap.put("A", 71.037114);
        aaMap.put("R", 156.101111);
        aaMap.put("N", 114.042927);
        aaMap.put("D", 115.026943);
        aaMap.put("C", 103.009185);
        aaMap.put("E", 129.042593);
        aaMap.put("Q", 128.058578);
        aaMap.put("G", 57.021464);
        aaMap.put("H", 137.058912);
        aaMap.put("I", 113.084064);
        aaMap.put("L", 113.084064);
        aaMap.put("K", 128.094963);
        aaMap.put("M", 131.040485);
        aaMap.put("F", 147.068414);
        aaMap.put("P", 97.052764);
        aaMap.put("S", 87.032028);
        aaMap.put("T", 101.047679);
        aaMap.put("U", 150.95363);
        aaMap.put("W", 186.079313);
        aaMap.put("Y", 163.06332);
        aaMap.put("V", 99.068414);
        aaMap.put("[", 0.0);
        aaMap.put("]", 0.0);
        /*        
         A=71.037114
         B=114.534940
         C=160.030649
         D=115.026943
         E=129.042593
         F=147.068414
         G=57.021464
         H=137.058912
         I=113.084064
         J=0.000000
         K=128.094963
         L=113.084064
         M=131.040485
         N=114.042927
         O=0.000000
         P=97.052764
         Q=128.058578
         R=156.101111
         S=87.032028
         T=101.047679
         U=150.953630
         V=99.068414
         W=186.079313
         X=111.000000
         Y=163.063329
         Z=128.550590
         Hydrogen=1.007825
         Carbon=12.000000
         Nitrogen=14.003074
         Oxygen=15.994915
         Electron=0.000549
         C_term=17.002740
         N_term=1.007825
         */

    }

    private void addFixedModsToAAMap(String aa, double modMass) {

        System.out.println("aa " + aa);
        double mass = aaMap.get(aa);
        double newMass = mass + modMass;
        aaMap.put(aa, newMass);
        System.out.println("Changed mass of " + aa + " to " + newMass);
    }

    private double getPepMr(String peptide) {
        double mr = 0.0;
        //double hmass = 1.007825;

        peptide = peptide.toUpperCase();
        peptide = peptide.trim();
        peptide = "[" + peptide + "]";
        //System.out.print("Peptide:" + peptide);

        String[] aas = peptide.split("");
        //double pepMass = 0.0;
        for (String aa : aas) {
            if (aa.length() == 1) {
                if (aaMap.get(aa) != null) {
                    mr += aaMap.get(aa);
                } else {
                    mr = -999;
                    break;
                }
            }
        }

        //System.out.println( " mr: " + mr); 
        return mr;

    }

    /*
     * This does first pass through MZID file to load SIIs into a hash, find low and high Mr values and alter aaMap with fixedMods
     */
    private List<Double> preProcessMZIDFile() {
        List<Double> lowAndHighMr = new ArrayList<>();

        int psmCount = 0;

        mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));
        mzIdentML = (MzIdentML) mzIdentMLUnmarshaller.unmarshal(MzIdentML.class);
        Iterator<Cv> iterCv = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
        while (iterCv.hasNext()) {
            Cv cv = iterCv.next();
            if (cv.getUri().toLowerCase().contains("psi")) {
                psiCV = cv;
            } else if (cv.getUri().toLowerCase().contains("unit")) {
                unitCV = cv;
            }
        }

        Iterator<SearchModification> iterSM = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SearchModification);
        //System.out.println("iterSM" + iterSM.toString() + " test: " + iterSM.hasNext());

        while (iterSM.hasNext()) {
            SearchModification searchMod = iterSM.next();
            //System.out.println("Search mod: " + searchMod.getMassDelta());
            if (searchMod.isFixedMod()) {
                float modMass = searchMod.getMassDelta();
                for (String res : searchMod.getResidues()) {

                    if (searchModIsNTerminal(searchMod)) {
                        addFixedModsToAAMap("[", 1.0 * modMass);
                    } else if (searchModIsCTerminal(searchMod)) {
                        addFixedModsToAAMap("]", 1.0 * modMass);
                    } else {
                        addFixedModsToAAMap(res, 1.0 * modMass);
                    }
                }
            }
        }

        Iterator<DBSequence> iterDBSeq = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        //System.out.println("iterSM" + iterSM.toString() + " test: " + iterSM.hasNext());

        while (iterDBSeq.hasNext()) {
            DBSequence dbSeq = iterDBSeq.next();
            idToDBSeq.put(dbSeq.getId(), dbSeq);
        }

        Iterator<SpectrumIdentificationItem> iterSII = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationItem);

        while (iterSII.hasNext()) {
            SpectrumIdentificationItem sii = iterSII.next();

            if (sii.isPassThreshold()) {
                idToSIIPassThresholdMap.put(sii.getId(), sii);
                psmCount++;
                double mz = sii.getExperimentalMassToCharge();
                int charge = sii.getChargeState();
                double mr = (sii.getExperimentalMassToCharge() * charge) - (HMASS * charge);
                if (mr < lowMr) {
                    lowMr = mr;
                } else if (mr > highMr) {
                    highMr = mr;
                }
            }
        }

        if (psmCount == 0) {

            System.out.println("Error: no PSMs with passThreshold = true, exiting");

        }

        System.out.println("Applying limits to emPAI calculation using mrs:" + lowMr + " to " + highMr);

        lowAndHighMr.add(lowMr);
        lowAndHighMr.add(highMr);
        return lowAndHighMr;
    }

    /*
     * Helper method to check if a search mod is a terminal mod
     */
    private boolean searchModIsNTerminal(SearchModification mod) {
        boolean isNTerminal = false;

        List<String> terminalModCvTerms = new ArrayList<>();
        terminalModCvTerms.add("MS:1001189");
        for (CvParam param : mod.getCvParam()) {
            if (terminalModCvTerms.contains(param.getAccession())) {
                isNTerminal = true;
            }
        }

        return isNTerminal;
    }

    /*
     * Helper method to check if a search mod is a terminal mod
     */
    private boolean searchModIsCTerminal(SearchModification mod) {
        boolean isCTerminal = false;

        List<String> terminalModCvTerms = new ArrayList<>();
        terminalModCvTerms.add("MS:1001190");
        for (CvParam param : mod.getCvParam()) {
            if (terminalModCvTerms.contains(param.getAccession())) {
                isCTerminal = true;
            }
        }

        return isCTerminal;
    }

    private void insertEmPAIIntoMzid() {

        Iterator<ProteinDetectionHypothesis> iterPDH = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinDetectionHypothesis);

        for (ProteinAmbiguityGroup pag : mzIdentML.getDataCollection().getAnalysisData().getProteinDetectionList().getProteinAmbiguityGroup()) {

            for (ProteinDetectionHypothesis pdh : pag.getProteinDetectionHypothesis()) {

                Map<Double, String> uniqueParentIonsMap = new HashMap<>();

                if (verbose) {

                    System.out.println("Calculating emPAI for " + pdh.getDBSequenceRef());
                }
                for (PeptideHypothesis ph : pdh.getPeptideHypothesis()) {
                    for (SpectrumIdentificationItemRef siiRef : ph.getSpectrumIdentificationItemRef()) {
                        String siiId = siiRef.getSpectrumIdentificationItemRef();
                        SpectrumIdentificationItem sii = idToSIIPassThresholdMap.get(siiId);
                        if (sii != null && sii.isPassThreshold()) {
                            int charge = sii.getChargeState();
                            //Double pepMr = (sii.getCalculatedMassToCharge() * charge) - (charge * HMASS);
                            uniqueParentIonsMap.put(sii.getCalculatedMassToCharge(), sii.getId());
                        }
                    }
                }
                int observed = uniqueParentIonsMap.size();

                if (verbose) {
                    for (Double unique : uniqueParentIonsMap.keySet()) {
                        System.out.println("\tUnique value: " + unique + " from (for example): " + uniqueParentIonsMap.get(unique));
                    }
                }

                String protAcc = idToDBSeq.get(pdh.getDBSequenceRef()).getAccession();
                if (protAcc != null) {
                    if (accToObservable.get(protAcc) != null) {
                        int observable = accToObservable.get(protAcc);
                        Double empai = Math.pow(10, (0.0 + observed) / (0.0 + observable)) - 1;

                        if (verbose) {
                            System.out.println(protAcc + "\t" + observed + "\t" + observable + "\t" + empai);
                        }
                        pdh.getCvParam().add(mzidLibUtils.makeCvParam("MS:1001905", "emPAI value", psiCV, "" + empai));

                        if (verbose) {
                            pdh.getUserParam().add(mzidLibUtils.makeUserParam("observed", "" + observed));
                            pdh.getUserParam().add(mzidLibUtils.makeUserParam("observable", "" + observable));
                        }
                    } else {
                        System.out.println("Protein: " + protAcc + " not found in fasta:" + inputFasta);
                    }
                }
            }

        }
    }

    /*
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
     Boolean alreadyPresent = false;
     for(CvParam cvParam : dbSeq.getCvParam()){
     if(cvParam.getAccession().equals("MS:1001088")){
     alreadyPresent = true;
     cvParam.setValue(defLine);
     }
     }
     if(!alreadyPresent){
     dbSeq.getCvParam().add(makeCvParam("MS:1001088","protein description",psiCV,defLine));
     }
     }
                
     }
            
     } catch (OutOfMemoryError error) {
     System.out.println("Out of Memory Error: " + error.getMessage());
           
     } catch (Exception ex) {
     ex.printStackTrace();
     }

     }
     */
    // Write the new data into a file
    private void writeToMzIdentMLFile() {
        try {
            MzIdentMLMarshaller m = new MzIdentMLMarshaller();

            m.marshal(mzIdentML, new FileOutputStream(outputFile));
        } catch (FileNotFoundException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 01 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

}
