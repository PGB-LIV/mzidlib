/**
 * @author Fawaz Ghali, Ritesh Krishna, University of Liverpool 2012
 */
package uk.ac.liv.mzidlib.multiplesearch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

public class MzIdentMLReader {

    String searchEngineIdentifierScore;
    String searchEngineIdentifierExpectation; // The value used in the original algorithm
    String inputMzIdentMLFile;
    String searchEngine;
    String mascotDecoyTag;
    Map<String, List<List<String>>> peptideModificationHash;
    Map<String, String> peptideIdAndSequenceHash;
    Map<String, List<List<Object>>> spectrumInformationHash;
    Map<String, String> dbReferenceHash;
    /*
     * Added by FG MzIdentML elements
     */
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;

    // The constructor takes the inputFile name and the name of the search engine
    public MzIdentMLReader(String xmlFileName, String searchEngine) {

        inputMzIdentMLFile = xmlFileName;
        this.searchEngine = searchEngine;

        if (searchEngine.equals("mascot")) {
            searchEngineIdentifierExpectation = "Mascot:expectation value";
            searchEngineIdentifierScore = "Mascot:score";
            mascotDecoyTag = null;
        } else if (searchEngine.equals("omssa")) {
            searchEngineIdentifierExpectation = "OMSSA:evalue";
            searchEngineIdentifierScore = "OMSSA:pvalue";
            mascotDecoyTag = null;
        } else if (searchEngine.equals("X!Tandem")) {
            searchEngineIdentifierExpectation = "X\\!Tandem:expect";
            searchEngineIdentifierScore = "X\\!Tandem:hyperscore";
            mascotDecoyTag = null;
        }

        peptideIdAndSequenceHash = new HashMap<String, String>();
        peptideModificationHash = new HashMap<String, List<List<String>>>();
        spectrumInformationHash = new HashMap<String, List<List<Object>>>();
        dbReferenceHash = new HashMap<String, String>();
    }

    // Overloaded constructor for Mascot, requires an identifier to distinguish between Decoy or normal sequence
    public MzIdentMLReader(String xmlFileName, String searchEngine, String mascotDecoyIdentifier) {

        inputMzIdentMLFile = xmlFileName;
        this.searchEngine = searchEngine;

        if (searchEngine.equals("mascot")) {
            searchEngineIdentifierExpectation = "Mascot:expectation value";
            searchEngineIdentifierScore = "Mascot:score";
            mascotDecoyTag = mascotDecoyIdentifier;
        } else if (searchEngine.equals("omssa")) {
            searchEngineIdentifierExpectation = "OMSSA:evalue";
            searchEngineIdentifierScore = "OMSSA:pvalue";
            mascotDecoyTag = null;
        } else if (searchEngine.equals("X!Tandem")) {
            searchEngineIdentifierExpectation = "X\\!Tandem:expect";
            searchEngineIdentifierScore = "X\\!Tandem:hyperscore";
            mascotDecoyTag = null;
        }

        peptideIdAndSequenceHash = new HashMap<String, String>();
        peptideModificationHash = new HashMap<String, List<List<String>>>();
        spectrumInformationHash = new HashMap<String, List<List<Object>>>();
        dbReferenceHash = new HashMap<String, String>();
    }

    //Updated by FG to use the mzIdentML API instead of the manual reading
    public void parseDocument() {
        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzIdentMLFile));
            // Create a hash of Peptide id / peptide sequence
            Iterator<Peptide> iterPeptide = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
            while (iterPeptide.hasNext()) {
                Peptide peptide = iterPeptide.next();
                String pepId = peptide.getId();
                String pepSeq = peptide.getPeptideSequence();
                List<List<String>> modArray = new ArrayList<List<String>>(); // There can be multiple modifications for a peptide
                List<String> tempSingleMod = new ArrayList<String>();
                List<Modification> modificationmList = peptide.getModification();
                for (int i = 0; i < modificationmList.size(); i++) {
                    Modification modification = modificationmList.get(i);
                    String location = modification.getLocation().toString();
                    String residue = "";
                    if(modification.getResidues() != null && !modification.getResidues().isEmpty())
                    	residue = modification.getResidues().get(0);
                    String monoisotopicMassDelta = modification.getMonoisotopicMassDelta().toString();
                    
                    tempSingleMod.add(0, location);
                    tempSingleMod.add(1, residue);
                    tempSingleMod.add(2, monoisotopicMassDelta);
                    
                    List<CvParam> cvParamList = modification.getCvParam();
                    if (cvParamList != null && cvParamList.get(0).getName() != null) {
                        tempSingleMod.add(3, cvParamList.get(0).getName());
                    }
                    modArray.add(new ArrayList<String>(tempSingleMod));
                    tempSingleMod.clear();
                }
                // Peptide ID and sequence pair; Peptide ID and Modification arrays
                peptideIdAndSequenceHash.put(pepId, pepSeq);
                peptideModificationHash.put(pepId, new ArrayList<List<String>>(modArray));
            }

            // Handle the DBSequence references in case of Mascot; in other cases we can ignore 
            // it and leave the hash empty
            if (searchEngine.equals("mascot")) {
                Iterator<DBSequence> iterDBSequence = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
                while (iterDBSequence.hasNext()) {
                    DBSequence dbSeq = iterDBSequence.next();
                    dbReferenceHash.put(dbSeq.getId(), dbSeq.getAccession());
                }
            }
            // SIRs
            List<SpectrumIdentificationResult> sirListTemp = new ArrayList<>();
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                sirListTemp.add(iterSpectrumIdentificationResult.next());

            }

            List<SpectrumIdentificationItem> siiListTemp = new ArrayList<>();
            // loop on SIRs
            for (int z = 0; z < sirListTemp.size(); z++) {
                SpectrumIdentificationResult spectrumIdentificationResult = sirListTemp.get(z);
                
                String spectrumResultId = spectrumIdentificationResult.getId();
                siiListTemp = spectrumIdentificationResult.getSpectrumIdentificationItem();
                List<List<Object>> multiSpectrum = new ArrayList<List<Object>>();
                List<Object> temp = new ArrayList<Object>();
                // loop on SIIs
                for (int i = 0; i < siiListTemp.size(); i++) {
                    String spectrumItemId = new String(), peptideRef = new String();
                    int rank = 0, chargedState = 0;
                    double calculatedMass = 0.0, experimentalMass = 0.0;
                    double expectationValue = 0.0, scoreValue = 0.0;
                    int pepStart = 0, pepEnd = 0;
                    String pepDBSequence_Ref = new String();
                    boolean isDecoy = false;
                    SpectrumIdentificationItem spectrumIdentItem = siiListTemp.get(i);
                    spectrumItemId = spectrumIdentItem.getId();
                    experimentalMass = Double.valueOf(spectrumIdentItem.getExperimentalMassToCharge());
                    calculatedMass = Double.valueOf(spectrumIdentItem.getCalculatedMassToCharge());
                    chargedState = spectrumIdentItem.getChargeState();
                    peptideRef = spectrumIdentItem.getPeptideRef();
                    rank = spectrumIdentItem.getRank();
                  
                    // NOTE - This also has to be a list, as there can be many peptides
                    List<List<Object>> peptideEvd = new ArrayList<List<Object>>();

                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentItem.getPeptideEvidenceRef();
                    // Process <PeptideEvidence>
                    List<Object> singlePeptideEvd = new ArrayList<Object>();
                    for (int j = 0; j < peptideEvidenceRefList.size(); j++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(j);
                        PeptideEvidence peptiedEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());

                        if (peptiedEvidence != null) {
                            pepStart = peptiedEvidence.getStart().intValue();
                            pepEnd = peptiedEvidence.getEnd().intValue();
                            pepDBSequence_Ref = peptiedEvidence.getDBSequenceRef();
                            isDecoy = peptiedEvidence.isIsDecoy();

                            // Special case for mascot to recognise decoy or non-decoy
                            if (searchEngine.equals("mascot")) {
                                if (pepDBSequence_Ref.contains(mascotDecoyTag)) {
                                    isDecoy = true;
                                }
                            }
                            singlePeptideEvd.add(0, pepStart);
                            singlePeptideEvd.add(1, pepEnd);
                            singlePeptideEvd.add(2, isDecoy);
                            if (searchEngine.equals("mascot")) { // Correct the DBSeq_ anomaly in Mascot
                                if (pepDBSequence_Ref.contains("DBSeq_")) {
                                    String dbAccn = dbReferenceHash.get(pepDBSequence_Ref);
                                    singlePeptideEvd.add(3, dbAccn);
                                    //singlePeptideEvd.add(3, pepDBSequence_Ref.replace("DBSeq_",""));
                                } else {
                                    singlePeptideEvd.add(3, pepDBSequence_Ref);
                                }
                            } else {
                                singlePeptideEvd.add(3, pepDBSequence_Ref);
                            }

                            peptideEvd.add(new ArrayList<Object>(singlePeptideEvd));
                            singlePeptideEvd.clear();
                        }

                    }
                    List<CvParam> cvParamList = spectrumIdentItem.getCvParam();

                    for (int j = 0; j < cvParamList.size(); j++) {
                        CvParam cvParam = cvParamList.get(j);
                        String score_expectationIdentifier = cvParam.getName();
                        if (score_expectationIdentifier != null) {
                            if (score_expectationIdentifier.equals(searchEngineIdentifierExpectation)) {
                                expectationValue = Double.parseDouble(cvParam.getValue());
                            }
                            if (score_expectationIdentifier.equals(searchEngineIdentifierScore)) {
                                scoreValue = Double.parseDouble(cvParam.getValue());
                            }
                        } else {
                            System.out.println("Fatal error : Score or Expectation value missing for " + spectrumItemId);
                        }

                    }
                    temp.add(0, spectrumItemId);
                    temp.add(1, peptideRef);
                    temp.add(2, rank);
                    temp.add(3, chargedState);
                    temp.add(4, calculatedMass);
                    temp.add(5, experimentalMass);
                    temp.add(6, expectationValue);
                    temp.add(7, scoreValue);

                    if (!peptideEvd.isEmpty()) {
                        temp.add(8, new ArrayList<List<Object>>(peptideEvd));
                        peptideEvd.clear();
                    }
                }
                if (!temp.isEmpty()) {
                    multiSpectrum.add(new ArrayList<Object>(temp));
                    temp.clear();
                }

                List<CvParam> cvList = spectrumIdentificationResult.getCvParam();
                for (int i = 0; i < cvList.size(); i++) {
                    CvParam cvParam = cvList.get(i);
                    if (cvParam.getName().equals("spectrum title")) {
                        if (searchEngine.equals("mascot")) {
                            for (int k = 0; k < multiSpectrum.size(); k++) {
                                multiSpectrum.get(k).set(0, cvParam.getValue());
                            }
                        }

                    }

                }

                spectrumInformationHash.put(spectrumResultId, new ArrayList<List<Object>>(multiSpectrum));
            }

            // If the search engine is Mascot, call the routine to do the necessary conversions
            //if (searchEngine.equals("mascot")) {
            //    dealWithMascotFile();
            //}

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Convert the MAscot IDs to generic item_xx format; and deal with the
     * DBRef_ strings in the accessions
     */
    private void dealWithMascotFile() {

        String[] specKeys = (String[]) spectrumInformationHash.keySet().toArray(new String[0]);

        for (int k = 0; k < specKeys.length; k++) {
            String specKey = specKeys[k];

            String dtaIdentifier = spectrumInformationHash.get(specKey).get(0).get(0).toString().trim();

            for (int i = 0; i < spectrumInformationHash.get(specKey).size(); i++) {
                String rank = spectrumInformationHash.get(specKey).get(i).get(2).toString();
                String[] tokens = specKey.split("=");
                String replacementSpecKey = "item_ref_" + tokens[1] + "_" + rank;
                spectrumInformationHash.get(specKey).get(i).set(0, replacementSpecKey);
            }

            List<List<Object>> contents = spectrumInformationHash.get(specKey);
            spectrumInformationHash.put(dtaIdentifier, contents);
            spectrumInformationHash.remove(specKey);
        }

    }

    public Map<String, List<List<String>>> getPeptideModification() {
        return peptideModificationHash;
    }

    public Map<String, String> getPeptideSequences() {
        return peptideIdAndSequenceHash;
    }

    public Map<String, List<List<Object>>> getSpectrumInfo() {
        return spectrumInformationHash;
    }

    /**
     * @param args
     */
    public static void main(String[] args) throws Exception {

        String xmlToRead = args[0];
        String searchEngine = args[1];
        String outFile = args[2];

        String mascotDecoyTag = "Rev";
        //String xmlToRead = "exampleMzIDFiles/7merge_mascot.mzid";
        //String searchEngine = "mascot";
        //String outFile = "output/mascot_mzid.txt";


        //String xmlToRead = "exampleMzIDFiles/125merge_omssa.mzid";
        //String searchEngine = "omssa";
        //String outFile = "output/omssa_mzid_junk.txt";

        //String xmlToRead = "exampleMzIDFiles/125merge_tandem.mzid";
        //String searchEngine = "X!Tandem";
        //String outFile = "output/tandem_mzid.txt"; 


        MzIdentMLReader mzRead;
        if (searchEngine.equals("mascot")) {
            mzRead = new MzIdentMLReader(xmlToRead, searchEngine, mascotDecoyTag);
        } else {
            mzRead = new MzIdentMLReader(xmlToRead, searchEngine);
        }

        mzRead.parseDocument();

        Writer out = new BufferedWriter(new FileWriter(outFile));

        Collection<String> spec = mzRead.spectrumInformationHash.keySet();
        Iterator<String> spec_it = spec.iterator();

        while (spec_it.hasNext()) {
            String specKey = spec_it.next().toString();
            out.write(specKey + "\n");
            out.write(mzRead.spectrumInformationHash.get(specKey).toString() + "\n");
        }
        out.write(mzRead.getPeptideSequences().toString());
        out.write("\n\n");
        out.write(mzRead.getPeptideModification().toString());
        out.close();

    }
}
