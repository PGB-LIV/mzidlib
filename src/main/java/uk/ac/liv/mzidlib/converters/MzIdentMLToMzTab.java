package uk.ac.liv.mzidlib.converters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.xml.bind.JAXBException;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.MzIdentMLObject;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author Fawaz Ghali, 2014
 */
public class MzIdentMLToMzTab {

    private String input = null;
    private String output = null;
    // private MzIdentML mzIdentML;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private String search_engine = null;
    private String search_engine_score = null;
    private HashMap<String, SearchDatabase> searchDatabaseHashMap = new HashMap();

    private List fixedMods = new ArrayList();
    private List variableMods = new ArrayList();

    public static void main(String args[]) {

        MzIdentMLToMzTab mzIdentMLToMzTab = new MzIdentMLToMzTab(args[0], args[1]);
    }

    public MzIdentMLToMzTab(String input, String output) {
        this.input = input;
        this.output = output;
        Writer out = null;

        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(input));
            System.out.println("Reading the input file: " + input);
            out = new BufferedWriter(new FileWriter(output));
            System.out.println("Writing the output file: " + output);
            out.write(writeHeader());
            out.write(writePRT());
            out.write(writePSM());
            out.close();
            System.out.println("Done.");
        } catch (IOException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } finally {
            try {
                out.close();
            } catch (IOException ex) {
                String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
                String className = this.getClass().getName();
                String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                        + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
                System.out.println(message);
            }
        }

    }

    private String writeHeader() {
        StringBuilder header = new StringBuilder();
        header.append("COM\tReport of a \"Summary Identification report\"\n");
        header.append("MTD\tmzTab-version\t1.0.0\n");
        header.append("MTD\tmzTab-mode\tSummary\n");
        header.append("MTD\tmzTab-type\tIdentification\n");
        header.append("MTD\tdescription	mzTab example file for reporting a summary report of identification data\n");

//             MTD	ms_run[1]-location	file://C:/User/name/data/xyc_heat_1.mzML																
        SpectraData spectraData = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SpectraData);
        if (spectraData != null && spectraData.getLocation() != null) {
            header.append("MTD\tms_run[1]-location\t file://").append(spectraData.getLocation().replace("\\", "/")).append("\n");

        }
        String fixed_mod = "";
        String variable_mod = "";
        int variable_mod_count = 0;
        int fixed_mod_count = 0;
        AnalysisProtocolCollection analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);

        List<SpectrumIdentificationProtocol> spectrumIdentificationProtocolList = analysisProtocolCollection.getSpectrumIdentificationProtocol();
        for (int z = 0; z < spectrumIdentificationProtocolList.size(); z++) {
            SpectrumIdentificationProtocol spectrumIdentificationProtocol = spectrumIdentificationProtocolList.get(z);
            ModificationParams modificationParams = spectrumIdentificationProtocol.getModificationParams();
            List<SearchModification> searchModificationList = modificationParams.getSearchModification();
            for (int i = 0; i < searchModificationList.size(); i++) {
                SearchModification searchModification = searchModificationList.get(i);
                List<CvParam> cvParamList = searchModification.getCvParam();
                for (int j = 0; j < cvParamList.size(); j++) {
                    CvParam cvParam = cvParamList.get(j);
                    if (searchModification.isFixedMod()) {
                        fixed_mod_count = fixed_mod_count + 1;
                        String newEntry = "MTD\tfixed_mod[" + fixed_mod_count + "]\t[UNIMOD, " + cvParam.getAccession() + ", " + cvParam.getName() + ", ]\n";
                        fixed_mod = fixed_mod + newEntry;
                        fixedMods.add(cvParam.getName());
                    } else {
                        variable_mod_count = variable_mod_count + 1;
                        String newEntry = "MTD\tvariable_mod[" + variable_mod_count + "]\t[UNIMOD, " + cvParam.getAccession() + ", " + cvParam.getName() + ", ]\n";
                        variable_mod = variable_mod + newEntry;
                        variableMods.add(cvParam.getName());
                    }

                }
            }
        }

        header.append(fixed_mod);
        header.append(variable_mod);

        Inputs inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
        List<SearchDatabase> searchDatabaseList = inputs.getSearchDatabase();
        for (int i = 0; i < searchDatabaseList.size(); i++) {
            SearchDatabase searchDatabase = searchDatabaseList.get(i);
            searchDatabaseHashMap.put(searchDatabase.getId(), searchDatabase);
        }

        AnalysisSoftwareList analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
        if (analysisSoftwareList != null) {
            List<AnalysisSoftware> allASList = analysisSoftwareList.getAnalysisSoftware();
            for (int i = 0; i < allASList.size(); i++) {
                AnalysisSoftware analysisSoftware = allASList.get(i);
                if (analysisSoftware.getSoftwareName() != null && analysisSoftware.getSoftwareName().getCvParam() != null) {
                    search_engine = analysisSoftware.getSoftwareName().getCvParam().getAccession() + ","
                            + analysisSoftware.getSoftwareName().getCvParam().getName() + ",";
                }
                break;
            }
        }

        SequenceCollection sequenceCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
        List<Peptide> peptideList = sequenceCollection.getPeptide();
        for (int i = 0; i < peptideList.size(); i++) {
            Peptide peptide = peptideList.get(i);

        }
        return header.toString();
    }

    private String writePRT() {
        StringBuilder prt = new StringBuilder();
        prt.append("PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tbest_search_engine_score\tsearch_engine_score_ms_run[1]\tambiguity_members\tmodifications\tprotein_coverage\n");
        try {
            String ambiguity_members = "";
            Iterator<MzIdentMLObject> iterProteinAmbiguityGroup = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {
                ProteinAmbiguityGroup proteinAmbiguityGroup = (ProteinAmbiguityGroup) iterProteinAmbiguityGroup.next();
                ProteinDetectionHypothesis anchorProteinDetectionHypothesis = null;
                boolean anchorProtein = false;
                for (int j = 0; j < proteinAmbiguityGroup.getProteinDetectionHypothesis().size(); j++) {
                    ProteinDetectionHypothesis proteinDetectionHypothesis = proteinAmbiguityGroup.getProteinDetectionHypothesis().get(j);
                    List<CvParam> cvParamListproteinDetectionHypothesis = proteinDetectionHypothesis.getCvParam();
                    for (int s = 0; s < cvParamListproteinDetectionHypothesis.size(); s++) {
                        CvParam cvParam = cvParamListproteinDetectionHypothesis.get(s);
                        String accession = cvParam.getAccession();
                        if (accession.equals("MS:1001591")) {
                            anchorProteinDetectionHypothesis = proteinDetectionHypothesis;
                            anchorProtein = true;
                            break;
                        }
                    }
                }
                if (!anchorProtein) {
                    anchorProteinDetectionHypothesis = proteinAmbiguityGroup.getProteinDetectionHypothesis().get(0);

                }

                for (int j = 0; j < proteinAmbiguityGroup.getProteinDetectionHypothesis().size(); j++) {
                    ProteinDetectionHypothesis proteinDetectionHypothesis = proteinAmbiguityGroup.getProteinDetectionHypothesis().get(j);
                    String dBSequenceRef = proteinDetectionHypothesis.getDBSequenceRef();
                    if (!dBSequenceRef.equals(anchorProteinDetectionHypothesis.getDBSequenceRef())) {
                        DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, dBSequenceRef);
                        ambiguity_members = ambiguity_members + dBSequence.getAccession() + ",";
                    }

                }
                if (!ambiguity_members.equals("")) {
                    ambiguity_members = ambiguity_members.substring(0, ambiguity_members.length() - 1);
                }else{
                    ambiguity_members = "null";
                }

                String modifications = "";
                List<PeptideHypothesis> peptideHypothesisList = anchorProteinDetectionHypothesis.getPeptideHypothesis();
                for (int i = 0; i < peptideHypothesisList.size(); i++) {
                    PeptideHypothesis peptideHypothesis = peptideHypothesisList.get(i);
                    String peptideEvidenceRef = peptideHypothesis.getPeptideEvidenceRef();
                    PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef);
                    int start = peptideEvidence.getStart().intValue();
                    String peptideRef = peptideEvidence.getPeptideRef();
                    Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, peptideRef);
                    List<Modification> modificationsList = peptide.getModification();
                    for (int j = 0; j < modificationsList.size(); j++) {
                        Modification modification = modificationsList.get(j);
                        int location = modification.getLocation().intValue();
                        List<CvParam> cvParamList = modification.getCvParam();
                        for (int k = 0; k < cvParamList.size(); k++) {
                            CvParam cvParam = cvParamList.get(k);
                            if (cvParam.getAccession().startsWith("UNIMOD")) {
                                if (!fixedMods.contains(cvParam.getName())) {
                                    int pos = start + location;
                                    modifications = modifications + pos + "-" + cvParam.getAccession() + ",";
                                }
                            }

                        }

                    }
                }
                if (!modifications.equals("")) {
                    modifications = modifications.substring(0, modifications.length() - 1);
                }

                String dBSequenceRef = anchorProteinDetectionHypothesis.getDBSequenceRef();
                DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, dBSequenceRef);
                String proteinDesc = null;
                List<CvParam> cvParamList = dBSequence.getCvParam();
                for (int i = 0; i < cvParamList.size(); i++) {
                    CvParam cvParam = cvParamList.get(i);
                    if (cvParam.getAccession().equals("MS:1001088")) {
                        proteinDesc = cvParam.getValue();
                    }
                }
                String searchDatabaseName = null;
                String searchDatabaseVersion = null;
                String searchDatabaseRef = dBSequence.getSearchDatabaseRef();
                SearchDatabase searchDatabase = searchDatabaseHashMap.get(searchDatabaseRef);
                if (searchDatabase != null) {
                    Param param = searchDatabase.getDatabaseName();
                    if (param != null) {
                        UserParam userParam = param.getUserParam();
                        if (userParam != null) {
                            searchDatabaseName = userParam.getValue();
                        }
                    }
                    if (searchDatabaseName == null) {
                        searchDatabaseName = searchDatabase.getId();
                    }
                    searchDatabaseVersion = searchDatabase.getVersion();
                }

                prt.append("PRT\t");
                prt.append(dBSequence.getAccession()).append("\t");
                prt.append(proteinDesc).append("\t");
                prt.append("null\t");//taxid
                prt.append("null\t");//species
                prt.append(searchDatabaseName).append("\t");//database
                prt.append(searchDatabaseVersion).append("\t");//database_version
                prt.append("[MS,").append(search_engine).append("]\t");//search_engine
                prt.append("null\t");//best_search_engine_score
                prt.append("null\t");//search_engine_score_ms_run[1]
                prt.append(ambiguity_members).append("\t");//ambiguity_members 
                prt.append(modifications).append("\t");//modifications
                prt.append("null\t");//protein_coverage
                prt.append("\n");
            }

        } catch (JAXBException ex) {
            ex.printStackTrace();
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

        return prt.toString();
    }

    private String writePSM() {
        StringBuilder psm = new StringBuilder();
        psm.append("PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score\tmodifications\tspectra_ref\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tpre\tpost\tstart\tend\n");
        int pSM_ID = 1;
        try {
            Iterator<MzIdentMLObject> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.next();
                List<SpectrumIdentificationItem> spectrumIdentificationItemList = spectrumIdentificationResult.getSpectrumIdentificationItem();
                for (int t = 0; t < spectrumIdentificationItemList.size(); t++) {
                    SpectrumIdentificationItem spectrumIdentificationItem = spectrumIdentificationItemList.get(t);
                    String peptideRef = spectrumIdentificationItem.getPeptideRef();
                    Peptide peptide = mzIdentMLUnmarshaller.unmarshal(Peptide.class, peptideRef);

                    String modifications = "";
                    List<PeptideEvidenceRef> peptideEvidenceRefList = spectrumIdentificationItem.getPeptideEvidenceRef();
                    for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
                        PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                        int start = peptideEvidence.getStart().intValue();
                        String peptideRef1 = peptideEvidence.getPeptideRef();
                        Peptide peptide1 = mzIdentMLUnmarshaller.unmarshal(Peptide.class, peptideRef1);
                        List<Modification> modificationsList = peptide1.getModification();
                        for (int j = 0; j < modificationsList.size(); j++) {
                            Modification modification = modificationsList.get(j);
                            int location = modification.getLocation().intValue();
                            List<CvParam> cvParamList = modification.getCvParam();
                            for (int k = 0; k < cvParamList.size(); k++) {
                                CvParam cvParam = cvParamList.get(k);
                                if (cvParam.getAccession().startsWith("UNIMOD")) {
                                    int pos = start + location;
                                    modifications = modifications + pos + "-" + cvParam.getAccession() + ",";
                                }

                            }

                        }
                    }
                    if (!modifications.equals("")) {
                        modifications = modifications.substring(0, modifications.length() - 1);
                    }
                    search_engine_score = "";
                    List<CvParam> cvParamList = spectrumIdentificationItem.getCvParam();
                    for (int i = 0; i < cvParamList.size(); i++) {
                        CvParam cvParam = cvParamList.get(i);
                        search_engine_score = search_engine_score + "[MS,";
                        search_engine_score = search_engine_score + cvParam.getAccession() + "," + cvParam.getName() + "," + cvParam.getValue();

                        if (i == cvParamList.size() - 1) {
                            search_engine_score = search_engine_score + "]";
                        } else {
                            search_engine_score = search_engine_score + "]|";
                        }
                    }
                    for (int i = 0; i < peptideEvidenceRefList.size(); i++) {
                        PeptideEvidenceRef peptideEvidenceRef = peptideEvidenceRefList.get(i);
                        PeptideEvidence peptideEvidence = mzIdentMLUnmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                        psm.append("PSM\t");
                        psm.append(peptide.getPeptideSequence()).append("\t");//Peptide Sequence
                        //psm.append(spectrumIdentificationItem.getId() + "\t");//SII ID
                        psm.append(pSM_ID).append("\t");//SII ID
                        pSM_ID = pSM_ID + 1;
                        String dbSeqRef = peptideEvidence.getDBSequenceRef();
                        DBSequence dBSequence = mzIdentMLUnmarshaller.unmarshal(DBSequence.class, dbSeqRef);

                        psm.append(dBSequence.getAccession()).append("\t");//Accession
                        psm.append("null\t");//unique
                        String searchDatabaseName = null;
                        String searchDatabaseVersion = null;
                        String searchDatabaseRef = dBSequence.getSearchDatabaseRef();
                        SearchDatabase searchDatabase = searchDatabaseHashMap.get(searchDatabaseRef);
                        if (searchDatabase != null) {
                            Param param = searchDatabase.getDatabaseName();
                            if (param != null) {
                                UserParam userParam = param.getUserParam();
                                if (userParam != null) {
                                    searchDatabaseName = userParam.getValue();
                                }
                            }
                            if (searchDatabaseName == null) {
                                searchDatabaseName = searchDatabase.getId();
                            }
                            searchDatabaseVersion = searchDatabase.getVersion();
                        }

                        psm.append(searchDatabaseName).append("\t");//searchDatabaseName
                        psm.append(searchDatabaseVersion).append("\t");//database version
                        psm.append("[MS,").append(search_engine).append("]\t");//search_engine
                        psm.append(search_engine_score).append("\t");//search_engine score
                        psm.append(modifications).append("\t");//modications  

                        psm.append("ms_run[1]:").append(spectrumIdentificationResult.getSpectrumID()).append("\t");//spectra ref

                        List<CvParam> cvParamListSir = spectrumIdentificationResult.getCvParam();
                        String rt = null;
                        for (int j = 0; j < cvParamListSir.size(); j++) {
                            CvParam cvParam = cvParamListSir.get(j);
                            if (cvParam.getAccession().equals("MS:1000016")) {
                                rt = cvParam.getValue();
                                break;
                            }

                        }
                        psm.append(rt).append("\t");//rt
                        psm.append(spectrumIdentificationItem.getChargeState()).append("\t");//ChargeState
                        psm.append(spectrumIdentificationItem.getCalculatedMassToCharge()).append("\t");//CalculatedMassToCharge
                        psm.append(spectrumIdentificationItem.getExperimentalMassToCharge()).append("\t");//ExperimentalMassToCharge
                        psm.append(peptideEvidence.getPre()).append("\t");//pre
                        psm.append(peptideEvidence.getPost()).append("\t");//post
                        psm.append(peptideEvidence.getEnd()).append("\t");//end
                        psm.append(peptideEvidence.getStart()).append("\t"); //start
                        psm.append("\n");
                    }

                }
            }
        } catch (JAXBException ex) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 04 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
        return psm.toString();
    }
}
