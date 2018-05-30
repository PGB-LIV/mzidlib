/*
 * Date: 14-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.fdr.FalseDiscoveryRateGlobalReader.java
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package uk.ac.liv.mzidlib.fdr;

import static uk.ac.liv.mzidlib.constants.CvConstants.MASCOT_EXPECTATION_VALUE;
import static uk.ac.liv.mzidlib.constants.CvConstants.MSGF_EVALUE;
import static uk.ac.liv.mzidlib.constants.CvConstants.OMSSA_EVALUE;
import static uk.ac.liv.mzidlib.constants.CvConstants.PROTEINPROSPECTOR_EXPECTATION_VALUE;
import static uk.ac.liv.mzidlib.constants.CvConstants.SEQUEST_EXPECTATION_VALUE;
import static uk.ac.liv.mzidlib.constants.CvConstants.XTANDEM_EXPECT;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.MzIdentMLObject;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinAmbiguityGroup;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionHypothesis;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 14-Aug-2017 11:11:37
 */
class FalseDiscoveryRateGlobalReader {

    private final MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private final String decoy;
    private String fdrLevel;
    private String proteinLevel;

    private AnalysisCollection analysisCollection;
    private AnalysisSoftwareList analysisSoftwareList;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private AuditCollection auditCollection;
    private CvList cvList;
    private Inputs inputs;
    private Provider provider;
    private String searchDatabase_Ref;

    private final Map<String, DBSequence> dBSequenceMap = new HashMap<>();
    private final Map<String, PeptideEvidence> peptideEvidenceMap
            = new HashMap<>();
    private final Map<String, Peptide> peptideMap = new HashMap<>();
    private final Map<String, String> peptideIdAndSequence = new HashMap<>();

    // Store all PSM related information to a specific peptide
    private final Map<String, List<String>> peptidePSMMap = new HashMap<>();
    private final Map<String, Integer> peptidePSMCount = new HashMap<>();
    private final Map<String, String> bestScorePSM = new HashMap<>();

    //TODO: Make this as an array of possible values than a long string
    private final String allowedEvalues = XTANDEM_EXPECT.getAccession() + ";"
            + MASCOT_EXPECTATION_VALUE.getAccession() + ";"
            + SEQUEST_EXPECTATION_VALUE.getAccession() + ";"
            + OMSSA_EVALUE.getAccession() + ";"
            + PROTEINPROSPECTOR_EXPECTATION_VALUE.getAccession() + ";"
            + MSGF_EVALUE.getAccession();

    private boolean usingFileDecoyAttribute = false;

    /*
     * For each Spectrum Result id we need to store :: - peptides associated-
     * scores - evalue - decoy = true/false
     */
    private final List<String> spectrumResult = new ArrayList<>();
    private final List<String> spectrumItem = new ArrayList<>();
    private final List<String> peptideNames = new ArrayList<>();
    private final List<Double> evalues = new ArrayList<>();
    private final List<String> decoyOrNot = new ArrayList<>();

    FalseDiscoveryRateGlobalReader(File file, String decoy, String fdrLevel,
                                   String proteinLevel) {
        this.decoy = decoy;
        this.fdrLevel = fdrLevel;
        this.proteinLevel = proteinLevel;
        if (decoy == null || decoy.isEmpty()) {
            this.usingFileDecoyAttribute = true;
        }
        this.mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(file);
        try {
            readMzIdentML();
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className
                    + "\" was not completed because of " + error.getMessage()
                    + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    /**
     * @return the mzIdentMLUnmarshaller
     */
    public MzIdentMLUnmarshaller getMzIdentMLUnmarshaller() {
        return mzIdentMLUnmarshaller;
    }

    /**
     * @return the decoy
     */
    public String getDecoy() {
        return decoy;
    }

    /**
     * @return the fdrLevel
     */
    public String getFdrLevel() {
        return fdrLevel;
    }

    /**
     * @param fdrLevel the fdrLevel to set
     */
    public void setFdrLevel(String fdrLevel) {
        this.fdrLevel = fdrLevel;
    }

    /**
     * @return the proteinLevel
     */
    public String getProteinLevel() {
        return proteinLevel;
    }

    /**
     * @param proteinLevel the proteinLevel to set
     */
    public void setProteinLevel(String proteinLevel) {
        this.proteinLevel = proteinLevel;
    }

    /**
     * @return the analysisCollection
     */
    public AnalysisCollection getAnalysisCollection() {
        return analysisCollection;
    }

    /**
     * @return the analysisSoftwareList
     */
    public AnalysisSoftwareList getAnalysisSoftwareList() {
        return analysisSoftwareList;
    }

    /**
     * @return the analysisProtocolCollection
     */
    public AnalysisProtocolCollection getAnalysisProtocolCollection() {
        return analysisProtocolCollection;
    }

    /**
     * @return the auditCollection
     */
    public AuditCollection getAuditCollection() {
        return auditCollection;
    }

    /**
     * @return the cvList
     */
    public CvList getCvList() {
        return cvList;
    }

    /**
     * @return the inputs
     */
    public Inputs getInputs() {
        return inputs;
    }

    /**
     * @return the provider
     */
    public Provider getProvider() {
        return provider;
    }

    /**
     * @return the searchDatabase_Ref
     */
    public String getSearchDatabase_Ref() {
        return searchDatabase_Ref;
    }

    /**
     * @return the dBSequenceMap
     */
    public Map<String, DBSequence> getdBSequenceMap() {
        return dBSequenceMap;
    }

    /**
     * @return the peptideMap
     */
    public Map<String, Peptide> getPeptideMap() {
        return peptideMap;
    }

    /**
     * @return the peptideIdAndSequence
     */
    public Map<String, String> getPeptideIdAndSequence() {
        return peptideIdAndSequence;
    }

    /**
     * @return the peptidePSMMap
     */
    public Map<String, List<String>> getPeptidePSMMap() {
        return peptidePSMMap;
    }

    /**
     * @return the peptidePSMCount
     */
    public Map<String, Integer> getPeptidePSMCount() {
        return peptidePSMCount;
    }

    /**
     * @return the bestScorePSM
     */
    public Map<String, String> getBestScorePSM() {
        return bestScorePSM;
    }

    /**
     * @return the allowedEvalues
     */
    public String getAllowedEvalues() {
        return allowedEvalues;
    }

    /**
     * @return the spectrumResult
     */
    public List<String> getSpectrumResult() {
        return spectrumResult;
    }

    /**
     * @return the spectrumItem
     */
    public List<String> getSpectrumItem() {
        return spectrumItem;
    }

    /**
     * @return the peptideNames
     */
    public List<String> getPeptideNames() {
        return peptideNames;
    }

    /**
     * @return the evalues
     */
    public List<Double> getEvalues() {
        return evalues;
    }

    /**
     * @return the decoyOrNot
     */
    public List<String> getDecoyOrNot() {
        return decoyOrNot;
    }

    /**
     * @return the peptideEvidenceMap
     */
    public Map<String, PeptideEvidence> getPeptideEvidenceMap() {
        return peptideEvidenceMap;
    }

    private void readMzIdentML() {
        readMetaData();
        readSequenceCollection();
        switch (getFdrLevel()) {
            case "Peptide":
                readMzIdentMLPeptide();
                break;
            case "ProteinGroup":
                readMzIdentMLProteinGroup();
                break;
            case "PSM":
                readMzIdentMLPSM();
                break;
            default:
                System.err.println("Unknown FDR Level:" + getFdrLevel());
                break;
        }
    }

    private void readMetaData() {
        cvList = getMzIdentMLUnmarshaller().unmarshal(MzIdentMLElement.CvList);
        analysisSoftwareList = getMzIdentMLUnmarshaller().unmarshal(
                MzIdentMLElement.AnalysisSoftwareList);
        auditCollection = getMzIdentMLUnmarshaller().unmarshal(
                MzIdentMLElement.AuditCollection);
        provider = getMzIdentMLUnmarshaller().unmarshal(MzIdentMLElement.Provider);
        analysisProtocolCollection = getMzIdentMLUnmarshaller().unmarshal(
                MzIdentMLElement.AnalysisProtocolCollection);
        analysisCollection = getMzIdentMLUnmarshaller().unmarshal(
                MzIdentMLElement.AnalysisCollection);
        inputs = getMzIdentMLUnmarshaller().unmarshal(MzIdentMLElement.Inputs);
        searchDatabase_Ref = getInputs().getSearchDatabase().get(0).getId();
    }

    private void readSequenceCollection() {
        Iterator<DBSequence> iterDBSequence = getMzIdentMLUnmarshaller().
                unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
        while (iterDBSequence.hasNext()) {
            DBSequence dBSequence = iterDBSequence.next();
            getdBSequenceMap().put(dBSequence.getId(), dBSequence);
        }
        getPeptideEvidenceMap().clear();
        Iterator<PeptideEvidence> iterPeptideEvidence = getMzIdentMLUnmarshaller().
                unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
        while (iterPeptideEvidence.hasNext()) {
            PeptideEvidence pe = iterPeptideEvidence.next();
            getPeptideEvidenceMap().put(pe.getId(), pe);
        }
        Iterator<Peptide> iterPeptide = getMzIdentMLUnmarshaller().
                unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
        while (iterPeptide.hasNext()) {
            Peptide peptide = iterPeptide.next();
            getPeptideMap().put(peptide.getId(), peptide);
            String pepId = peptide.getId();
            String pepSeq = peptide.getPeptideSequence();
            getPeptideIdAndSequence().put(pepId, pepSeq);
        }
    }

    private void readMzIdentMLPeptide() {
        try {

            Iterator<MzIdentMLObject> iterSpectrumIdentificationResult
                    = getMzIdentMLUnmarshaller().unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult
                        = (SpectrumIdentificationResult) iterSpectrumIdentificationResult.
                        next();
                String spectrumResultId = spectrumIdentificationResult.
                        getSpectrumID();

                for (SpectrumIdentificationItem spectrumIdentItem
                        : spectrumIdentificationResult.
                        getSpectrumIdentificationItem()) {
                    List<CvParam> cvParamListSpectrumIdentificationItem
                            = spectrumIdentItem.getCvParam();
                    List<PeptideEvidenceRef> peptideEvidenceRefList
                            = spectrumIdentItem.getPeptideEvidenceRef();
                    boolean isdecoy = false;
                    for (PeptideEvidenceRef peptideEvidenceRef1
                            : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence1
                                = getPeptideEvidenceMap().get(
                                        peptideEvidenceRef1.
                                        getPeptideEvidenceRef());
                        if (peptideEvidence1 != null) {
                            isdecoy = peptideEvidence1.isIsDecoy();
                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence1.isIsDecoy()) {
                                    isdecoy = true;
                                    break;
                                }
                            } else if (getDecoy() != null && !decoy.equals("")
                                    && (getDecoy().length() > 1)) {
                                if (getdBSequenceMap().get(peptideEvidence1.
                                        getDBSequenceRef()).getAccession().
                                        contains(getDecoy())) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                System.out.println(
                                        "Error - no decoy value set, need to use alternative constructor");
                            }
                        }
                    }
                    double combinedFDR = 0;
                    for (CvParam cvParam : cvParamListSpectrumIdentificationItem) {
                        String accession = cvParam.getAccession();

                        if (getAllowedEvalues().contains(accession)) {
                            combinedFDR = Double.valueOf(cvParam.getValue());

                        }
                    }

                    String peptideRef = spectrumIdentItem.getPeptideRef();

                    if (!peptidePSMMap.containsKey(peptideRef)) {
                        List<String> tmp = new ArrayList<>();
                        tmp.add(spectrumIdentItem.getId());
                        getPeptidePSMMap().put(peptideRef, tmp);
                        getPeptidePSMCount().put(peptideRef, 1);
                    } else if (getPeptidePSMMap().containsKey(peptideRef)
                            && !peptidePSMMap.get(peptideRef).contains(
                            spectrumIdentItem.getId())) {
                        getPeptidePSMMap().get(peptideRef).add(spectrumIdentItem.
                                getId());
                        int oldCount = getPeptidePSMCount().get(peptideRef);
                        oldCount = oldCount + 1;
                        getPeptidePSMCount().put(peptideRef, oldCount);
                    }

                    if (!bestScorePSM.containsKey(peptideRef)) {
                        String newKey = spectrumIdentItem.getId() + ":_:"
                                + combinedFDR + ":_:"
                                + Boolean.toString(isdecoy) + ":_:"
                                + spectrumResultId;
                        getBestScorePSM().put(peptideRef, newKey);
                    } else {
                        String value = getBestScorePSM().get(peptideRef);
                        if (value != null && value.contains(":_:")) {
                            String[] psm_fdr = value.split(":_:");
                            if (combinedFDR < Double.valueOf(psm_fdr[1])) {
                                getBestScorePSM().remove(peptideRef);
                                String newKey = spectrumIdentItem.getId()
                                        + ":_:" + combinedFDR + ":_:" + Boolean.
                                        toString(isdecoy) + ":_:"
                                        + spectrumResultId;
                                getBestScorePSM().put(peptideRef, newKey);

                            }
                        }
                    }
                }
            }

            getSpectrumResult().clear();
            getSpectrumItem().clear();
            getPeptideNames().clear();
            getEvalues().clear();
            getDecoyOrNot().clear();

            for (Map.Entry<String, String> entry : getBestScorePSM().entrySet()) {
                String peptideRef = entry.getKey();
                String value = entry.getValue();
                String[] values = value.split(":_:");
                getPeptideNames().add(peptideRef);
                getSpectrumItem().add(values[0]);
                getEvalues().add(Double.valueOf(values[1]));
                getDecoyOrNot().add(values[2]);
                getSpectrumResult().add(values[3]);

            }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className
                    + "\" was not completed because of " + error.getMessage()
                    + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    private void readMzIdentMLProteinGroup() {
        try {

            Iterator<ProteinAmbiguityGroup> iterProteinAmbiguityGroup
                    = getMzIdentMLUnmarshaller()
                    .unmarshalCollectionFromXpath(
                            MzIdentMLElement.ProteinAmbiguityGroup);
            while (iterProteinAmbiguityGroup.hasNext()) {
                ProteinAmbiguityGroup proteinAmbiguityGroup
                        = iterProteinAmbiguityGroup.next();
                boolean anchorProtein = false;
                ProteinDetectionHypothesis anchorProteinDetectionHypothesis
                        = null;

                String proteinAmbiguityGroupId = proteinAmbiguityGroup.getId();
                double combinedFDR = 0;
                boolean isdecoy = false;

                String dbSeqRef;

                for (int j = 0; j < proteinAmbiguityGroup.
                        getProteinDetectionHypothesis().size(); j++) {
                    ProteinDetectionHypothesis proteinDetectionHypothesis
                            = proteinAmbiguityGroup
                            .getProteinDetectionHypothesis().get(j);
                    List<CvParam> cvParamListproteinDetectionHypothesis
                            = proteinDetectionHypothesis.getCvParam();
                    for (CvParam cvParam : cvParamListproteinDetectionHypothesis) {
                        String accession = cvParam.getAccession();
                        if (accession.equals("MS:1001591")) {
                            anchorProteinDetectionHypothesis
                                    = proteinDetectionHypothesis;
                            anchorProtein = true;
                            break;
                        }
                    }
                }
                if (!anchorProtein) {
                    anchorProteinDetectionHypothesis = proteinAmbiguityGroup.
                            getProteinDetectionHypothesis().get(0);

                }
                if (anchorProteinDetectionHypothesis != null) {
                    List<PeptideHypothesis> peptideHypothesisList
                            = anchorProteinDetectionHypothesis.
                            getPeptideHypothesis();
                    dbSeqRef = anchorProteinDetectionHypothesis.
                            getDBSequenceRef();

                    for (PeptideHypothesis peptideHypothesis
                            : peptideHypothesisList) {
                        PeptideEvidence peptideEvidence1
                                = getPeptideEvidenceMap()
                                .get(peptideHypothesis.getPeptideEvidenceRef());
                        if (peptideEvidence1 != null) {
                            isdecoy = peptideEvidence1.isIsDecoy();
                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence1.isIsDecoy()) {
                                    isdecoy = true;
                                    break;
                                }
                            } else if (getDecoy() != null && !decoy.equals("")
                                    && (getDecoy().length() > 1)) {
                                if (getdBSequenceMap().get(peptideEvidence1.
                                        getDBSequenceRef()).getAccession()
                                        .contains(getDecoy())) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                System.out.println(
                                        "Error - no decoy value set, need to use alternative constructor");
                            }
                        }
                    }

                    if (getProteinLevel().equals("PAG")) {
                        for (CvParam cvParam : proteinAmbiguityGroup.
                                getCvParam()) {
                            String accession = cvParam.getAccession();
                            if (getAllowedEvalues().contains(accession)) {
                                combinedFDR = Double.valueOf(cvParam.getValue());
                            }
                        }
                    }
                    if (getProteinLevel().equals("PDH")) {
                        for (CvParam cvParam : anchorProteinDetectionHypothesis.
                                getCvParam()) {
                            String accession = cvParam.getAccession();
                            if (getAllowedEvalues().contains(accession)) {
                                combinedFDR = Double.valueOf(cvParam.getValue());
                            }
                        }

                    }

                    getSpectrumResult().add(proteinAmbiguityGroupId);
                    getSpectrumItem().add(anchorProteinDetectionHypothesis.getId());
                    getPeptideNames().add(dbSeqRef);
                    getEvalues().add(combinedFDR);
                    getDecoyOrNot().add(String.valueOf(isdecoy));
                }
            }

        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className
                    + "\" was not completed because of " + error.getMessage()
                    + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }

    }

    private void readMzIdentMLPSM() {
        try {
            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult
                    = getMzIdentMLUnmarshaller().unmarshalCollectionFromXpath(
                            MzIdentMLElement.SpectrumIdentificationResult);
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult spectrumIdentificationResult
                        = iterSpectrumIdentificationResult.next();
                String spectrumResultId = spectrumIdentificationResult.
                        getSpectrumID();
                for (SpectrumIdentificationItem spectrumIdentItem
                        : spectrumIdentificationResult.
                        getSpectrumIdentificationItem()) {
                    String peptideRef = spectrumIdentItem.getPeptideRef();
                    List<CvParam> cvParamListSpectrumIdentificationItem
                            = spectrumIdentItem.getCvParam();
                    List<PeptideEvidenceRef> peptideEvidenceRefList
                            = spectrumIdentItem.getPeptideEvidenceRef();
                    boolean isdecoy = false;
                    for (PeptideEvidenceRef peptideEvidenceRef1
                            : peptideEvidenceRefList) {
                        PeptideEvidence peptideEvidence1
                                = getPeptideEvidenceMap().get(
                                        peptideEvidenceRef1.
                                        getPeptideEvidenceRef());
                        if (peptideEvidence1 != null) {
                            isdecoy = peptideEvidence1.isIsDecoy();
                            if (usingFileDecoyAttribute) {
                                if (peptideEvidence1.isIsDecoy()) {
                                    isdecoy = true;
                                    break;
                                }
                            } else if (getDecoy() != null && !decoy.equals("")
                                    && (getDecoy().length() > 1)) {
                                if (getdBSequenceMap().get(peptideEvidence1.
                                        getDBSequenceRef()).getAccession().
                                        contains(getDecoy())) {
                                    isdecoy = true;
                                    break;
                                }
                            } else {
                                System.out.println(
                                        "Error - no decoy value set, need to use alternative constructor");
                            }
                        }
                    } // end of for each

                    if (isdecoy) {
                        getDecoyOrNot().add("true");
                    } else {
                        getDecoyOrNot().add("false");

                    }
                    getSpectrumResult().add(spectrumResultId);
                    getPeptideNames().add(peptideRef);
                    getSpectrumItem().add(spectrumIdentItem.getId());
                    for (CvParam cvParam : cvParamListSpectrumIdentificationItem) {
                        String accession = cvParam.getAccession();
                        if (getAllowedEvalues().contains(accession)) {
                            getEvalues().add(new Double(cvParam.getValue()));
                        }

                    }
                }
            }
        } catch (OutOfMemoryError error) {
            String methodName = Thread.currentThread().getStackTrace()[1].
                    getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className
                    + "\" was not completed because of " + error.getMessage()
                    + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

}
