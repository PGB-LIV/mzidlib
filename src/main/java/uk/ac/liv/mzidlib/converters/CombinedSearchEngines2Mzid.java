package uk.ac.liv.mzidlib.converters;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.model.mzidml.*;

/**
 *
 * @author Fawaz Ghali, University of Liverpool, 2012
 */
public class CombinedSearchEngines2Mzid {

    //***************************************** Mascot *****************************************//
    //key = spectrumID e.g. spectrumID="index=4196")
    private Map<String, SpectrumIdentificationResult> mascotSpectrumIdentificationResultHashMap = new HashMap<>();
    //(key = SII)
    private Map<String, SpectrumIdentificationItem> mascotSpectrumIdentificationItemHashMap = new HashMap<>();
    //(key = DBSeq)
    private Map<String, DBSequence> mascotDbSequenceHashMap = new HashMap<>();
    //(key = PeptideEvidence ID)
    private Map<String, PeptideEvidence> mascotPeptideEvidenceHashMap = new HashMap<>();
    //(key = Peptide ID)
    private Map<String, Peptide> mascotPeptideHashMap = new HashMap<>();
    //***************************************** Omssa *****************************************//
    //key = spectrumID e.g. spectrumID="index=4196")
    private Map<String, SpectrumIdentificationResult> omssaSpectrumIdentificationResultHashMap = new HashMap<>();
    //(key = SII)
    private Map<String, SpectrumIdentificationItem> omssaSpectrumIdentificationItemHashMap = new HashMap<>();
    //(key = DBSeq)
    private Map<String, DBSequence> omssaDbSequenceHashMap = new HashMap<>();
    //(key = PeptideEvidence ID)
    private Map<String, PeptideEvidence> omssaPeptideEvidenceHashMap = new HashMap<>();
    //(key = Peptide ID)
    private Map<String, Peptide> omssaPeptideHashMap = new HashMap<>();
    //***************************************** Tandem *****************************************//
    //key = spectrumID e.g. spectrumID="index=4196")
    private Map<String, SpectrumIdentificationResult> tandemSpectrumIdentificationResultHashMap = new HashMap<>();
    //(key = SII)
    private Map<String, SpectrumIdentificationItem> tandemSpectrumIdentificationItemHashMap = new HashMap<>();
    //(key = DBSeq)
    private Map<String, DBSequence> tandemDbSequenceHashMap = new HashMap<>();
    //(key = PeptideEvidence ID)
    private Map<String, PeptideEvidence> tandemPeptideEvidenceHashMap = new HashMap<>();
    //(key = Peptide ID)
    private Map<String, Peptide> tandemPeptideHashMap = new HashMap<>();
    //***************************************** Combined *****************************************//
    //key = spectrumID e.g. spectrumID="index=4196")
    private Map<String, SpectrumIdentificationResult> combinedSpectrumIdentificationResultHashMap = new HashMap<>();
    //(key = SII)
    private Map<String, SpectrumIdentificationItem> combinedSpectrumIdentificationItemHashMap = new HashMap<>();
    //(key = DBSeq)
    private Map<String, DBSequence> combinedDbSequenceHashMap = new HashMap<>();
    //(key = PeptideEvidence ID)
    private Map<String, PeptideEvidence> combinedPeptideEvidenceHashMap = new HashMap<>();
    //(key = Peptide ID)
    private Map<String, Peptide> combinedPeptideHashMap = new HashMap<>();
    private Map<String, String> unimodHashmap = new HashMap<>();

    ;

    //***************************************** setter and getter *****************************************//

    public Map<String, SpectrumIdentificationResult> getMascotSpectrumIdentificationResultHashMap() {
        return mascotSpectrumIdentificationResultHashMap;
    }

    public void setMascotSpectrumIdentificationResultHashMap(Map<String, SpectrumIdentificationResult> mascotSpectrumIdentificationResultHashMap) {
        this.mascotSpectrumIdentificationResultHashMap = mascotSpectrumIdentificationResultHashMap;
    }

    public Map<String, SpectrumIdentificationItem> getMascotSpectrumIdentificationItemHashMap() {
        return mascotSpectrumIdentificationItemHashMap;
    }

    public void setMascotSpectrumIdentificationItemHashMap(Map<String, SpectrumIdentificationItem> mascotSpectrumIdentificationItemHashMap) {
        this.mascotSpectrumIdentificationItemHashMap = mascotSpectrumIdentificationItemHashMap;
    }

    public Map<String, DBSequence> getMascotDbSequenceHashMap() {
        return mascotDbSequenceHashMap;
    }

    public void setMascotDbSequenceHashMap(Map<String, DBSequence> mascotDbSequenceHashMap) {
        this.mascotDbSequenceHashMap = mascotDbSequenceHashMap;
    }

    public Map<String, PeptideEvidence> getMascotPeptideEvidenceHashMap() {
        return mascotPeptideEvidenceHashMap;
    }

    public void setMascotPeptideEvidenceHashMap(Map<String, PeptideEvidence> mascotPeptideEvidenceHashMap) {
        this.mascotPeptideEvidenceHashMap = mascotPeptideEvidenceHashMap;
    }

    public Map<String, Peptide> getMascotPeptideHashMap() {
        return mascotPeptideHashMap;
    }

    public void setMascotPeptideHashMap(Map<String, Peptide> mascotPeptideHashMap) {
        this.mascotPeptideHashMap = mascotPeptideHashMap;
    }

    public Map<String, SpectrumIdentificationResult> getOmssaSpectrumIdentificationResultHashMap() {
        return omssaSpectrumIdentificationResultHashMap;
    }

    public void setOmssaSpectrumIdentificationResultHashMap(Map<String, SpectrumIdentificationResult> omssaSpectrumIdentificationResultHashMap) {
        this.omssaSpectrumIdentificationResultHashMap = omssaSpectrumIdentificationResultHashMap;
    }

    public Map<String, SpectrumIdentificationItem> getOmssaSpectrumIdentificationItemHashMap() {
        return omssaSpectrumIdentificationItemHashMap;
    }

    public void setOmssaSpectrumIdentificationItemHashMap(Map<String, SpectrumIdentificationItem> omssaSpectrumIdentificationItemHashMap) {
        this.omssaSpectrumIdentificationItemHashMap = omssaSpectrumIdentificationItemHashMap;
    }

    public Map<String, DBSequence> getOmssaDbSequenceHashMap() {
        return omssaDbSequenceHashMap;
    }

    public void setOmssaDbSequenceHashMap(Map<String, DBSequence> omssaDbSequenceHashMap) {
        this.omssaDbSequenceHashMap = omssaDbSequenceHashMap;
    }

    public Map<String, PeptideEvidence> getOmssaPeptideEvidenceHashMap() {
        return omssaPeptideEvidenceHashMap;
    }

    public void setOmssaPeptideEvidenceHashMap(Map<String, PeptideEvidence> omssaPeptideEvidenceHashMap) {
        this.omssaPeptideEvidenceHashMap = omssaPeptideEvidenceHashMap;
    }

    public Map<String, Peptide> getOmssaPeptideHashMap() {
        return omssaPeptideHashMap;
    }

    public void setOmssaPeptideHashMap(Map<String, Peptide> omssaPeptideHashMap) {
        this.omssaPeptideHashMap = omssaPeptideHashMap;
    }

    public Map<String, SpectrumIdentificationResult> getTandemSpectrumIdentificationResultHashMap() {
        return tandemSpectrumIdentificationResultHashMap;
    }

    public void setTandemSpectrumIdentificationResultHashMap(Map<String, SpectrumIdentificationResult> tandemSpectrumIdentificationResultHashMap) {
        this.tandemSpectrumIdentificationResultHashMap = tandemSpectrumIdentificationResultHashMap;
    }

    public Map<String, SpectrumIdentificationItem> getTandemSpectrumIdentificationItemHashMap() {
        return tandemSpectrumIdentificationItemHashMap;
    }

    public void setTandemSpectrumIdentificationItemHashMap(Map<String, SpectrumIdentificationItem> tandemSpectrumIdentificationItemHashMap) {
        this.tandemSpectrumIdentificationItemHashMap = tandemSpectrumIdentificationItemHashMap;
    }

    public Map<String, DBSequence> getTandemDbSequenceHashMap() {
        return tandemDbSequenceHashMap;
    }

    public void setTandemDbSequenceHashMap(Map<String, DBSequence> tandemDbSequenceHashMap) {
        this.tandemDbSequenceHashMap = tandemDbSequenceHashMap;
    }

    public Map<String, PeptideEvidence> getTandemPeptideEvidenceHashMap() {
        return tandemPeptideEvidenceHashMap;
    }

    public void setTandemPeptideEvidenceHashMap(Map<String, PeptideEvidence> tandemPeptideEvidenceHashMap) {
        this.tandemPeptideEvidenceHashMap = tandemPeptideEvidenceHashMap;
    }

    public Map<String, Peptide> getTandemPeptideHashMap() {
        return tandemPeptideHashMap;
    }

    public void setTandemPeptideHashMap(Map<String, Peptide> tandemPeptideHashMap) {
        this.tandemPeptideHashMap = tandemPeptideHashMap;
    }

    public Map<String, SpectrumIdentificationResult> getCombinedSpectrumIdentificationResultHashMap() {
        return combinedSpectrumIdentificationResultHashMap;
    }

    public void setCombinedSpectrumIdentificationResultHashMap(Map<String, SpectrumIdentificationResult> combinedSpectrumIdentificationResultHashMap) {
        this.combinedSpectrumIdentificationResultHashMap = combinedSpectrumIdentificationResultHashMap;
    }

    public Map<String, SpectrumIdentificationItem> getCombinedSpectrumIdentificationItemHashMap() {
        return combinedSpectrumIdentificationItemHashMap;
    }

    public void setCombinedSpectrumIdentificationItemHashMap(Map<String, SpectrumIdentificationItem> combinedSpectrumIdentificationItemHashMap) {
        this.combinedSpectrumIdentificationItemHashMap = combinedSpectrumIdentificationItemHashMap;
    }

    public Map<String, DBSequence> getCombinedDbSequenceHashMap() {
        return combinedDbSequenceHashMap;
    }

    public void setCombinedDbSequenceHashMap(Map<String, DBSequence> combinedDbSequenceHashMap) {
        this.combinedDbSequenceHashMap = combinedDbSequenceHashMap;
    }

    public Map<String, PeptideEvidence> getCombinedPeptideEvidenceHashMap() {
        return combinedPeptideEvidenceHashMap;
    }

    public void setCombinedPeptideEvidenceHashMap(Map<String, PeptideEvidence> combinedPeptideEvidenceHashMap) {
        this.combinedPeptideEvidenceHashMap = combinedPeptideEvidenceHashMap;
    }

    public Map<String, Peptide> getCombinedPeptideHashMap() {
        return combinedPeptideHashMap;
    }

    public void setCombinedPeptideHashMap(Map<String, Peptide> combinedPeptideHashMap) {
        this.combinedPeptideHashMap = combinedPeptideHashMap;
    }

    public CvParam makeCvParam(String accession, String name, String value) {
        CvParam cvParam = new CvParam();
        cvParam.setAccession(accession);
        cvParam.setName(name);
        Cv psiCV = new Cv();
        psiCV.setUri("https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo");
        psiCV.setId("PSI-MS");
//        psiCV.setVersion("2.25.0");
        psiCV.setFullName("PSI-MS");
        cvParam.setCv(psiCV);
        cvParam.setValue(value);
        return cvParam;
    }

    public double getScoreFromSII(SpectrumIdentificationItem sii, String cvParamAccForScore) {

        double score = 0.0;

        for (CvParam cvParam : sii.getCvParam()) {
            if (cvParam.getAccession().equals(cvParamAccForScore)) {
                score = Double.parseDouble(cvParam.getValue());
            }
        }

        return score;

    }
    /*
     * Accepts the SIIs associated with any SIR and inserts the correct rank
     * values
     *
     */

    public void addSIIToListAndSetRank(List<SpectrumIdentificationItem> siiList, SpectrumIdentificationItem sii, String scoreCvParamNameToOrderBy, Boolean orderLowToHigh, String sirID) {

        double scoreOfNewSII = getScoreFromSII(sii, scoreCvParamNameToOrderBy);
        final String cvParamScore = scoreCvParamNameToOrderBy;
        final boolean lowToHigh = orderLowToHigh;

        siiList.add(sii);
        Collections.sort(siiList, new Comparator<SpectrumIdentificationItem>() {

            @Override
            public int compare(SpectrumIdentificationItem sii1, SpectrumIdentificationItem sii2) {
                double sii1Score = getScoreFromSII(sii1, cvParamScore);
                double sii2Score = getScoreFromSII(sii2, cvParamScore);
                int i = 0;

                if (lowToHigh) {
                    if (sii1Score < sii2Score) {
                        i = -1;
                    } else if (sii1Score > sii2Score) {
                        i = +1;
                    } else {
                        i = 0;
                    }
                } else {
                    if (sii2Score < sii1Score) {
                        i = -1;
                    } else if (sii2Score > sii1Score) {
                        i = +1;
                    } else {
                        i = 0;
                    }
                }

                //write code to compare next by alphabetical order ... complete putting PDHs into correct PAGs etc...PAGs
                return i;
            }
        });

        int rank = 0;
        int innerRankCounter = 2;
        double lastScore = -999.0;

        for (SpectrumIdentificationItem spii : siiList) {
            double score = getScoreFromSII(spii, scoreCvParamNameToOrderBy);
            boolean sameScore = false;
            if (score != lastScore) { //Otherwise set the same rank as previous
                rank++;
                sameScore = true;
                spii.setId(sirID + "_SII_" + rank);
                spii.setRank(rank);
                innerRankCounter = 2;
            } else {
                spii.setId(sirID + "SII_" + rank + "_" + innerRankCounter);
                spii.setRank(rank);
                innerRankCounter++;
            }

        }

    }

    public Map<String, String> getUnimodHashmap() {
        return unimodHashmap;
    }
}
