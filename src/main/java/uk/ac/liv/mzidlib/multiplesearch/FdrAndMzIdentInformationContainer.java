package uk.ac.liv.mzidlib.multiplesearch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FdrAndMzIdentInformationContainer {

    String xmlFileName;
    String searchEngine;
    Map<String, List<List<String>>> peptideModifications;
    Map<String, String> peptideSequence;
    Map<String, List<List<Object>>> spectrumInfo;
    List<String> sorted_spectrumID;
    List<String> sorted_peptideID;
    List<Double> sorted_evalues;
    List<Double> sorted_scores;
    List<String> sorted_decoyornot;
    List<Double> sortedSimpleFDR;
    List<Double> sorted_qValues;
    List<Double> sorted_estimatedFDR;

    FdrAndMzIdentInformationContainer() {
        xmlFileName = new String("");
        searchEngine = new String("");
        peptideModifications = new HashMap<>(1);
        peptideSequence = new HashMap<>(1);
        spectrumInfo = new HashMap<>(1);
        sorted_spectrumID = new ArrayList<>(1);
        sorted_peptideID = new ArrayList<>(1);
        sorted_evalues = new ArrayList<>(1);
        sorted_scores = new ArrayList<>(1);
        sorted_decoyornot = new ArrayList<>(1);
        sortedSimpleFDR = new ArrayList<>(1);
        sorted_qValues = new ArrayList<>(1);
        sorted_estimatedFDR = new ArrayList<>(1);
    }

    public void populateData(String xmlToRead, String searchEngine,
            Map<String, List<List<String>>> pepMod,
            Map<String, String> pepSeq,
            Map<String, List<List<Object>>> specInfo,
            List<String> sorted_spec, List<String> sorted_pepID,
            List<Double> sorted_eval, List<Double> sorted_score,
            List<String> sorted_decoy, List<Double> sorted_FDR,
            List<Double> sorted_qValues, List<Double> sorted_estFDR) {

        this.xmlFileName = xmlToRead;
        this.searchEngine = searchEngine;
        this.peptideModifications = new HashMap<>(pepMod);
        this.peptideSequence = new HashMap<>(pepSeq);
        this.spectrumInfo = new HashMap<>(specInfo);
        this.sorted_spectrumID = new ArrayList<>(sorted_spec);
        this.sorted_peptideID = new ArrayList<>(sorted_pepID);
        this.sorted_evalues = new ArrayList<>(sorted_eval);
        this.sorted_scores = new ArrayList<>(sorted_score);
        this.sorted_decoyornot = new ArrayList<>(sorted_decoy);
        this.sortedSimpleFDR = new ArrayList<>(sorted_FDR);
        this.sorted_qValues = new ArrayList<>(sorted_qValues);
        this.sorted_estimatedFDR = new ArrayList<>(sorted_estFDR);
    }

}
