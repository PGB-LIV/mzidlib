package uk.ac.liv.mzidlib.multiplesearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SortedMap;
/**
 * Use the TreeMap to sort the values, and return their "indices"
 *
 * @authour Ritesh
 * @date May 5, 2010
 */
import java.util.TreeMap;

public class TreeSortForIndices {

    /*
     * Sort out the information which we extracted based on the e-value/Score
     * The function should sort the data column provided and also return the
     * "original" indices of the sorted items. The indices column can be used to
     * associate corresponding information from spectrumResult, peptideNames
     * etc.
     */
    public Integer[] sortTheValueColumn(Object[] computedValues,
            boolean betterScoresAreLower) {

        SortedMap<Double, List<Integer>> map = new TreeMap<Double, List<Integer>>();
        for (int i = 0; i < computedValues.length; i++) {
            List<Integer> ind = map.get(computedValues[i]);
            if (ind == null) {
                ind = new ArrayList<Integer>();
                map.put((Double) computedValues[i], ind);

            }
            ind.add(i);
        }

        // Flatten the list
        List<Integer> indices = new ArrayList<Integer>();
        for (List<Integer> arr : map.values()) {
            indices.addAll(arr);
        }

        // 14/08/2013 Added by FG to check if scores are lower or not
        if (!betterScoresAreLower) {
            Collections.reverse(indices);
        }

        return (indices.toArray(new Integer[0]));
    }

}
