/**
 *
 */
package bgi.ipeak.percolator.util;

/**
 * @author Administrator
 *
 */
public interface PSMInf {

    /**
     * Return the index of the PSM, the same as the index in MS/MS spectrum.
     *
     * @return integer the index of PSM.
     */
    public int getPSMId();

    /**
     * Return the matched peptide.
     *
     * @return String for the matched peptide.
     */
    public PeptideInfor getPeptide();

    /**
     * Return the rank of this match.
     *
     * @return integer of the rank of this match.
     */
    public int getRank();

    /**
     * Does the PSM map to the target database?
     *
     * @return true if the PSM maps to the target database.
     */
    public boolean isTarget();

    /**
     * Return the score of the PSM.
     *
     * @return double for the score of the PSM.
     */
    public double getScore();

    /**
     * Return the missed cleavages of the PSM.
     *
     * @return integer for the missed cleavages.
     */
    public int getMissedCleavages();

    /**
     * Return the number of matched ions.
     *
     * @return integer for the number of matched ions.
     */
    public int getNumberOfIonsMatched();

    /**
     * Return the performance parameters of MS/MS search.
     *
     * @return the search parameters.
     */
    public Parameters getPerformParameters();

    /**
     * Return the matched ions.
     *
     * @return the matched ions.
     */
    public FragmentIon[] getFragmentIons();
}
