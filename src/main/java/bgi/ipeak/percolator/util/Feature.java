package bgi.ipeak.percolator.util;

public interface Feature {
	/**
	 * Return the score of a PSM.
	 * @return the score of a PSM.
	 */
	public double getScore();
	/**
	 * Return the precursor delta mass of a PSM in dalton.
	 * @return the precursor delta mass of a PSM in dalton.
	 */
	public double getDeltaMass();
	/**
	 * Return the absolute precursor delta mass of a PSM in dalton.
	 * @return the absolute precursor delta mass of a PSM in dalton.
	 */
	public double getAbsDeltaMass();
	/**
	 * Return the isotopic precursor delta mass of a PSM in dalton.
	 * @return the isotopic precursor delta mass of a PSM in dalton.
	 */
	public double getIsoDeltaMass();
	/**
	 * Return the precursor delta mass of a PSM in ppm.
	 * @return the precursor delta mass of a PSM in ppm.
	 */
	public double getDeltaMassPPM();
	/**
	 * Return the absolute precursor delta mass of a PSM in ppm.
	 * @return the absolute precursor delta mass of a PSM in ppm.
	 */
	public double getAbsDeltaMassPPM();
	/**
	 * Return the isotopic precursor delta mass of a PSM in ppm.
	 * @return the isotopic precursor delta mass of a PSM in ppm.
	 */
	public double getIsoDeltaMassPPM();
	/**
	 * Return the ratio of sites with variable modification against potential variable modification sites in a peptide.
	 * @return the ratio of sites with variable modification against potential variable modification sites in a peptide.
	 */
	public double getVarMods();
	/**
	 * The neutral logarithm of total intensity.
	 * @return double of neutral logarithm of total intensity.
	 */
	public double getLnTotalIntensity();
	/**
	 * Return the maximum intensity of all matched fragment ions.
	 * @return the maximum intensity of all matched fragment ions.
	 */
	public double getMaxMatchedIonIntensity();
	/**
	 * Return the log-transformed total intensity of all matched fragment ions.
	 * @return the log-transformed total intensity of all matched fragment ions.
	 */
	public double getMatchedIonTotalInte();
	/**
	 * Return the relative total intensity of all matched fragment ions.
	 * matchedIonTotalInte divided by total intensity of the spectrum.
	 * @return the relative total intensity of all matched fragment ions.
	 */
	public double getRelMatchedIonTotalInte();
	/**
	 * Return the mean delta mass of all matched fragment ions.
	 * @return the mean delta mass of all matched fragment ions.
	 */
	public double getFragDmError();
	/**
	 * Return the median delta mass of all matched fragment ions.
	 * @return the median delta mass of all matched fragment ions.
	 */
	public double getFragDmMedian();
	/**
	 * Return the interquartile range of delta masses of all matched fragment ions.
	 * @return the interquartile range of delta masses of all matched fragment ions.
	 */
	public double getFragDmIqr();
	/**
	 * Return the median delta mass of all matched fragment ions in ppm.
	 * @return the median delta mass of all matched fragment ions in ppm.
	 */
	public double getFragDmMedianPPM();
	/**
	 * Return the interquartile range of delta masses of all matched fragment ions in ppm.
	 * @return the interquartile range of delta masses of all matched fragment ions in ppm.
	 */
	public double getFragDmIqrPPM();
	/**
	 * Return the longest matched ion series.
	 * @return the longest matched ion series.
	 */
	public int getLongest();
	/**
	 * Return an integer array of the fraction of all matched ion series.
	 * The fraction equals to the number of matched fragment ions derived from the same ion series, e.g. b, y, 
	 * divided by the peptide length minus one.
	 * @return an integer array of the fraction of all matched ion series.
	 */
	public int[] getFracIonSeries();
	/**
	 * Return an integer array of the percentage of intensity of all matched ion series.
	 * The total intensity of each ion series is divided by that of all matched fragment ions.
	 * @return an integer array of the percentage of intensity of all matched ion series.
	 */
	public int[] getInteFracIonSeries();
}
