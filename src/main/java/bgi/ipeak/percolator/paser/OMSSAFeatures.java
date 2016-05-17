/**
 *
 */
package bgi.ipeak.percolator.paser;

import bgi.ipeak.percolator.util.Feature;
import bgi.ipeak.percolator.util.PSMImpl;
import bgi.ipeak.percolator.util.SpectrumBasicInfo;

/**
 * @author Administrator
 *
 */
public class OMSSAFeatures extends SpectrumBasicInfo implements Feature {

    /**
     * Source file name without postfix and path.
     */
    private String fileName = null;
    /**
     * PSM information.
     */
    private PSMImpl iPSM = null;
    /**
     * The neutral logarithm of total intensity.
     */
    private double lnTotalIntensity = 0.0;
    /**
     * Negative logarithm of P-value.
     */
    private double pvalue = 0.0;
    /**
     * The longest matched ion series.
     */
    private int longest = 0;
    /**
     * does the peptide have a C-terminal enzymatic (tryptic) site?
     */
    private boolean entryC = false;
    /**
     * Is the peptide preceded by an enzymatic (tryptic) site?
     */
    private boolean entryN = false;
    /**
     * The number of enzymatic (tryptic) sites.
     */
    private int enzInt = 0;
    /**
     * Score, negative logarithm of E-value.
     */
    private double score = 0.0;
    /**
     * the precursor delta mass of a PSM in dalton.
     */
    private double deltaMass = 0.0;
    /**
     * the absolute precursor delta mass of a PSM in dalton.
     */
    private double absDeltaMass = 0.0;
    /**
     * the isotopic precursor delta mass of a PSM in dalton.
     */
    private double isoDeltaMass = 0.0;
    /**
     * the absolute precursor delta mass of a PSM in ppm.
     */
    private double absDeltaMassPPM = 0.0;
    /**
     * the precursor delta mass of a PSM in ppm.
     */
    private double deltaMassPPM = 0.0;
    /**
     * the isotopic precursor delta mass of a PSM in ppm.
     */
    private double isoDeltaMassPPM = 0.0;
    /**
     * the ratio of sites with variable modification against potential variable
     * modification sites in a peptide.
     */
    private double varMods = 0.0;
    /**
     * the maximum intensity of all matched fragment ions.
     */
    private double maxMatchedIonIntensity = 0.0;
    /**
     * the log-transformed total intensity of all matched fragment ions.
     */
    private double matchedIonTotalInte = 0.0;
    /**
     * the relative total intensity of all matched fragment ions.
     */
    private double relMatchedIonTotalInte = 0.0;
    /**
     * the mean delta mass of all matched fragment ions.
     */
    private double fragDmError = 0.0;
    /**
     * the median delta mass of all matched fragment ions.
     */
    private double fragDmMedian = 0.0;
    /**
     * the interquartile range of delta masses of all matched fragment ions.
     */
    private double fragDmIqr = 0.0;
    /**
     * the median delta mass of all matched fragment ions in ppm.
     */
    private double fragDmMedianPPM = 0.0;
    /**
     * the interquartile range of delta masses of all matched fragment ions in
     * ppm.
     */
    private double fragDmIqrPPM = 0.0;
    /**
     * The number of spectra matched the peptide.
     */
    private int qMatch = 0;
    /**
     * an integer array of the fraction of all matched ion series.
     */
    private int[] fracIonSeries = null;
    /**
     * an integer array of the percentage of intensity of all matched ion
     * series.
     */
    private int[] inteFracIonSeries = null;

    /**
     * @return the fileName
     */
    public String getFileName() {
        return fileName;
    }

    /**
     * @return the iPSM
     */
    public PSMImpl getiPSM() {
        return iPSM;
    }

    /**
     * @param fileName the fileName to set
     */
    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    /**
     * @param iPSM the iPSM to set
     */
    public void setiPSM(PSMImpl iPSM) {
        this.iPSM = iPSM;
    }

    /**
     * @param score the score to set
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**
     * @param deltaMass the deltaMass to set
     */
    public void setDeltaMass(double deltaMass) {
        this.deltaMass = deltaMass;
    }

    /**
     * @param absDeltaMass the absDeltaMass to set
     */
    public void setAbsDeltaMass(double absDeltaMass) {
        this.absDeltaMass = absDeltaMass;
    }

    /**
     * @param isoDeltaMass the isoDeltaMass to set
     */
    public void setIsoDeltaMass(double isoDeltaMass) {
        this.isoDeltaMass = isoDeltaMass;
    }

    /**
     * @param absDeltaMassPPM the absDeltaMassPPM to set
     */
    public void setAbsDeltaMassPPM(double absDeltaMassPPM) {
        this.absDeltaMassPPM = absDeltaMassPPM;
    }

    /**
     * @param deltaMassPPM the deltaMassPPM to set
     */
    public void setDeltaMassPPM(double deltaMassPPM) {
        this.deltaMassPPM = deltaMassPPM;
    }

    /**
     * @param isoDeltaMassPPM the isoDeltaMassPPM to set
     */
    public void setIsoDeltaMassPPM(double isoDeltaMassPPM) {
        this.isoDeltaMassPPM = isoDeltaMassPPM;
    }

    /**
     * @param varMods the varMods to set
     */
    public void setVarMods(double varMods) {
        this.varMods = varMods;
    }

    /**
     * @param maxMatchedIonIntensity the maxMatchedIonIntensity to set
     */
    public void setMaxMatchedIonIntensity(double maxMatchedIonIntensity) {
        this.maxMatchedIonIntensity = maxMatchedIonIntensity;
    }

    /**
     * @param matchedIonTotalInte the matchedIonTotalInte to set
     */
    public void setMatchedIonTotalInte(double matchedIonTotalInte) {
        this.matchedIonTotalInte = matchedIonTotalInte;
    }

    /**
     * @param relMatchedIonTotalInte the relMatchedIonTotalInte to set
     */
    public void setRelMatchedIonTotalInte(double relMatchedIonTotalInte) {
        this.relMatchedIonTotalInte = relMatchedIonTotalInte;
    }

    /**
     * @param fragDmError the fragDmError to set
     */
    public void setFragDmError(double fragDmError) {
        this.fragDmError = fragDmError;
    }

    /**
     * @param fragDmMedian the fragDmMedian to set
     */
    public void setFragDmMedian(double fragDmMedian) {
        this.fragDmMedian = fragDmMedian;
    }

    /**
     * @param fragDmIqr the fragDmIqr to set
     */
    public void setFragDmIqr(double fragDmIqr) {
        this.fragDmIqr = fragDmIqr;
    }

    /**
     * @param fragDmMedianPPM the fragDmMedianPPM to set
     */
    public void setFragDmMedianPPM(double fragDmMedianPPM) {
        this.fragDmMedianPPM = fragDmMedianPPM;
    }

    /**
     * @param fragDmIqrPPM the fragDmIqrPPM to set
     */
    public void setFragDmIqrPPM(double fragDmIqrPPM) {
        this.fragDmIqrPPM = fragDmIqrPPM;
    }

    /**
     * @param fracIonSeries the fracIonSeries to set
     */
    public void setFracIonSeries(int[] fracIonSeries) {
        this.fracIonSeries = fracIonSeries;
    }

    /**
     * @param inteFracIonSeries the inteFracIonSeries to set
     */
    public void setInteFracIonSeries(int[] inteFracIonSeries) {
        this.inteFracIonSeries = inteFracIonSeries;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getScore()
     */
    @Override
    public double getScore() {
        return score;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getDeltaMass()
     */
    @Override
    public double getDeltaMass() {
        return deltaMass;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getAbsDeltaMass()
     */
    @Override
    public double getAbsDeltaMass() {
        return absDeltaMass;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getIsoDeltaMass()
     */
    @Override
    public double getIsoDeltaMass() {
        return isoDeltaMass;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getDeltaMassPPM()
     */
    @Override
    public double getDeltaMassPPM() {
        return deltaMassPPM;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getAbsDeltaMassPPM()
     */
    @Override
    public double getAbsDeltaMassPPM() {
        return absDeltaMassPPM;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getIsoDeltaMassPPM()
     */
    @Override
    public double getIsoDeltaMassPPM() {
        return isoDeltaMassPPM;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getVarMods()
     */
    @Override
    public double getVarMods() {
        return varMods;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getMaxMatchedIonIntensity()
     */
    @Override
    public double getMaxMatchedIonIntensity() {
        return maxMatchedIonIntensity;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getMatchedIonTotalInte()
     */
    @Override
    public double getMatchedIonTotalInte() {
        return matchedIonTotalInte;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getRelMatchedIonTotalInte()
     */
    @Override
    public double getRelMatchedIonTotalInte() {
        return relMatchedIonTotalInte;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFragDmError()
     */
    @Override
    public double getFragDmError() {
        return fragDmError;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFragDmMedian()
     */
    @Override
    public double getFragDmMedian() {
        return fragDmMedian;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFragDmIqr()
     */
    @Override
    public double getFragDmIqr() {
        return fragDmIqr;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFragDmMedianPPM()
     */
    @Override
    public double getFragDmMedianPPM() {
        return fragDmMedianPPM;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFragDmIqrPPM()
     */
    @Override
    public double getFragDmIqrPPM() {
        return fragDmIqrPPM;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getFracIonSeries()
     */
    @Override
    public int[] getFracIonSeries() {
        return fracIonSeries;
    }

    /* (non-Javadoc)
     * @see bgi.org.cn.ipeak.interfaces.Feature#getInteFracIonSeries()
     */
    @Override
    public int[] getInteFracIonSeries() {
        return inteFracIonSeries;
    }

    @Override
    public int getLongest() {
        return longest;
    }

    /**
     * @return the pvalue
     */
    public double getPvalue() {
        return pvalue;
    }

    /**
     * @param pvalue the pvalue to set
     */
    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    /**
     * @param longest the longest to set
     */
    public void setLongest(int longest) {
        this.longest = longest;
    }

    /**
     * @return the entryC
     */
    public boolean isEntryC() {
        return entryC;
    }

    /**
     * @return the entryN
     */
    public boolean isEntryN() {
        return entryN;
    }

    /**
     * @param entryC the entryC to set
     */
    public void setEntryC(boolean entryC) {
        this.entryC = entryC;
    }

    /**
     * @param entryN the entryN to set
     */
    public void setEntryN(boolean entryN) {
        this.entryN = entryN;
    }

    @Override
    public double getLnTotalIntensity() {
        return this.lnTotalIntensity;
    }

    /**
     * @param lnTotalIntensity the lnTotalIntensity to set
     */
    public void setLnTotalIntensity(double lnTotalIntensity) {
        this.lnTotalIntensity = lnTotalIntensity;
    }

    /**
     * @param qMatch the qMatch to set
     */
    public void setqMatch(int qMatch) {
        this.qMatch = qMatch;
    }

    /**
     * @return the qMatch
     */
    public int getqMatch() {
        return qMatch;
    }

    /**
     * @param enzInt the enzInt to set
     */
    public void setEnzInt(int enzInt) {
        this.enzInt = enzInt;
    }

    /**
     * @return the enzInt
     */
    public int getEnzInt() {
        return enzInt;
    }
}
