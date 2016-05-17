/**
 *
 */
package bgi.ipeak.percolator.paser;

import bgi.ipeak.percolator.util.PSMImpl;
import bgi.ipeak.percolator.util.SpectrumBasicInfo;

/**
 * @author Administrator
 *
 */
public class MSGFFeatures extends SpectrumBasicInfo {

    /**
     * Source file name without postfix and path.
     */
    private String fileName = null;
    /**
     * PSM information.
     */
    private PSMImpl iPSM = null;
    /**
     * Raw match quality score of MS-GF+.
     */
    private double rawScore = 0.0;
    /**
     * Maximum possible raw match quality score to this spectrum.
     */
    private double deNovoScore = 0.0;
    /**
     * ScoreRatio divided by DeNovoScore
     */
    private double scoreRatio = 0.0;
    /**
     * Difference between RawScore and DeNovoScore
     */
    private double energy = 0.0;
    /**
     * Negative logarithm of the MS-GF+ E Value.
     */
    private double lnEvalue = 0.0;
    /**
     * Negative logarithm of the MS-GF+ Spectral E Value.
     */
    private double lnSpecEvalue = 0.0;
    /**
     * Number of additional neutrons in peptide.
     */
    private double IsotopeError = 0.0;
    /**
     * Summed intensity of identified fragment ions, divided by that of all
     * fragment ions, logged
     */
    private double lnExplainedIonCurrent = 0.0;
    /**
     * Summed intensity of identified N-terminal fragments, divided by that of
     * all identified fragments, logged
     */
    private double lnNTermIonCurrentRatio = 0.0;
    /**
     * Summed intensity of identified N-terminal fragments, divided by that of
     * all identified fragments, logged
     */
    private double lnCTermIonCurrentRatio = 0.0;
    /**
     * Summed intensity of all observed fragment ions, logged.
     */
    private double lnMS2IonCurrent = 0.0;
    /**
     * Difference between theoretical and experimental mass.
     */
    private double deltaMass = 0.0;
    /**
     * Absolute value of the difference between theoretical and experimental
     * mass
     */
    private double absDeltaMass = 0.0;
    /**
     * Mean of mass errors of the seven fragment ion peaks with the highest
     * intensities.
     */
    private double meanErrorTop7 = 0.0;
    /**
     * Squared MeanErrorTop7.
     */
    private double sqMeanErrorTop7 = 0.0;
    /**
     * Standard deviation of mass errors of the seven fragment ion peaks with
     * the highest intensities
     */
    private double stdevErrorTop7 = 0.0;
    /**
     * Boolean, N-terminal agrees with enzymatic cleavage rules.
     */
    private boolean enzN = false;
    /**
     * Boolean, C-terminal agrees with enzymatic cleavage rules.
     */
    private boolean enzC = false;
    /**
     * Number of internal cleavage sites.
     */
    private int enzInt = 0;

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
     * @return the rawScore
     */
    public double getRawScore() {
        return rawScore;
    }

    /**
     * @return the deNovoScore
     */
    public double getDeNovoScore() {
        return deNovoScore;
    }

    /**
     * @return the scoreRatio
     */
    public double getScoreRatio() {
        return scoreRatio;
    }

    /**
     * @return the energy
     */
    public double getEnergy() {
        return energy;
    }

    /**
     * @return the lnEvalue
     */
    public double getLnEvalue() {
        return lnEvalue;
    }

    /**
     * @return the lnSpecEvalue
     */
    public double getLnSpecEvalue() {
        return lnSpecEvalue;
    }

    /**
     * @return the isotopeError
     */
    public double getIsotopeError() {
        return IsotopeError;
    }

    /**
     * @return the lnExplainedIonCurrent
     */
    public double getLnExplainedIonCurrent() {
        return lnExplainedIonCurrent;
    }

    /**
     * @return the lnNTermIonCurrentRatio
     */
    public double getLnNTermIonCurrentRatio() {
        return lnNTermIonCurrentRatio;
    }

    /**
     * @return the lnCTermIonCurrentRatio
     */
    public double getLnCTermIonCurrentRatio() {
        return lnCTermIonCurrentRatio;
    }

    /**
     * @return the lnMS2IonCurrent
     */
    public double getLnMS2IonCurrent() {
        return lnMS2IonCurrent;
    }

    /**
     * @return the deltaMass
     */
    public double getDeltaMass() {
        return deltaMass;
    }

    /**
     * @return the absDeltaMass
     */
    public double getAbsDeltaMass() {
        return absDeltaMass;
    }

    /**
     * @return the meanErrorTop7
     */
    public double getMeanErrorTop7() {
        return meanErrorTop7;
    }

    /**
     * @return the sqMeanErrorTop7
     */
    public double getSqMeanErrorTop7() {
        return sqMeanErrorTop7;
    }

    /**
     * @return the stdevErrorTop7
     */
    public double getStdevErrorTop7() {
        return stdevErrorTop7;
    }

    /**
     * @return the enzN
     */
    public boolean isEnzN() {
        return enzN;
    }

    /**
     * @return the enzC
     */
    public boolean isEnzC() {
        return enzC;
    }

    /**
     * @return the enzInt
     */
    public int getEnzInt() {
        return enzInt;
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
     * @param rawScore the rawScore to set
     */
    public void setRawScore(double rawScore) {
        this.rawScore = rawScore;
    }

    /**
     * @param deNovoScore the deNovoScore to set
     */
    public void setDeNovoScore(double deNovoScore) {
        this.deNovoScore = deNovoScore;
    }

    /**
     * @param scoreRatio the scoreRatio to set
     */
    public void setScoreRatio(double scoreRatio) {
        this.scoreRatio = scoreRatio;
    }

    /**
     * @param energy the energy to set
     */
    public void setEnergy(double energy) {
        this.energy = energy;
    }

    /**
     * @param lnEvalue the lnEvalue to set
     */
    public void setLnEvalue(double lnEvalue) {
        this.lnEvalue = lnEvalue;
    }

    /**
     * @param lnSpecEvalue the lnSpecEvalue to set
     */
    public void setLnSpecEvalue(double lnSpecEvalue) {
        this.lnSpecEvalue = lnSpecEvalue;
    }

    /**
     * @param isotopeError the isotopeError to set
     */
    public void setIsotopeError(double isotopeError) {
        IsotopeError = isotopeError;
    }

    /**
     * @param lnExplainedIonCurrent the lnExplainedIonCurrent to set
     */
    public void setLnExplainedIonCurrent(double lnExplainedIonCurrent) {
        this.lnExplainedIonCurrent = lnExplainedIonCurrent;
    }

    /**
     * @param lnNTermIonCurrentRatio the lnNTermIonCurrentRatio to set
     */
    public void setLnNTermIonCurrentRatio(double lnNTermIonCurrentRatio) {
        this.lnNTermIonCurrentRatio = lnNTermIonCurrentRatio;
    }

    /**
     * @param lnCTermIonCurrentRatio the lnCTermIonCurrentRatio to set
     */
    public void setLnCTermIonCurrentRatio(double lnCTermIonCurrentRatio) {
        this.lnCTermIonCurrentRatio = lnCTermIonCurrentRatio;
    }

    /**
     * @param lnMS2IonCurrent the lnMS2IonCurrent to set
     */
    public void setLnMS2IonCurrent(double lnMS2IonCurrent) {
        this.lnMS2IonCurrent = lnMS2IonCurrent;
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
     * @param meanErrorTop7 the meanErrorTop7 to set
     */
    public void setMeanErrorTop7(double meanErrorTop7) {
        this.meanErrorTop7 = meanErrorTop7;
    }

    /**
     * @param sqMeanErrorTop7 the sqMeanErrorTop7 to set
     */
    public void setSqMeanErrorTop7(double sqMeanErrorTop7) {
        this.sqMeanErrorTop7 = sqMeanErrorTop7;
    }

    /**
     * @param stdevErrorTop7 the stdevErrorTop7 to set
     */
    public void setStdevErrorTop7(double stdevErrorTop7) {
        this.stdevErrorTop7 = stdevErrorTop7;
    }

    /**
     * @param enzN the enzN to set
     */
    public void setEnzN(boolean enzN) {
        this.enzN = enzN;
    }

    /**
     * @param enzC the enzC to set
     */
    public void setEnzC(boolean enzC) {
        this.enzC = enzC;
    }

    /**
     * @param enzInt the enzInt to set
     */
    public void setEnzInt(int enzInt) {
        this.enzInt = enzInt;
    }
}
