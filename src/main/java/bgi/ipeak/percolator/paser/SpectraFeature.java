package bgi.ipeak.percolator.paser;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

/**
 * The features used for OMSSAPercolator.
 *
 * @author Guilin Li.
 *
 */
public class SpectraFeature extends SpectraBaseInfo {

    //Mascot
    /**
     * Spectrum scan number
     */
    private String scan = null;
    /**
     * Traget or decoy library peptide
     */
    private int targetDecoy = -1;
    /**
     * Calculated mass minus observed peptide mass in Da
     */
    private double dM = 0;
    /**
     * @deprecated Charge state
     */
    private double charge1 = 0.0;
    /**
     * @deprecated
     */
    private double charge2 = 0.0;
    /**
     * @deprecated
     */
    private double charge3 = 0.0;
    /**
     * @deprecated
     */
    private double charge4 = 0.0;
    /**
     * @deprecated Mascot score (always on), IonScore
     */
    private double mScore = 0;
    /**
     * @deprecated Mascot score minus Mascot score of next best non-isobaric
     * peptide hit, rank1-rank2 with inequal seq
     */
    private double lgDScore = 0;
    /**
     * Calculated Mr
     */
    private double mrCalc = 0;
    /**
     * Calculated minus observed peptide mass in ppm
     */
    private double dMppm = 0;
    /**
     * Absolute value of calculated minus observed peptide mass in Da
     */
    private double absDM = 0;
    /**
     * Absolute value of calculated minus observed peptide mass in ppm
     */
    private double absDMppm = 0;
    /**
     * Calculated minus observed peptide mass, after eliminating possible
     * isotope errors up to 2 Da, in Da
     */
    private double isoDM = 0;
    /**
     * Calculated minus observed peptide mass, after eliminating possible
     * isotope errors up to 2 Da, in ppm
     */
    private double isoDMppm = 0;
    /**
     * @deprecated Number of missed cleavages (always 0 if no enzyme)
     */
    private int mcNum = 0;
    /**
     * Number of modified sites divided by number of modifiable sites
     */
    private double varmodsRatio = 0;
    /**
     * The number of variable mods used in the peptide. That is, if there are 10
     * Met and 5 of these are oxidised, this counts as 1. A peptide with Met-OX,
     * phosphoS, deamidation, and acetylation, would count as 5.
     */
    private int varmodsCount = 0;
    /**
     * Total number of modifiable sites
     */
    private int modifiable = 0;
    /**
     * Total number of modified residues and terminii
     */
    private int modified = 0;
    /**
     * Log total ion intensity. The 20 most intense peaks in each 100 Da bin are
     * used for all features, and totInt reports this value
     */
    private double totIne = 0;
    /**
     * Log total matched ion intensity
     */
    private double intMatchedTot = 0;
    /**
     * Total matched ion intensity divided by total ion intensity as a
     * percentage (no logs involved)
     */
    private double relIntMatchTot = 0;
    /**
     * Median value of all matched fragment errors in Da
     */
    private double fragDeltaMed = 0;
    /**
     * Interquartile range value of all matched fragment errors in Da
     */
    private double fragDeltaIqr = 0;
    /**
     * Median value of all matched fragment errors in ppm
     */
    private double fragDeltaMedPPM = 0;
    /**
     * Interquartile range value of all matched fragment errors in ppm
     */
    private double fragDeltaIqrPPM = 0;
    /**
     * @deprecated 2nd order polynomial fit to m/z vs delta. Result is RSquared
     * multiplied by the number of points divided by 100
     */
    private double fragDeltaPolyFit = 0;
    /**
     * Longest sequence matched ions, reported separately for each ion series
     * (backbone only), as with fracIonsMatched
     */
    private int longestSeq = 1;
    /**
     * The number of peptide matches for which an ms-ms match was attempted
     */
    private int qmatch = 0;
    /**
     * The peptide string that was matched
     */
    private String pepString = null;
    /**
     * A tab separated list of accessions of proteins that contain this peptide.
     * Must be last feature in list
     */
    private String proteins = null;
    //XTandem
    /**
     * The maximum fragment ion intensity
     */
    private double maxFragIonInt = 0;
    /**
     * @deprecated log10-value of the sum of all the fragment ion intensities.
     */
    private double summedScore = 0;
    /**
     * The log10 value of the expectation value for the peptide identification.
     */
    private double pepEvalue = 0;
    /**
     * @deprecated XTandem's score for the peptide Identification
     */
    private double hyperScore = 0;
    /**
     * @deprecated The HyperScore of the second best peptide match of the
     * spectrum
     */
    private double nextScore = 0;
    /**
     * @deprecated The difference of HyperScore between the best and the second
     * best peptide matches
     */
    private double deltaScore = 0;
    /**
     * The average mass error of all the fragment ions
     */
    private double fragError = 0;
    /**
     * The summation of square of each matched fragment ion m/z error.
     */
    private double sumSquareFragError = 0;
    /**
     * The standard deviation of matched fragment ion m/z error.
     */
    private double stdevFragError = 0;
    /**
     * @deprecated The summed intensities of b-ions
     */
    private double bScore = 0;
    /**
     * @deprecated The summed intensities of y-ions
     */
    private double yScore = 0;
    /**
     * @deprecated The number of peaks that matched between the theoretical and
     * the test mass spectrum, b-ions
     */
    private int bions = 0;
    /**
     * @deprecated The number of peaks that matched between the theoretical and
     * the test mass spectrum, y-ions
     */
    private int yions = 0;
    /**
     * @deprecated The fraction of fragment ions being matched in an ion series
     * (y, b ion series)
     */
    private double bIonFraction = 0;
    /**
     * @deprecated
     */
    private double yIonFraction = 0;
    /**
     * Is the peptide preceded by a N-terminal enzymatic(tryptic) site
     */
    private boolean nTermTryptic = false;
    /**
     * Does the peptide have a C-terminal enzymatic(tryptic) site
     */
    private boolean cTermTryptic = false;
    /**
     * the trypsin site in this peptide.
     */
    private int enzN = 0;
    /**
     * PeptideInfor identification length
     */
    private int pepLen = 0;
    private String mod = null;
    /**
     * Constructed by father class
     *
     * @param spb
     */
    /**
     * @deprecated Matched fragment ion intensity
     */
    private double matchedB = 0.0;
    /**
     * @deprecated
     */
    private double matchedY = 0.0;
    /**
     * @deprecated
     */
    private double matchedX = 0.0;
    /**
     * @deprecated
     */
    private double matchedA = 0.0;
    /**
     * @deprecated
     */
    private double matchedC = 0.0;
    /**
     * @deprecated
     */
    private double matchedZ = 0.0;
    /**
     * P-value for OMSSA
     *
     * @param spb
     */
    private double pvalue = 0.0;

    private String hitID;
    private HashMap<String, Double> ionSeriesInt = new HashMap<String, Double>();
    private HashMap<String, Integer> ionSeriesFrac = new HashMap<String, Integer>();
    private HashMap<String, Integer> outIonInte = new HashMap<String, Integer>();

    public SpectraFeature(SpectraBaseInfo spb) {
        this.setRt(spb.getRt());
        this.setCharge(spb.getCharge());
        this.setIntensity(spb.getIntensity());
        this.setMz(spb.getMz());
    }

    public SpectraFeature() {
    }

    /**
     * calculate the total intensity of the spectrum. bin size=100Da, the 20
     * most intense peak.
     *
     * @return Log total intensity.
     */
    public double calLogTotInt() {
        HashMap<Double, Double> mzToIntensity = new HashMap<Double, Double>();
        for (int j = 0; j < getMz().length; j++) {
            mzToIntensity.put(getMz()[j], getIntensity()[j]);
        }
        double[] mzArray = getMz();
        Arrays.sort(mzArray);
        Vector<Double> intensityMeta = new Vector<Double>();
        double minMZ = mzArray[0];
        double totalInt = 0;
        for (int j = 0; j < mzArray.length; j++) {
            if (mzArray[j] < minMZ + 100) {
                intensityMeta.add(mzToIntensity.get(mzArray[j]));
            } else {
                minMZ += 100;
                if (intensityMeta.size() != 0) {
                    Double[] inteMeta = intensityMeta.toArray(new Double[0]);
                    Arrays.sort(inteMeta);
                    for (int i = inteMeta.length - 1; i >= inteMeta.length - 20 && i >= 0; i--) {
                        totalInt += inteMeta[i];
                    }
                    intensityMeta.clear();
                    j--;
                }
            }
        }
        if (intensityMeta.size() != 0) {
            Double[] inteMeta = intensityMeta.toArray(new Double[0]);
            Arrays.sort(inteMeta);
            for (int i = inteMeta.length - 1; i >= inteMeta.length - 20 && i >= 0; i--) {
                totalInt += inteMeta[i];
            }
            intensityMeta.clear();
        }
        if (totalInt >= Math.E) {
            return Math.log(totalInt) / Math.log(Math.E);
        } else {
            return -1.0;
        }
    }

    public void setScan(String scan) {
        this.scan = scan;
    }

    public String getScan() {
        return scan;
    }

    public double getdM() {
        return dM;
    }

    public void setTargetDecoy(int targetDecoy) {
        this.targetDecoy = targetDecoy;
    }

    public int getTargetDecoy() {
        return targetDecoy;
    }

    public void setdM(double dM) {
        this.dM = dM;
    }

    public double getmScore() {
        return mScore;
    }

    public void setmScore(double mScore) {
        this.mScore = mScore;
    }

    public double getLgDScore() {
        return lgDScore;
    }

    public void setLgDScore(double lgDScore) {
        this.lgDScore = lgDScore;
    }

    public double getMrCalc() {
        return mrCalc;
    }

    public void setMrCalc(double mrCalc) {
        this.mrCalc = mrCalc;
    }

    public double getdMppm() {
        return dMppm;
    }

    public void setdMppm(double dMppm) {
        this.dMppm = dMppm;
    }

    public double getAbsDM() {
        return absDM;
    }

    public void setAbsDM(double absDM) {
        this.absDM = absDM;
    }

    public double getAbsDMppm() {
        return absDMppm;
    }

    public void setAbsDMppm(double absDMppm) {
        this.absDMppm = absDMppm;
    }

    public double getIsoDM() {
        return isoDM;
    }

    public void setIsoDM(double isoDM) {
        this.isoDM = isoDM;
    }

    public double getIsoDMppm() {
        return isoDMppm;
    }

    public void setIsoDMppm(double isoDMppm) {
        this.isoDMppm = isoDMppm;
    }

    public int getMcNum() {
        return mcNum;
    }

    public void setMcNum(int mcNum) {
        this.mcNum = mcNum;
    }

    public double getVarmodsRatio() {
        return varmodsRatio;
    }

    public void setVarmodsRatio(double varmodsNum) {
        this.varmodsRatio = varmodsNum;
    }

    public int getVarmodsCount() {
        return varmodsCount;
    }

    public void setVarmodsCount(int varmodsCount) {
        this.varmodsCount = varmodsCount;
    }

    public int getModifiable() {
        return modifiable;
    }

    public void setModifiable(int modifiable) {
        this.modifiable = modifiable;
    }

    public void setModified(int modified) {
        this.modified = modified;
    }

    public int getModified() {
        return modified;
    }

    public double getTotIne() {
        return totIne;
    }

    public void setTotIne(double totIne) {
        this.totIne = totIne;
    }

    public double getIntMatchedTot() {
        return intMatchedTot;
    }

    public void setIntMatchedTot(double intMatchedTot) {
        this.intMatchedTot = intMatchedTot;
    }

    public double getRelIntMatchTot() {
        return relIntMatchTot;
    }

    public void setRelIntMatchTot(double relIntMatchTot) {
        this.relIntMatchTot = relIntMatchTot;
    }

    public double getFragDeltaMed() {
        return fragDeltaMed;
    }

    public void setFragDeltaMed(double fragDeltaMed) {
        this.fragDeltaMed = fragDeltaMed;
    }

    public double getFragDeltaIqr() {
        return fragDeltaIqr;
    }

    public void setFragDeltaIqr(double fragDeltaIqr) {
        this.fragDeltaIqr = fragDeltaIqr;
    }

    public double getFragDeltaMedPPM() {
        return fragDeltaMedPPM;
    }

    public void setFragDeltaMedPPM(double fragDeltaMedPPM) {
        this.fragDeltaMedPPM = fragDeltaMedPPM;
    }

    public double getFragDeltaIqrPPM() {
        return fragDeltaIqrPPM;
    }

    public void setFragDeltaIqrPPM(double fragDeltaIqrPPM) {
        this.fragDeltaIqrPPM = fragDeltaIqrPPM;
    }

    public double getFragDeltaPolyFit() {
        return fragDeltaPolyFit;
    }

    public void setFragDeltaPolyFit(double fragDeltaPolyFit) {
        this.fragDeltaPolyFit = fragDeltaPolyFit;
    }

    public int getLongestSeq() {
        return longestSeq;
    }

    public void setLongestSeq(int longestSeq) {
        this.longestSeq = longestSeq;
    }

    public int getQmatch() {
        return qmatch;
    }

    public void setQmatch(int qmatch) {
        this.qmatch = qmatch;
    }

    public String getPepString() {
        return pepString;
    }

    public void setPepString(String pepString) {
        this.pepString = pepString;
    }

    public String getProteins() {
        return proteins;
    }

    public void setProteins(String proList) {
        this.proteins = proList;
    }

    public void setMaxFragIonInt(double maxFragIonInt) {
        this.maxFragIonInt = maxFragIonInt;
    }

    public double getMaxFragIonInt() {
        return maxFragIonInt;
    }

    public void setSummedScore(double summedScore) {
        this.summedScore = summedScore;
    }

    public double getSummedScore() {
        return summedScore;
    }

    public void setPepEvalue(double pepEvalue) {
        this.pepEvalue = pepEvalue;
    }

    public double getPepEvalue() {
        return pepEvalue;
    }

    public void setHyperScore(double hyperScore) {
        this.hyperScore = hyperScore;
    }

    public double getHyperScore() {
        return hyperScore;
    }

    public void setNextScore(double nextScore) {
        this.nextScore = nextScore;
    }

    public double getNextScore() {
        return nextScore;
    }

    public void setDeltaScore(double deltaScore) {
        this.deltaScore = deltaScore;
    }

    public double getDeltaScore() {
        return deltaScore;
    }

    public void setFragError(double fragError) {
        this.fragError = fragError;
    }

    public double getFragError() {
        return fragError;
    }

    public void setbScore(double bScore) {
        this.bScore = bScore;
    }

    public double getbScore() {
        return bScore;
    }

    public void setyScore(double yScore) {
        this.yScore = yScore;
    }

    public double getyScore() {
        return yScore;
    }

    public void setBions(int bions) {
        this.bions = bions;
    }

    public int getBions() {
        return bions;
    }

    public void setYions(int yions) {
        this.yions = yions;
    }

    public int getYions() {
        return yions;
    }

    public void setbIonFraction(double bIonFraction) {
        this.bIonFraction = bIonFraction;
    }

    public double getbIonFraction() {
        return bIonFraction;
    }

    public void setyIonFraction(double yIonFraction) {
        this.yIonFraction = yIonFraction;
    }

    public double getyIonFraction() {
        return yIonFraction;
    }

    public void setcTermTryptic(boolean cTermTryptic) {
        this.cTermTryptic = cTermTryptic;
    }

    public boolean iscTermTryptic() {
        return cTermTryptic;
    }

    public void setPepLen(int pepLen) {
        this.pepLen = pepLen;
    }

    public int getPepLen() {
        return pepLen;
    }

    public void setCharge1(double charge1) {
        this.charge1 = charge1;
    }

    public double getCharge1() {
        return charge1;
    }

    public void setCharge2(double charge2) {
        this.charge2 = charge2;
    }

    public double getCharge2() {
        return charge2;
    }

    public void setCharge3(double charge3) {
        this.charge3 = charge3;
    }

    public double getCharge3() {
        return charge3;
    }

    public void setCharge4(double charge4) {
        this.charge4 = charge4;
    }

    public double getCharge4() {
        return charge4;
    }

    public void setMatchedB(double matchedB) {
        this.matchedB = matchedB;
    }

    public double getMatchedB() {
        return matchedB;
    }

    public void setMatchedY(double matchedY) {
        this.matchedY = matchedY;
    }

    public double getMatchedY() {
        return matchedY;
    }

    public void setMatchedX(double matchedX) {
        this.matchedX = matchedX;
    }

    public double getMatchedX() {
        return matchedX;
    }

    public void setMatchedA(double matchedA) {
        this.matchedA = matchedA;
    }

    public double getMatchedA() {
        return matchedA;
    }

    public void setMatchedC(double matchedC) {
        this.matchedC = matchedC;
    }

    public double getMatchedC() {
        return matchedC;
    }

    public void setMatchedZ(double matchedZ) {
        this.matchedZ = matchedZ;
    }

    public double getMatchedZ() {
        return matchedZ;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    public double getPvalue() {
        return pvalue;
    }

    public void setMod(String mod) {
        this.mod = mod;
    }

    public String getMod() {
        return mod;
    }

    /**
     * @param enzN the enzN to set
     */
    public void setEnzN(int enzN) {
        this.enzN = enzN;
    }

    /**
     * @return the enzN
     */
    public int getEnzN() {
        return enzN;
    }

    /**
     * @param sumSquareFragError the sumSquareFragError to set
     */
    public void setSumSquareFragError(double sumSquareFragError) {
        this.sumSquareFragError = sumSquareFragError;
    }

    /**
     * @return the sumSquareFragError
     */
    public double getSumSquareFragError() {
        return sumSquareFragError;
    }

    /**
     * @param stdevFragError the stdevFragError to set
     */
    public void setStdevFragError(double stdevFragError) {
        this.stdevFragError = stdevFragError;
    }

    /**
     * @return the stdevFragError
     */
    public double getStdevFragError() {
        return stdevFragError;
    }

    /**
     * @param nTermTryptic the nTermTryptic to set
     */
    public void setnTermTryptic(boolean nTermTryptic) {
        this.nTermTryptic = nTermTryptic;
    }

    /**
     * @return the nTermTryptic
     */
    public boolean isnTermTryptic() {
        return nTermTryptic;
    }

    public void setHitID(String hitID) {
        this.hitID = hitID;
    }

    public String getHitID() {
        return hitID;
    }

    public void setIonSeriesFrac(HashMap<String, Integer> ionSeriesFrac) {
        this.ionSeriesFrac = ionSeriesFrac;
    }

    public void setIonSeriesInt(HashMap<String, Double> ionSeriesInt) {
        this.ionSeriesInt = ionSeriesInt;
    }

    public HashMap<String, Integer> getIonSeriesFrac() {
        return ionSeriesFrac;
    }

    public HashMap<String, Double> getIonSeriesInt() {
        return ionSeriesInt;
    }

    public void setOutIonInte(HashMap<String, Integer> outIonInte) {
        this.outIonInte = outIonInte;
    }

    public HashMap<String, Integer> getOutIonInte() {
        return outIonInte;
    }
}
