/**
 * 
 */
package bgi.ipeak.percolator.util;

import java.util.ArrayList;

/**
 * @author Administrator
 *
 */
public class PeptideInfor {
	/**
	 * Does the peptide map to decoy database?
	 */
	private boolean decoy=false;
	/**
	 * PeptideInfor sequence.
	 */
	private String peptideSequence=null;
	/**
	 * The theoretical mass of the peptide.
	 */
	private double theoreticalMass=0.0;
	/**
	 * The experimental mass of the peptide.
	 */
	private double experimentalMass=0.0;
	/**
	 * The charge of the peptide.
	 */
	private byte charge=0;
	/**
	 * The variable modifications of the peptide.
	 */
	@SuppressWarnings("deprecation")
	private Modification[] modList=null;
	/**
	 * The modified peptide sequence.
	 */
	private String modifiedSequence=null;
	/**
	 * The proteins accession mapped by the peptide.
	 */
	private ArrayList<String> proteins=null;
	/**
	 * The length of the peptide.
	 */
	private int peptideLength=0;
	/**
	 * Is the peptide is unique for the protein(s)?
	 */
	private boolean unique=false;
	/**
	 * The number of spectra matched to the peptide.
	 */
	private int qMatch=0;
	/**
	 * The amino acids before the peptide.
	 */
	private String postAA=null;
	/**
	 * The amino acids after the peptide.
	 */
	private String preAA=null;
	/**
     * The identity threshold, when the peptideHitScore is above this threshold it will tagged as an identification.
     */
    private double identityThreshold = 0;
    /**
     * The homology threshold.
     */
    private double homologyThreshold = 0;
    /**
     * @deprecated
     * The start position of the peptide in the mapped protein, 1-based.
     */
    private int proStart=0;
    /**
     * @deprecated
     * The end position of the peptide in the mapped protein, 1-based.
     */
    private int proEnd=0;
    /**
     * Empty constructor.
     */
    public PeptideInfor(){}
	/**
	 * @return the peptideSequence
	 */
	public String getPeptideSequence() {
		return peptideSequence;
	}
	/**
	 * @return the theoriticalMass
	 */
	public double getTheoriticalMass() {
		return theoreticalMass;
	}
	/**
	 * @return the experimentalMass
	 */
	public double getExperimentalMass() {
		return experimentalMass;
	}
	/**
	 * @return the charge
	 */
	public byte getCharge() {
		return charge;
	}
	/**
	 * @return the modifiedSequence
	 */
	public String getModifiedSequence() {
		return modifiedSequence;
	}
	/**
	 * @return the proteins
	 */
	public ArrayList<String> getProteins() {
		return proteins;
	}
	/**
	 * @return the peptideLength
	 */
	public int getPeptideLength() {
		return peptideLength;
	}
	/**
	 * @return the unique
	 */
	public boolean isUnique() {
		return unique;
	}
	/**
	 * @return the qMatch
	 */
	public int getqMatch() {
		return qMatch;
	}
	/**
	 * @return the postAA
	 */
	public String getPostAA() {
		return postAA;
	}
	/**
	 * @return the preAA
	 */
	public String getPreAA() {
		return preAA;
	}
	/**
	 * @return the identityThreshold
	 */
	public double getIdentityThreshold() {
		return identityThreshold;
	}
	/**
	 * @return the homologyThreshold
	 */
	public double getHomologyThreshold() {
		return homologyThreshold;
	}
	/**
	 * @return the proStart
	 */
	public int getProStart() {
		return proStart;
	}
	/**
	 * @return the proEnd
	 */
	public int getProEnd() {
		return proEnd;
	}
	/**
	 * @param peptideSequence the peptideSequence to set
	 */
	public void setPeptideSequence(String peptideSequence) {
		this.peptideSequence = peptideSequence;
	}
	/**
	 * @param theoreticalMass the theoriticalMass to set
	 */
	public void setTheoriticalMass(double theoreticalMass) {
		this.theoreticalMass = theoreticalMass;
	}
	/**
	 * @param experimentalMass the experimentalMass to set
	 */
	public void setExperimentalMass(double experimentalMass) {
		this.experimentalMass = experimentalMass;
	}
	/**
	 * @param charge the charge to set
	 */
	public void setCharge(byte charge) {
		this.charge = charge;
	}
	/**
	 * @param modifiedSequence the modifiedSequence to set
	 */
	public void setModifiedSequence(String modifiedSequence) {
		this.modifiedSequence = modifiedSequence;
	}
	/**
	 * @param proteins the proteins to set
	 */
	public void setProteins(ArrayList<String> proteins) {
		this.proteins = proteins;
	}
	/**
	 * @param peptideLength the peptideLength to set
	 */
	public void setPeptideLength(int peptideLength) {
		this.peptideLength = peptideLength;
	}
	/**
	 * @param unique the unique to set
	 */
	public void setUnique(boolean unique) {
		this.unique = unique;
	}
	/**
	 * @param qMatch the qMatch to set
	 */
	public void setqMatch(int qMatch) {
		this.qMatch = qMatch;
	}
	/**
	 * @param postAA the postAA to set
	 */
	public void setPostAA(String postAA) {
		this.postAA = postAA;
	}
	/**
	 * @param preAA the preAA to set
	 */
	public void setPreAA(String preAA) {
		this.preAA = preAA;
	}
	/**
	 * @param identityThreshold the identityThreshold to set
	 */
	public void setIdentityThreshold(double identityThreshold) {
		this.identityThreshold = identityThreshold;
	}
	/**
	 * @param homologyThreshold the homologyThreshold to set
	 */
	public void setHomologyThreshold(double homologyThreshold) {
		this.homologyThreshold = homologyThreshold;
	}
	/**
	 * @param proStart the proStart to set
	 */
	public void setProStart(int proStart) {
		this.proStart = proStart;
	}
	/**
	 * @param proEnd the proEnd to set
	 */
	public void setProEnd(int proEnd) {
		this.proEnd = proEnd;
	}
	/**
	 * @param modList the modList to set
	 */
	@SuppressWarnings("deprecation")
	public void setModList(Modification[] modList) {
		this.modList = modList;
	}
	/**
	 * @return the modList
	 */
	@SuppressWarnings("deprecation")
	public Modification[] getModList() {
		return modList;
	}
	/**
	 * @param decoy the decoy to set
	 */
	public void setDecoy(boolean decoy) {
		this.decoy = decoy;
	}
	/**
	 * @return the decoy
	 */
	public boolean isDecoy() {
		return decoy;
	}
}
