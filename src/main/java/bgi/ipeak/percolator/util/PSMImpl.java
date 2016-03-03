/**
 * 
 */
package bgi.ipeak.percolator.util;


/**
 * @author Administrator
 *
 */
public class PSMImpl implements PSMInf {
	/**
	 * The index of the PSM, the same as the index in MS/MS spectrum.
	 */
	private int PSMId=0;
	/**
	 * The matched peptide.
	 */
	private PeptideInfor peptide=null;
	/**
	 * The rank of this match.
	 */
	private int rank=0;
	/**
	 * True if the PSM map to the target database.
	 */
	private boolean target=false;
	/**
	 * The score of the PSM.
	 */
	private double score=0.0;
	/**
	 * The missed cleavages of the PSM.
	 */
	private int missedCleavages=0;
	/**
	 * The number of matched ions.
	 */
	private int numberOfIonsMatched=0;
	/**
	 * The performance parameters of MS/MS search.
	 */
	private Parameters performParameters=null;
	/**
	 * The matched ions.
	 */
	private FragmentIon[] fragmentIons=null;
	/**
	 * Empty constructor.
	 */
	public PSMImpl(){}
	/**
	 * Constructor for PSMImpl from a PeptideInfor instance and a FragmentIon instance array.
	 * @param aPeptide Peptide.
	 * @param aFragmentIons Fragment ions.
	 */
	public PSMImpl(PeptideInfor aPeptide, FragmentIon[] aFragmentIons) {
		setPeptide(aPeptide);
		setFragmentIons(aFragmentIons);
		setNumberOfIonsMatched(aFragmentIons.length);
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getPSMId()
	 */
	@Override
	public int getPSMId() {
		return this.PSMId;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getPeptide()
	 */
	@Override
	public PeptideInfor getPeptide() {
		return this.peptide;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getRank()
	 */
	@Override
	public int getRank() {
		return this.rank;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#isTarget()
	 */
	@Override
	public boolean isTarget() {
		return this.target;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getScore()
	 */
	@Override
	public double getScore() {
		return this.score;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getMissedCleavages()
	 */
	@Override
	public int getMissedCleavages() {
		return this.missedCleavages;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getNumberOfIonsMatched()
	 */
	@Override
	public int getNumberOfIonsMatched() {
		return this.numberOfIonsMatched;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getPerformParameters()
	 */
	@Override
	public Parameters getPerformParameters() {
		return this.performParameters;
	}

	/* (non-Javadoc)
	 * @see bgi.org.cn.ipeak.interfaces.PSMInf#getFragmentIons()
	 */
	@Override
	public FragmentIon[] getFragmentIons() {
		return this.fragmentIons;
	}

	/**
	 * @param pSMId the pSMId to set
	 */
	public void setPSMId(int pSMId) {
		PSMId = pSMId;
	}

	/**
	 * @param peptide the peptide to set
	 */
	public void setPeptide(PeptideInfor peptide) {
		this.peptide = peptide;
	}

	/**
	 * @param rank the rank to set
	 */
	public void setRank(int rank) {
		this.rank = rank;
	}

	/**
	 * @param target the target to set
	 */
	public void setTarget(boolean target) {
		this.target = target;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @param missedCleavages the missedCleavages to set
	 */
	public void setMissedCleavages(int missedCleavages) {
		this.missedCleavages = missedCleavages;
	}

	/**
	 * @param numberOfIonsMatched the numberOfIonsMatched to set
	 */
	public void setNumberOfIonsMatched(int numberOfIonsMatched) {
		this.numberOfIonsMatched = numberOfIonsMatched;
	}

	/**
	 * @param performParameters the performParameters to set
	 */
	public void setPerformParameters(Parameters performParameters) {
		this.performParameters = performParameters;
	}

	/**
	 * @param fragmentIons the fragmentIons to set
	 */
	public void setFragmentIons(FragmentIon[] fragmentIons) {
		this.fragmentIons = fragmentIons;
	}

}
