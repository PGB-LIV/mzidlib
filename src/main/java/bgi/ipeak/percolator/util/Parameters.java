/**
 * 
 */
package bgi.ipeak.percolator.util;


/**
 * @author Administrator
 *
 */
public class Parameters {
	/**
	 * Parent mass tolerance, default in dalton.
	 */
	private double parentTolerance=0.0;
	/**
	 * If true, use ppm as the unit of parent mass tolerance.
	 */
	private boolean parentTolUnitPPM=false;
	/**
	 * Fragment ion mass tolerance, default in dalton.
	 */
	private double fragmentTolerance=0.0;
	/**
	 * If true, use ppm as the unit of fragment ion mass tolerance.
	 */
	private boolean fragmentTolUnitPPM=false;
	/**
	 * Enzyme name.
	 */
	private String enzymeName=null;
	/**
	 * Enzyme id, based on OMSSA.
	 */
	private String enzymeAccession=null;
	/**
	 * Mass type of parent mass tolerance.
	 */
	private String parentMassType=null;
	/**
	 * Mass type of fragment ion mass tolerance.
	 */
	private String fragmentMassType=null;
	/**
	 * Variable modifications.
	 */
	@SuppressWarnings("deprecation")
	private Modification[] mods=null;
	/**
	 * Is the search performed by automated decoy database search?
	 */
	private boolean autoDecoy=false;
	/**
	 * False discovery rate to filter the PSMs after percolator or search.
	 */
	private double falseDiscoveryRate=0;
	/**
	 * Is a ion series will be searched?
	 */
	private boolean aIon=false;
	/**
	 * Is b ion series will be searched?
	 */
	private boolean bIon=false;
	/**
	 * Is c ion series will be searched?
	 */
	private boolean cIon=false;
	/**
	 * Is x ion series will be searched?
	 */
	private boolean xIon=false;
	/**
	 * Is y ion series will be searched?
	 */
	private boolean yIon=false;
	/**
	 * Is z ion series will be searched?
	 */
	private boolean zIon=false;
	/**
	 * Empty constructor.
	 */
	public Parameters(){}
	/**
	 * @return the parentTolerance
	 */
	public double getParentTolerance() {
		return parentTolerance;
	}
	/**
	 * @return the parentTolUnitPPM
	 */
	public boolean isParentTolUnitPPM() {
		return parentTolUnitPPM;
	}
	/**
	 * @return the fragmentTolerance
	 */
	public double getFragmentTolerance() {
		return fragmentTolerance;
	}
	/**
	 * @return the fragmentTolUnitPPM
	 */
	public boolean isFragmentTolUnitPPM() {
		return fragmentTolUnitPPM;
	}
	/**
	 * @return the enzymeName
	 */
	public String getEnzymeName() {
		return enzymeName;
	}
	/**
	 * @return the parentMassType
	 */
	public String getParentMassType() {
		return parentMassType;
	}
	/**
	 * @return the fragmentMassType
	 */
	public String getFragmentMassType() {
		return fragmentMassType;
	}
	/**
	 * @return the autoDecoy
	 */
	public boolean isAutoDecoy() {
		return autoDecoy;
	}
	/**
	 * @return the falseDiscoveryRate
	 */
	public double getFalseDiscoveryRate() {
		return falseDiscoveryRate;
	}
	/**
	 * @return the aIon
	 */
	public boolean isaIon() {
		return aIon;
	}
	/**
	 * @return the bIon
	 */
	public boolean isbIon() {
		return bIon;
	}
	/**
	 * @return the cIon
	 */
	public boolean iscIon() {
		return cIon;
	}
	/**
	 * @return the xIon
	 */
	public boolean isxIon() {
		return xIon;
	}
	/**
	 * @return the yIon
	 */
	public boolean isyIon() {
		return yIon;
	}
	/**
	 * @return the zIon
	 */
	public boolean iszIon() {
		return zIon;
	}
	/**
	 * @param parentTolerance the parentTolerance to set
	 */
	public void setParentTolerance(double parentTolerance) {
		this.parentTolerance = parentTolerance;
	}
	/**
	 * @param parentTolUnitPPM the parentTolUnitPPM to set
	 */
	public void setParentTolUnitPPM(boolean parentTolUnitPPM) {
		this.parentTolUnitPPM = parentTolUnitPPM;
	}
	/**
	 * @param fragmentTolerance the fragmentTolerance to set
	 */
	public void setFragmentTolerance(double fragmentTolerance) {
		this.fragmentTolerance = fragmentTolerance;
	}
	/**
	 * @param fragmentTolUnitPPM the fragmentTolUnitPPM to set
	 */
	public void setFragmentTolUnitPPM(boolean fragmentTolUnitPPM) {
		this.fragmentTolUnitPPM = fragmentTolUnitPPM;
	}
	/**
	 * @param enzymeName the enzymeName to set
	 */
	public void setEnzymeName(String enzymeName) {
		this.enzymeName = enzymeName;
	}
	/**
	 * @param parentMassType the parentMassType to set
	 */
	public void setParentMassType(String parentMassType) {
		this.parentMassType = parentMassType;
	}
	/**
	 * @param fragmentMassType the fragmentMassType to set
	 */
	public void setFragmentMassType(String fragmentMassType) {
		this.fragmentMassType = fragmentMassType;
	}
	/**
	 * @param autoDecoy the autoDecoy to set
	 */
	public void setAutoDecoy(boolean autoDecoy) {
		this.autoDecoy = autoDecoy;
	}
	/**
	 * @param falseDiscoveryRate the falseDiscoveryRate to set
	 */
	public void setFalseDiscoveryRate(double falseDiscoveryRate) {
		this.falseDiscoveryRate = falseDiscoveryRate;
	}
	/**
	 * @param aIon the aIon to set
	 */
	public void setaIon(boolean aIon) {
		this.aIon = aIon;
	}
	/**
	 * @param bIon the bIon to set
	 */
	public void setbIon(boolean bIon) {
		this.bIon = bIon;
	}
	/**
	 * @param cIon the cIon to set
	 */
	public void setcIon(boolean cIon) {
		this.cIon = cIon;
	}
	/**
	 * @param xIon the xIon to set
	 */
	public void setxIon(boolean xIon) {
		this.xIon = xIon;
	}
	/**
	 * @param yIon the yIon to set
	 */
	public void setyIon(boolean yIon) {
		this.yIon = yIon;
	}
	/**
	 * @param zIon the zIon to set
	 */
	public void setzIon(boolean zIon) {
		this.zIon = zIon;
	}
	/**
	 * @param mods the mods to set
	 */
	@SuppressWarnings("deprecation")
	public void setMods(Modification[] mods) {
		this.mods = mods;
	}
	/**
	 * @return the mods
	 */
	@SuppressWarnings("deprecation")
	public Modification[] getMods() {
		return mods;
	}
	/**
	 * @param enzymeAccession the enzymeAccession to set
	 */
	public void setEnzymeAccession(String enzymeAccession) {
		this.enzymeAccession = enzymeAccession;
	}
	/**
	 * @return the enzymeAccession
	 */
	public String getEnzymeAccession() {
		return enzymeAccession;
	}
}
