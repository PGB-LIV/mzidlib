/**
 * 
 */
package bgi.ipeak.percolator.util;


/**
 * @author Administrator
 *
 */
public class FragmentIon implements Ion{
	
	/**
	 * The M/Z of the fragment ion.
	 */
	private double mz=0.0;
	/**
	 * The fragment ion intensity if it is a match.
	 */
	private double intensity=0.0;
	/**
	 * The soft type ID of the fragment ion. 
	 */
	private int ID=-1;
	/**
	 * The type of the fragment ion.
	 */
	private String type=null;
	/**
	 * The sub index of the fragment ion.
	 */
	private int number=-1;
	/**
	 * The mass error between the theoretical fragment ion and the matched peak in dalton.
	 */
	private double massError=0.0;
	/**
	 * If this fragment ion is double charged, this boolean will be set true.
	 */
	private boolean doubleCharged=false;

	public boolean setIDAccordeingToType(){
		if(this.type==null || this.type.isEmpty())
			return false;
		if(this.type.indexOf("A|a")!=-1)
			this.setID(this.isDoubleCharged()?A_DOUBLE_ION:A_ION);
		if(this.type.indexOf("B|b")!=-1)
			this.setID(this.isDoubleCharged()?B_DOUBLE_ION:B_ION);
		if(this.type.indexOf("C|c")!=-1)
			this.setID(this.isDoubleCharged()?C_DOUBLE_ION:C_ION);
		if(this.type.indexOf("X|x")!=-1)
			this.setID(this.isDoubleCharged()?X_DOUBLE_ION:X_ION);
		if(this.type.indexOf("Y|y")!=-1)
			this.setID(this.isDoubleCharged()?Y_DOUBLE_ION:Y_ION);
		if(this.type.indexOf("Z|z")!=-1)
			this.setID(this.isDoubleCharged()?Z_DOUBLE_ION:Z_ION);
		return true;
	}
	/**
	 * @param mz the mz to set
	 */
	public void setMz(double mz) {
		this.mz = mz;
	}

	/**
	 * @param intensity the intensity to set
	 */
	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}

	/**
	 * @param iD the iD to set
	 */
	public void setID(int iD) {
		ID = iD;
	}

	/**
	 * @param type the type to set
	 */
	public void setType(String type) {
		this.type = type;
	}

	/**
	 * @param number the number to set
	 */
	public void setNumber(int number) {
		this.number = number;
	}

	/**
	 * @param massError the massError to set
	 */
	public void setMassError(double massError) {
		this.massError = massError;
	}

	/**
	 * @param doubleCharged the doubleCharged to set
	 */
	public void setDoubleCharged(boolean doubleCharged) {
		this.doubleCharged = doubleCharged;
	}

	@Override
	public double getMZ() {
		return this.mz;
	}

	@Override
	public double getIntensity() {
		return this.intensity;
	}

	@Override
	public int getID() {
		return this.ID;
	}

	@Override
	public String getType() {
		return this.type;
	}

	@Override
	public int getNumber() {
		return this.number;
	}

	@Override
	public boolean isDoubleCharged() {
		return this.doubleCharged;
	}

	@Override
	public double getMassError() {
		return this.massError;
	}
}
