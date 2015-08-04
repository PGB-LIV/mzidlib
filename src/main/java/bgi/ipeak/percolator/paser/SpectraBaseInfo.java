package bgi.ipeak.percolator.paser;
/**
 * The basic information of the spectrum.
 * @author Guilin Li.
 */
public class SpectraBaseInfo {
	/**
	 * The M/Z array.
	 */
	private double[] mz=null;
	/**
	 * The intensity array.
	 */
	private double[] intensity=null;
	/**
	 * The retention time.
	 */
	private double rt=0;
	/**
	 * The charge.
	 */
	private byte charge=0;
	
	public double[] getMz() {
		return mz;
	}
	public void setMz(double[] mz) {
		this.mz = mz;
	}
	public double[] getIntensity() {
		return intensity;
	}
	public void setIntensity(double[] intensity) {
		this.intensity = intensity;
	}
	public double getRt() {
		return rt;
	}
	public void setRt(double rt) {
		this.rt = rt;
	}
	public void setCharge(byte charge) {
		this.charge = charge;
	}
	public byte getCharge() {
		return charge;
	}
}
