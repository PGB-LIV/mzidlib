package bgi.ipeak.percolator.util;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class SpectrumBasicInfo {
	/**
	 * The title of the MS/MS spectrum.
	 */
	private String title=null;
	/**
	 * The retention time of the MS/MS spectrum.
	 */
	private double retentionTime=0.0;
	/**
	 * The precursor m/z.
	 */
	private double precursorMZ=0.0;
	/**
	 * The precursor mass.
	 */
	private double precursorMass=0.0;
	/**
	 * The precursor charge.
	 */
	private byte charge=0;
	/**
	 * The precursor intensity.
	 */
	private double precursorIntensity=0.0;
	/**
	 * The spectrum index, 0-based coordinates.
	 */
	private int index=0;
	/**
	 * The Peak array that includes all the peaks found in the masspectrum.
	 */
	private Peak[] peakArray=null;
	
	/**
	 * Return the title of the MS/MS spectrum.
	 * @return the title
	 */
	public String getTitle() {
		return title;
	}
	/**
	 * Return the retention time of the MS/MS spectrum.
	 * @return the retentionTime
	 */
	public double getRetentionTime() {
		return retentionTime;
	}
	/**
	 * Return the precursor m/z of the MS/MS spectrum.
	 * @return the precursorMZ
	 */
	public double getPrecursorMZ() {
		return precursorMZ;
	}
	/**
	 * Return the precursor mass of the MS/MS spectrum.
	 * @return the precursorMass
	 */
	public double getPrecursorMass() {
		return precursorMass;
	}
	/**
	 * Return the precursor charge.
	 * @return the charge
	 */
	public byte getCharge() {
		return charge;
	}
	/**
	 * Return the precursor intensity.
	 * @return the precursorIntensity
	 */
	public double getPrecursorIntensity() {
		return precursorIntensity;
	}
	/**
	 * Return the index of the MS/MS spectrum, 0-based.
	 * @return the index
	 */
	public int getIndex() {
		return index;
	}
	/**
	 * Return the Peak arrray including all peaks found in the spectrum.
	 * @return the peakArray
	 */
	public Peak[] getPeakArray() {
		return peakArray;
	}
	/**
	 * Set the title of the MS/MS spectrum.
	 * @param title the title to set
	 */
	public void setTitle(String title) {
		this.title = title;
	}
	/**
	 * Set the retention time of the MS/MS spectrum.
	 * @param retentionTime the retentionTime to set
	 */
	public void setRetentionTime(double retentionTime) {
		this.retentionTime = retentionTime;
	}
	/**
	 * Set the precursor m/z of the MS/MS spectrum.
	 * @param precursorMZ the precursorMZ to set
	 */
	public void setPrecursorMZ(double precursorMZ) {
		this.precursorMZ = precursorMZ;
	}
	/**
	 * Set the precursor mass.
	 * @param precursorMass the precursorMass to set
	 */
	public void setPrecursorMass(double precursorMass) {
		this.precursorMass = precursorMass;
	}
	/**
	 * Set the precursor charge.
	 * @param charge the charge to set
	 */
	public void setCharge(byte charge) {
		this.charge = charge;
	}
	/**
	 * Set the precursor intensity.
	 * @param precursorIntensity the precursorIntensity to set
	 */
	public void setPrecursorIntensity(double precursorIntensity) {
		this.precursorIntensity = precursorIntensity;
	}
	/**
	 * Set the index of the MS/MS spectrum.
	 * @param index the index to set
	 */
	public void setIndex(int index) {
		this.index = index;
	}
	/**
	 * Set the Peak array that includes all peaks found in the spectrum.
	 * @param peakArray the peakArray to set
	 */
	public void setPeakArray(Peak[] peakArray) {
		this.peakArray = peakArray;
	}
	/**
     * This method returns all the Intensity values of the peaks in a double[].
     * @return double[] Peak intensity.
     */
    public double[] getIntensityArray() {
        double[] lIntensityArray = new double[peakArray.length];
        for (int i = 0; i < peakArray.length; i++) {
            lIntensityArray[i] = peakArray[i].getIntensity();
        }
        return lIntensityArray;
    }
    /**
     * This method returns all the MZ values of the peaks in a double[].
     * @return double[] Peak mass.
     */
    public double[] getMZArray() {
        double[] lMZArray = new double[peakArray.length];
        for (int i = 0; i < peakArray.length; i++) {
            lMZArray[i] = peakArray[i].getMass();
        }
        return lMZArray;
    }
    /**
     * This method calculates the log-transformed total intensity of a spectrum.
     * The spectrum is segmented to 100Da-fixed-size bins. The 20 most intense peaks
     * are picked from each bin. Then sum up them as the total intensity of the spectrum.
     * @return the log-transformed total intensity of a spectrum.
     */
    public double calLogTotInt(){
		HashMap<Double, Double> mzToIntensity=new HashMap<Double, Double>();
		for(int j=0;j<getMZArray().length;j++)
			mzToIntensity.put(getMZArray()[j], getIntensityArray()[j]);
		double[] mzArray=getMZArray();
		Arrays.sort(mzArray);
		//The intensity of peaks in each bin are pushed into a vector.
		Vector<Double> intensityMeta=new Vector<Double>(); 
		double minMZ=mzArray[0];
		double totalInt=0;
		for(int j=0;j<mzArray.length;j++)
			if(mzArray[j]<minMZ+100){
				intensityMeta.add(mzToIntensity.get(mzArray[j]));
			}else{
				minMZ+=100;
				if(intensityMeta.size()!=0){
					Double[] inteMeta=intensityMeta.toArray(new Double[0]);
					Arrays.sort(inteMeta);
					//Only the 20 most intense peak are taken into account.
					for(int i=inteMeta.length-1;i>=inteMeta.length-20 && i>=0;i--)
						totalInt+=inteMeta[i];
					intensityMeta.clear();
					j--;
				}
			}
		if(intensityMeta.size()!=0){
			Double[] inteMeta=intensityMeta.toArray(new Double[0]);
			Arrays.sort(inteMeta);
			for(int i=inteMeta.length-1;i>=inteMeta.length-20 && i>=0;i--)
				totalInt+=inteMeta[i];
			intensityMeta.clear();
		}
		if(totalInt>=Math.E)
			return Math.log(totalInt)/Math.log(Math.E);
		else
			return -1.0;
	}
}
