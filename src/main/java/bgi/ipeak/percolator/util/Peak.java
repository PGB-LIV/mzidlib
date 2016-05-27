package bgi.ipeak.percolator.util;

public class Peak {

    /**
     * The mass of the peak.
     */
    private double mass = 0;
    /**
     * The intensity of the peak.
     */
    private double intensity = 0;
    /**
     * The peak index in a spectrum, 0-based.
     */
    private int index = 0;

    /**
     * Constructor that makes a peak.
     *
     * @param aMass The mass of the peak.
     * @param aIntensity The intensity of the peak.
     * @param aIndex	The index of the peak.
     */
    public Peak(double aMass, double aIntensity, int aIndex) {
        mass = aMass;
        intensity = aIntensity;
        index = aIndex;
    }

    public Peak() {
    }

    /**
     * This method returns a String of the mass and intensity of a Peak
     * instance.
     *
     * @return String with peak data
     */
    public String toString() {
        return "mass:" + mass + "\nintensity:" + intensity;
    }

    /**
     * Return the mass of the peak.
     *
     * @return the mass
     */
    public double getMass() {
        return mass;
    }

    /**
     * Return the intensity of the peak.
     *
     * @return the intensity
     */
    public double getIntensity() {
        return intensity;
    }

    /**
     * Return the index of the peak.
     *
     * @return the index
     */
    public int getIndex() {
        return index;
    }

    /**
     * Set the mass of the peak.
     *
     * @param mass the mass to set
     */
    public void setMass(double mass) {
        this.mass = mass;
    }

    /**
     * Set the intensity of the peak.
     *
     * @param intensity the intensity to set
     */
    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }

    /**
     * Set the index of the peak.
     *
     * @param index the index to set
     */
    public void setIndex(int index) {
        this.index = index;
    }
}
