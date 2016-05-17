/**
 *
 */
package bgi.ipeak.percolator.util;

@Deprecated
/**
 * @author Administrator
 *
 */
public class Modification {

    /**
     * The standard name of the modification, based on UNIMOD.
     */
    public String type = null;
    /**
     * The mass of the modification.
     */
    public double mass = 0.0;
    /**
     * The modified residue.
     */
    public String residue = null;
    /**
     * The ID of the modification, based on UNIMOD.
     */
    public String modificationID = null;
    /**
     * The location of the modification. From N-term to C-term, 0-based
     * coordinates.
     */
    public int location = -1;
    /**
     * Is the modification fixed?
     */
    private boolean fixed = false;

    /**
     * Set the type of the modification
     *
     * @param type the type to set
     */
    public void setType(String type) {
        this.type = type;
    }

    /**
     * Set the mass of the modification.
     *
     * @param mass the mass to set
     */
    public void setMass(double mass) {
        this.mass = mass;
    }

    /**
     * Set the modified residue.
     *
     * @param residue the residue to set
     */
    public void setResidue(String residue) {
        this.residue = residue;
    }

    /**
     * Set the UNIMOD id of the modification.
     *
     * @param modificationID the modificationID to set
     */
    public void setModificationID(String modificationID) {
        this.modificationID = modificationID;
    }

    /**
     * Set the location of the modification.
     *
     * @param location the location to set
     */
    public void setLocation(int location) {
        this.location = location;
    }

    /**
     * Empty constructor.
     */
    public Modification() {
    }

    /**
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * @return the mass
     */
    public double getMass() {
        return mass;
    }

    /**
     * @return the residue
     */
    public String getResidue() {
        return residue;
    }

    /**
     * @return the modificationID
     */
    public String getModificationID() {
        return modificationID;
    }

    /**
     * @return the location
     */
    public int getLocation() {
        return location;
    }

    /**
     * @return the fixed
     */
    public boolean isFixed() {
        return fixed;
    }

    /**
     * @param fixed the fixed to set
     */
    public void setFixed(boolean fixed) {
        this.fixed = fixed;
    }
}
