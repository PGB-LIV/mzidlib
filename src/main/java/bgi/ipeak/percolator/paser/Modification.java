package bgi.ipeak.percolator.paser;

import java.util.Vector;

public class Modification {

    /**
     * The location, -1 if n-term, peplength if c-term.
     */
    private int location = -2;
    /**
     * The residues. NULL if no definition.
     */
    private Vector<String> residues = null;
    /**
     * Compact residues. NULL if no definition.
     */
    private String compactResidues = null;
    /**
     * The monoisotopic mass.
     */
    private double monoMass = 0;
    /**
     * The average mass
     */
    private double aveMass = 0;
    /**
     * The name in OMSSA.
     */
    private String name = null;
    /**
     * The type in OMSSA.
     */
    private int modType = -1;
    /**
     * Is a fixed modification?
     */
    private boolean fixed = false;
    /**
     * Is a peptide n-term modification?
     */
    private boolean pepnterm = false;
    /**
     * Is a peptide c-term modification?
     */
    private boolean pepcterm = false;
    /**
     * Is a protein n-term modification?
     */
    private boolean pronterm = false;
    /**
     * Is a protein n-term modification?
     */
    private boolean procterm = false;

    /**
     * @return the location
     */
    public int getLocation() {
        return location;
    }

    /**
     * @return the residues
     */
    public Vector<String> getResidues() {
        return residues;
    }

    /**
     * @return the compactResidues
     */
    public String getCompactResidues() {
        return compactResidues;
    }

    /**
     * @return the monoMass
     */
    public double getMonoMass() {
        return monoMass;
    }

    /**
     * @return the aveMass
     */
    public double getAveMass() {
        return aveMass;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @return the modType
     */
    public int getModType() {
        return modType;
    }

    /**
     * @return the fixed
     */
    public boolean isFixed() {
        return fixed;
    }

    /**
     * @param location the location to set
     */
    public void setLocation(int location) {
        this.location = location;
    }

    /**
     * @param residues the residues to set
     */
    public void setResidues(Vector<String> residues) {
        this.residues = residues;
    }

    /**
     * @param compactResidues the compactResidues to set
     */
    public void setCompactResidues(String compactResidues) {
        this.compactResidues = compactResidues;
    }

    /**
     * @param monoMass the monoMass to set
     */
    public void setMonoMass(double monoMass) {
        this.monoMass = monoMass;
    }

    /**
     * @param aveMass the aveMass to set
     */
    public void setAveMass(double aveMass) {
        this.aveMass = aveMass;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @param modType the modType to set
     */
    public void setModType(int modType) {
        this.modType = modType;
    }

    /**
     * @param fixed the fixed to set
     */
    public void setFixed(boolean fixed) {
        this.fixed = fixed;
    }

    /**
     * @return the pepnterm
     */
    public boolean isPepnterm() {
        return pepnterm;
    }

    /**
     * @return the pepcterm
     */
    public boolean isPepcterm() {
        return pepcterm;
    }

    /**
     * @return the pronterm
     */
    public boolean isPronterm() {
        return pronterm;
    }

    /**
     * @return the procterm
     */
    public boolean isProcterm() {
        return procterm;
    }

    /**
     * @param pepnterm the pepnterm to set
     */
    public void setPepnterm(boolean pepnterm) {
        this.pepnterm = pepnterm;
    }

    /**
     * @param pepcterm the pepcterm to set
     */
    public void setPepcterm(boolean pepcterm) {
        this.pepcterm = pepcterm;
    }

    /**
     * @param pronterm the pronterm to set
     */
    public void setPronterm(boolean pronterm) {
        this.pronterm = pronterm;
    }

    /**
     * @param procterm the procterm to set
     */
    public void setProcterm(boolean procterm) {
        this.procterm = procterm;
    }
}
