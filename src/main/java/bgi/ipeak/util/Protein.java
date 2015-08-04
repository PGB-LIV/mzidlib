package bgi.ipeak.util;

import java.util.ArrayList;

public class Protein {
	/**
	 * Accession of the protein.
	 */
	private String accession=null;
	/**
	 * The peptide of the protein.
	 */
	private ArrayList<String> peps=null;
	/**
	 * The cluster ID of the protein.
	 */
	private int clusterID=-1;

	private String groupID;
	/**
	 * The protein accessions of the cluster or the same set. dual purpose.
	 */	
	private ArrayList<String> sameset=null;
	/**
	 * If a protein contains a unique peptide, it makes sense.
	 */
	private boolean nonsense=true;
	/**
	 * The description of the protein.
	 */
	private String description=null;
	/**
	 * The length of the protein.
	 */
	private int length=0;
	/**
	 * The mass of the protein.
	 */
	private double mass=-1;

	private String protein_sequence="";
	
	private String behalf_protein;
	
	public Protein(){
		setPeps(new ArrayList<String>());
		setSameset(new ArrayList<String>());
		setClusterID(-1);
		setNonsense(true);
	}
	
	public String getAccession() {
		return accession;
	}
	public void setAccession(String accession) {
		this.accession = accession;
	}
	public ArrayList<String> getPeps() {
		return peps;
	}
	public void setPeps(ArrayList<String> peps) {
		this.peps = peps;
	}
	public int getClusterID() {
		return clusterID;
	}
	public void setClusterID(int clusterID) {
		this.clusterID = clusterID;
	}
	public ArrayList<String> getSameset() {
		return sameset;
	}
	public void setSameset(ArrayList<String> sameset) {
		this.sameset = sameset;
	}
	public boolean isNonsense() {
		return nonsense;
	}
	public void setNonsense(boolean nonsense) {
		this.nonsense = nonsense;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public String getDescription() {
		return description;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public int getLength() {
		return length;
	}
	public void setMass(double mass) {
		this.mass = mass;
	}

	public double getMass() {
		return mass;
	}
	public String getProtein_sequence() {
		return protein_sequence;
	}
	public void setProtein_sequence(String protein_sequence) {
		this.protein_sequence = protein_sequence;
	}	
	public String getGroupID() {
		return groupID;
	}
	public String getBehalf_protein() {		
		return behalf_protein;
	}
	public void setBehalf_protein(String behalf_protein) {
		this.behalf_protein = behalf_protein;
	}
	public void setGroupID(String groupID) {
		this.groupID = groupID;
	}
}
