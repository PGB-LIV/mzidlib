package bgi.ipeak.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Vector;

public class ProteinSummaryInformation {

    private String protein_group_id;
    private String protein_acc;
    private String reference;
    private Vector<String> uniq_peptides = new Vector<String>();
    private Vector<String> razor_peptides = new Vector<String>();

    private Vector<String> sameset_proteins;
    private Integer peptide_number;
    private Integer uniq_pep_number;
    private Integer psm_number;
    private Integer uniq_psm_number;
    private Double mass;

    public ProteinSummaryInformation() {

    }

    public void print_protein(BufferedWriter writed_file) throws IOException {
        writed_file.write(protein_group_id + "\t" + protein_acc + "\t" + mass + "\t" + peptide_number + "\t" + uniq_pep_number + "\t"
                + psm_number + "\t" + uniq_psm_number + "\t");

        if (uniq_peptides.isEmpty()) {
            writed_file.write("-\t");
        } else {
            for (int i = 0; i < uniq_peptides.size() - 1; i++) {
                writed_file.write(uniq_peptides.get(i) + ";");
            }
            writed_file.write(uniq_peptides.get(uniq_peptides.size() - 1) + "\t");
        }
        if (razor_peptides.isEmpty()) {
            writed_file.write("-\t");
        } else {
            for (int i = 0; i < razor_peptides.size() - 1; i++) {
                writed_file.write(razor_peptides.get(i) + ";");
            }
            writed_file.write(razor_peptides.get(razor_peptides.size() - 1) + "\t");
        }

        if (sameset_proteins.isEmpty()) {
            writed_file.write("-\t");
        } else {
            for (int i = 0; i < sameset_proteins.size() - 1; i++) {
                writed_file.write(sameset_proteins.get(i) + ";");
            }
            writed_file.write(sameset_proteins.get(sameset_proteins.size() - 1) + "\t");
        }
        writed_file.write(reference);
        writed_file.newLine();
        writed_file.flush();

    }

    public void setMass(Double mass) {
        this.mass = mass;
    }

    public void setPeptide_number(Integer peptide_number) {
        this.peptide_number = peptide_number;
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {

        return super.clone();
    }

    public void setProtein_group_id(String protein_group_id) {
        this.protein_group_id = protein_group_id;
    }

    public void setPsm_number(Integer psm_number) {
        this.psm_number = psm_number;
    }

    public void setRazor_peptides(Vector<String> razor_peptides) {
        this.razor_peptides = razor_peptides;
    }

    public void setReference(String reference) {
        this.reference = reference;
    }

    public void setSameset_proteins(Vector<String> sameset_proteins) {
        this.sameset_proteins = sameset_proteins;
    }

    public void setUniq_pep_number(Integer uniq_pep_number) {
        this.uniq_pep_number = uniq_pep_number;
    }

    public void setUniq_peptides(Vector<String> uniq_peptides) {
        this.uniq_peptides = uniq_peptides;
    }

    public void setUniq_psm_number(Integer uniq_psm_number) {
        this.uniq_psm_number = uniq_psm_number;
    }

    public void setProtein_acc(String protein_acc) {
        this.protein_acc = protein_acc;
    }
}
