package uk.ac.liv.mzidlib.gff;

/**
 * Create this object from the result files obtained after running the pipeline.
 * We can create such objects to feed into the GFF routines to retrieve the
 * relevant information for mapping.
 *
 * @author riteshk
 *
 */
public class ProteinResults {

    private String accession;
    private boolean decoyOrNot;
    private long start;
    private long end;
    private String peptideEvidenceID;
    private String fdrScore;
    private String peptideSeq;

    public ProteinResults(String accession, boolean decoyOrNot, long start, long end, String peptideEvidenceID, String fdrScore, String peptideSeq) {
        this.accession = accession;
        this.decoyOrNot = decoyOrNot;
        this.start = start;
        this.end = end;
        this.peptideEvidenceID = peptideEvidenceID;
        this.peptideSeq = peptideSeq;
        this.fdrScore = fdrScore;
    }

    public String getAccession() {
        return this.accession;
    }

    public String getPeptideSeq() {
        return this.peptideSeq;
    }

    public String getfdrScore() {
        return this.fdrScore;
    }

    public boolean getDecoyOrNot() {
        return this.decoyOrNot;
    }

    public long getstart() {
        return this.start;
    }

    public long getEnd() {
        return this.end;
    }

    public String getPeptideEvidenceID() {
        return peptideEvidenceID;
    }
}
