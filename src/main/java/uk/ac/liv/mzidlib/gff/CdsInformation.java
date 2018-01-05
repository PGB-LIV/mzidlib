package uk.ac.liv.mzidlib.gff;

public class CdsInformation {

    String seqid;
    String source;
    long start;
    long end;
    String strand;
    String phase;
    String attribute;

    public CdsInformation(String seqid, String source, long start, long end,
            String strand, String phase, String attr) {
        this.seqid = seqid;
        this.source = source;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.phase = phase;
        this.attribute = attr;
    }

    public String getSeqId() {
        return this.seqid;
    }

    public String getSource() {
        return this.source;
    }

    public long getStart() {
        return this.start;
    }

    public long getEnd() {
        return this.end;
    }

    public String getStrand() {
        return this.strand;
    }

    public String getSePhase() {
        return this.phase;
    }

    public String getAttribute() {
        return this.attribute;
    }

}
