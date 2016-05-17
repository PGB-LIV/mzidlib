package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "mod_t", propOrder = {
    "specificity",
    "delta",
    "ignore",
    "altName",
    "xref",
    "miscNotes"
})
public class ModT {

    @XmlElement(required = true)
    protected List<SpecificityT> specificity;
    @XmlElement(required = true)
    protected CompositionT delta;
    @XmlElement(name = "Ignore")
    protected List<CompositionT> ignore;
    @XmlElement(name = "alt_name")
    protected List<String> altName;
    protected List<XrefT> xref;
    @XmlElement(name = "misc_notes")
    protected String miscNotes;
    @XmlAttribute(required = true)
    protected String title;
    @XmlAttribute(name = "full_name", required = true)
    protected String fullName;
    @XmlAttribute(name = "username_of_poster", required = true)
    protected String usernameOfPoster;
    @XmlAttribute(name = "group_of_poster")
    protected String groupOfPoster;
    @XmlAttribute(name = "date_time_posted", required = true)
    protected String dateTimePosted;
    @XmlAttribute(name = "date_time_modified", required = true)
    protected String dateTimeModified;
    @XmlAttribute
    protected Boolean approved;
    @XmlAttribute(name = "ex_code_name")
    protected String exCodeName;
    @XmlAttribute(name = "record_id")
    protected Long recordId;

    public List<SpecificityT> getSpecificity() {
        if (specificity == null) {
            specificity = new ArrayList<SpecificityT>();
        }
        return this.specificity;
    }

    public CompositionT getDelta() {
        return delta;
    }

    public void setDelta(CompositionT value) {
        this.delta = value;
    }

    public List<CompositionT> getIgnore() {
        if (ignore == null) {
            ignore = new ArrayList<CompositionT>();
        }
        return this.ignore;
    }

    public List<String> getAltName() {
        if (altName == null) {
            altName = new ArrayList<String>();
        }
        return this.altName;
    }

    public List<XrefT> getXref() {
        if (xref == null) {
            xref = new ArrayList<XrefT>();
        }
        return this.xref;
    }

    public String getMiscNotes() {
        return miscNotes;
    }

    public void setMiscNotes(String value) {
        this.miscNotes = value;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String value) {
        this.title = value;
    }

    public String getFullName() {
        return fullName;
    }

    public void setFullName(String value) {
        this.fullName = value;
    }

    public String getUsernameOfPoster() {
        return usernameOfPoster;
    }

    public void setUsernameOfPoster(String value) {
        this.usernameOfPoster = value;
    }

    public String getGroupOfPoster() {
        return groupOfPoster;
    }

    public void setGroupOfPoster(String value) {
        this.groupOfPoster = value;
    }

    public String getDateTimePosted() {
        return dateTimePosted;
    }

    public void setDateTimePosted(String value) {
        this.dateTimePosted = value;
    }

    public String getDateTimeModified() {
        return dateTimeModified;
    }

    public void setDateTimeModified(String value) {
        this.dateTimeModified = value;
    }

    public Boolean isApproved() {
        return approved;
    }

    public void setApproved(Boolean value) {
        this.approved = value;
    }

    public String getExCodeName() {
        return exCodeName;
    }

    public void setExCodeName(String value) {
        this.exCodeName = value;
    }

    public Long getRecordId() {
        return recordId;
    }

    public void setRecordId(Long value) {
        this.recordId = value;
    }

}
