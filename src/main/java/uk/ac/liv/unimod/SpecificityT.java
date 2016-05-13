

package uk.ac.liv.unimod;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "specificity_t", propOrder = {
    "neutralLoss",
    "pepNeutralLoss",
    "miscNotes"
})
public class SpecificityT {

    @XmlElement(name = "NeutralLoss")
    protected List<NeutralLossT> neutralLoss;
    @XmlElement(name = "PepNeutralLoss")
    protected List<PepNeutralLossT> pepNeutralLoss;
    @XmlElement(name = "misc_notes")
    protected String miscNotes;
    @XmlAttribute
    protected Boolean hidden;
    @XmlAttribute(required = true)
    protected String site;
    @XmlAttribute(required = true)
    protected PositionT position;
    @XmlAttribute(required = true)
    protected String classification;
    @XmlAttribute(name = "spec_group")
    protected BigInteger specGroup;

   
    public List<NeutralLossT> getNeutralLoss() {
        if (neutralLoss == null) {
            neutralLoss = new ArrayList<NeutralLossT>();
        }
        return this.neutralLoss;
    }


    public List<PepNeutralLossT> getPepNeutralLoss() {
        if (pepNeutralLoss == null) {
            pepNeutralLoss = new ArrayList<PepNeutralLossT>();
        }
        return this.pepNeutralLoss;
    }


    public String getMiscNotes() {
        return miscNotes;
    }

    public void setMiscNotes(String value) {
        this.miscNotes = value;
    }


    public boolean isHidden() {
        if (hidden == null) {
            return false;
        } else {
            return hidden;
        }
    }

  
    public void setHidden(Boolean value) {
        this.hidden = value;
    }


    public String getSite() {
        return site;
    }

 
    public void setSite(String value) {
        this.site = value;
    }

  
    public PositionT getPosition() {
        return position;
    }

 
    public void setPosition(PositionT value) {
        this.position = value;
    }


    public String getClassification() {
        return classification;
    }

 
    public void setClassification(String value) {
        this.classification = value;
    }

    
    public BigInteger getSpecGroup() {
        if (specGroup == null) {
            return new BigInteger("1");
        } else {
            return specGroup;
        }
    }

  
    public void setSpecGroup(BigInteger value) {
        this.specGroup = value;
    }

}
