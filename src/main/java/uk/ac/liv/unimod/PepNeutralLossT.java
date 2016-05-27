package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "PepNeutralLoss_t")
public class PepNeutralLossT
        extends CompositionT {

    @XmlAttribute
    protected Boolean required;

    public boolean isRequired() {
        if (required == null) {
            return false;
        } else {
            return required;
        }
    }

    public void setRequired(Boolean value) {
        this.required = value;
    }

}
