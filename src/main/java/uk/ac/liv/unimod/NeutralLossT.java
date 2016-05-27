package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "NeutralLoss_t")
public class NeutralLossT
        extends CompositionT {

    @XmlAttribute
    protected Boolean flag;

    public boolean isFlag() {
        if (flag == null) {
            return false;
        } else {
            return flag;
        }
    }

    public void setFlag(Boolean value) {
        this.flag = value;
    }

}
