package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "elem_t")
public class ElemT {

    @XmlAttribute(required = true)
    protected String title;
    @XmlAttribute(name = "full_name", required = true)
    protected String fullName;
    @XmlAttribute(name = "avge_mass", required = true)
    protected double avgeMass;
    @XmlAttribute(name = "mono_mass", required = true)
    protected double monoMass;

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

    public double getAvgeMass() {
        return avgeMass;
    }

    public void setAvgeMass(double value) {
        this.avgeMass = value;
    }

    public double getMonoMass() {
        return monoMass;
    }

    public void setMonoMass(double value) {
        this.monoMass = value;
    }

}
