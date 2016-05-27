package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "aa_t", propOrder = {
    "element"
})
public class AaT {

    protected List<ElemRefT> element;
    @XmlAttribute
    protected String title;
    @XmlAttribute(name = "three_letter")
    protected String threeLetter;
    @XmlAttribute(name = "full_name")
    protected String fullName;
    @XmlAttribute(name = "mono_mass")
    protected Double monoMass;
    @XmlAttribute(name = "avge_mass")
    protected Double avgeMass;

    public List<ElemRefT> getElement() {
        if (element == null) {
            element = new ArrayList<ElemRefT>();
        }
        return this.element;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String value) {
        this.title = value;
    }

    public String getThreeLetter() {
        return threeLetter;
    }

    public void setThreeLetter(String value) {
        this.threeLetter = value;
    }

    public String getFullName() {
        return fullName;
    }

    public void setFullName(String value) {
        this.fullName = value;
    }

    public Double getMonoMass() {
        return monoMass;
    }

    public void setMonoMass(Double value) {
        this.monoMass = value;
    }

    public Double getAvgeMass() {
        return avgeMass;
    }

    public void setAvgeMass(Double value) {
        this.avgeMass = value;
    }

}
