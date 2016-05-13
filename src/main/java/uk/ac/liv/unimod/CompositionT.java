


package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlSeeAlso;
import javax.xml.bind.annotation.XmlType;



@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "composition_t", propOrder = {
    "element"
})
@XmlSeeAlso({
    NeutralLossT.class,
    PepNeutralLossT.class
})
public class CompositionT {

    protected List<ElemRefT> element;
    @XmlAttribute(required = true)
    protected String composition;
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

    public String getComposition() {
        return composition;
    }

    public void setComposition(String value) {
        this.composition = value;
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
