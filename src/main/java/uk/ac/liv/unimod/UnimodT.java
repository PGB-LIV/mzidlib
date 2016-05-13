

package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlSchemaType;
import javax.xml.bind.annotation.XmlType;


@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "unimod_t", propOrder = {
    "elements",
    "modifications",
    "aminoAcids",
    "modBricks"
})
public class UnimodT {

    protected ElementsT elements;
    protected ModificationsT modifications;
    @XmlElement(name = "amino_acids")
    protected AminoAcidsT aminoAcids;
    @XmlElement(name = "mod_bricks")
    protected ModBricksT modBricks;
    @XmlAttribute(required = true)
    @XmlSchemaType(name = "unsignedShort")
    protected int majorVersion;
    @XmlAttribute(required = true)
    protected int minorVersion;

  
    public ElementsT getElements() {
        return elements;
    }

 
    public void setElements(ElementsT value) {
        this.elements = value;
    }

  
    public ModificationsT getModifications() {
        return modifications;
    }


    public void setModifications(ModificationsT value) {
        this.modifications = value;
    }

 
    public AminoAcidsT getAminoAcids() {
        return aminoAcids;
    }


    public void setAminoAcids(AminoAcidsT value) {
        this.aminoAcids = value;
    }

    
    public ModBricksT getModBricks() {
        return modBricks;
    }

 
    public void setModBricks(ModBricksT value) {
        this.modBricks = value;
    }

   
    public int getMajorVersion() {
        return majorVersion;
    }

 
    public void setMajorVersion(int value) {
        this.majorVersion = value;
    }

    
    public int getMinorVersion() {
        return minorVersion;
    }

  
    public void setMinorVersion(int value) {
        this.minorVersion = value;
    }

}
