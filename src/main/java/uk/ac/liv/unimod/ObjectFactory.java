package uk.ac.liv.unimod;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.annotation.XmlElementDecl;
import javax.xml.bind.annotation.XmlRegistry;
import javax.xml.namespace.QName;

@XmlRegistry
public class ObjectFactory {

    private final static QName _Unimod_QNAME = new QName("http://www.unimod.org/xmlns/schema/unimod_2", "unimod");

    public ObjectFactory() {
    }

    public XrefT createXrefT() {
        return new XrefT();
    }

    public BrickT createBrickT() {
        return new BrickT();
    }

    public NeutralLossT createNeutralLossT() {
        return new NeutralLossT();
    }

    public ElementsT createElementsT() {
        return new ElementsT();
    }

    public ElemT createElemT() {
        return new ElemT();
    }

    public CompositionT createCompositionT() {
        return new CompositionT();
    }

    public SpecificityT createSpecificityT() {
        return new SpecificityT();
    }

    public PepNeutralLossT createPepNeutralLossT() {
        return new PepNeutralLossT();
    }

    public ModT createModT() {
        return new ModT();
    }

    public ModBricksT createModBricksT() {
        return new ModBricksT();
    }

    public AminoAcidsT createAminoAcidsT() {
        return new AminoAcidsT();
    }

    public ElemRefT createElemRefT() {
        return new ElemRefT();
    }

    public UnimodT createUnimodT() {
        return new UnimodT();
    }

    public AaT createAaT() {
        return new AaT();
    }

    public ModificationsT createModificationsT() {
        return new ModificationsT();
    }

    @XmlElementDecl(namespace = "http://www.unimod.org/xmlns/schema/unimod_2", name = "unimod")
    public JAXBElement<UnimodT> createUnimod(UnimodT value) {
        return new JAXBElement<UnimodT>(_Unimod_QNAME, UnimodT.class, null, value);
    }

}
