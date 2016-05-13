


package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlEnum;
import javax.xml.bind.annotation.XmlEnumValue;
import javax.xml.bind.annotation.XmlType;



@XmlType(name = "position_t")
@XmlEnum
public enum PositionT {

    @XmlEnumValue("Anywhere")
    ANYWHERE("Anywhere"),
    @XmlEnumValue("Any N-term")
    ANY_N_TERM("Any N-term"),
    @XmlEnumValue("Any C-term")
    ANY_C_TERM("Any C-term"),
    @XmlEnumValue("Protein N-term")
    PROTEIN_N_TERM("Protein N-term"),
    @XmlEnumValue("Protein C-term")
    PROTEIN_C_TERM("Protein C-term");
    private final String value;

    PositionT(String v) {
        value = v;
    }

    public String value() {
        return value;
    }

    public static PositionT fromValue(String v) {
        for (PositionT c: PositionT.values()) {
            if (c.value.equals(v)) {
                return c;
            }
        }
        throw new IllegalArgumentException(v);
    }

}
