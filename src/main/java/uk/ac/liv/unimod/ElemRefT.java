package uk.ac.liv.unimod;

import java.math.BigInteger;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "elem_ref_t")
public class ElemRefT {

    @XmlAttribute(required = true)
    protected String symbol;
    @XmlAttribute
    protected BigInteger number;

    public String getSymbol() {
        return symbol;
    }

    public void setSymbol(String value) {
        this.symbol = value;
    }

    public BigInteger getNumber() {
        if (number == null) {
            return new BigInteger("1");
        } else {
            return number;
        }
    }

    public void setNumber(BigInteger value) {
        this.number = value;
    }

}
