package uk.ac.liv.unimod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "xref_t", propOrder = {
    "text",
    "source",
    "url"
})
public class XrefT {

    @XmlElement(required = true)
    protected String text;
    @XmlElement(required = true)
    protected String source;
    protected String url;

    public String getText() {
        return text;
    }

    public void setText(String value) {
        this.text = value;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String value) {
        this.source = value;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String value) {
        this.url = value;
    }

}
