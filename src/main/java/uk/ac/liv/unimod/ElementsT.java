

package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "elements_t", propOrder = {
    "elem"
})
public class ElementsT {

    protected List<ElemT> elem;


    public List<ElemT> getElem() {
        if (elem == null) {
            elem = new ArrayList<ElemT>();
        }
        return this.elem;
    }

}
