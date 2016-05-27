package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "modifications_t", propOrder = {
    "mod"
})
public class ModificationsT {

    protected List<ModT> mod;

    public List<ModT> getMod() {
        if (mod == null) {
            mod = new ArrayList<ModT>();
        }
        return this.mod;
    }

}
