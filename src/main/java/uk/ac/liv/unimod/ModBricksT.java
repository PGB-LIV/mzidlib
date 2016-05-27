package uk.ac.liv.unimod;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlType;

@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "mod_bricks_t", propOrder = {
    "brick"
})
public class ModBricksT {

    protected List<BrickT> brick;

    public List<BrickT> getBrick() {
        if (brick == null) {
            brick = new ArrayList<BrickT>();
        }
        return this.brick;
    }

}
