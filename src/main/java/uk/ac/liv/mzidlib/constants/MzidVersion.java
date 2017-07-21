/*
 * Date: 21-Jul-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.constants.MzidVersion.java
 *
 */

package uk.ac.liv.mzidlib.constants;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 21-Jul-2017 12:35:02
 */
public enum MzidVersion {

    Version1_1("1.1"),
    Version1_2("1.2");

    private final String version;

    private MzidVersion(String ver) {
        this.version = ver;
    }

    @Override
    public String toString() {
        return "Mzid version:" + this.version;
    }

    public String getVersionString() {
        return this.version;
    }

    /**
     * Get MzidVersion object from input string.
     * Valid input string are "1.1" and "1.2",
     * other input string will result in NULL.
     *
     * @param ver version string
     *
     * @return MzidVersion
     */
    public static MzidVersion getVersion(String ver) {
        if (ver.equals("1.1")) {
            return Version1_1;
        } else if (ver.equals("1.2")) {
            return Version1_2;
        } else {
            return null;
        }
    }

}
