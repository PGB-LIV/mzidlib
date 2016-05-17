/**
 *
 */
package bgi.ipeak.percolator.util;

/**
 * @author Administrator
 *
 */
public interface Ion {

    /**
     * This variable is a hard coded identifier for a A-ion. Aion is the b-ion
     * minus CO.
     */
    public static final int A_ION = 0;
    /**
     * This variable is a hard coded identifier for a A++-ion. This is a double
     * charged a-ion.
     */
    public static final int A_DOUBLE_ION = 1;
    /**
     * This variable is a hard coded identifier for a B-ion.
     */
    public static final int B_ION = 2;
    /**
     * This variable is a hard coded identifier for a B++-ion. This is a double
     * charged b-ion.
     */
    public static final int B_DOUBLE_ION = 3;
    /**
     * This variable is a hard coded identifier for a C-ion. Cion is the b-ion
     * plus 'NH3'.
     */
    public static final int C_ION = 4;
    /**
     * This variable is a hard coded identifier for a C++-ion. This is a double
     * charged c-ion.
     */
    public static final int C_DOUBLE_ION = 5;
    /**
     * This variable is a hard coded identifier for a X-ion. Xion is the y-ion
     * plus 'CO', minus 'H2'.
     */
    public static final int X_ION = 6;
    /**
     * This variable is a hard coded identifier for a X++-ion. This is a double
     * charged x-ion.
     */
    public static final int X_DOUBLE_ION = 7;
    /**
     * This variable is a hard coded identifier for a Y-ion.
     */
    public static final int Y_ION = 8;
    /**
     * This variable is a hard coded identifier for a Y++-ion. This is a double
     * charged y-ion.
     */
    public static final int Y_DOUBLE_ION = 9;
    /**
     * This variable is a hard coded identifier for a Z-ion. Zion is the y-ion
     * minus 'NH3'.
     */
    public static final int Z_ION = 10;
    /**
     * This variable is a hard coded identifier for a Z++-ion. This is a double
     * charged z-ion.
     */
    public static final int Z_DOUBLE_ION = 11;

    /**
     * Return the M/Z of the fragment ion.
     *
     * @return double with the M/Z.
     */
    public double getMZ();

    /**
     * Return the fragment ion intensity if it is a match.
     *
     * @return the fragment ion intensity if it is a match.
     */
    public double getIntensity();

    /**
     * This method returns the soft type int ID of the fragment ion. e.g. an
     * A-ion has an integer ID = 0.
     * <br><b>A-ion : 0</B>
     * <br><b>A++-ion : 1</B>
     * <br><b>B-ion : 2</B>
     * <br><b>B++-ion : 3</B>
     * <br><b>C-ion : 4</B>
     * <br><b>C++-ion : 5</B>
     * <br><b>X-ion : 6</B>
     * <br><b>X++-ion : 7</B>
     * <br><b>Y-ion : 8</B>
     * <br><b>Y++-ion : 9</B>
     * <br><b>Z-ion : 10</B>
     * <br><b>Z++-ion : 11</B>
     *
     * @return int ID with the ion number.
     */
    public int getID();

    /**
     * Return the type of the fragment ion. (e.g. a b3 fragment ion will have
     * String iType = "b")
     *
     * @return String Return a String with the type of fragment ion.
     */
    public String getType();

    /**
     * Return the sub index of the fragment ion. (e.g. a b3 fragment ion will
     * have Number = 3)
     *
     * @return the sub index of the fragment ion. (e.g. a b3 fragment ion will
     * have Number = 3)
     */
    public int getNumber();

    /**
     * This method returns a boolean if this fragment ion is double charged.
     *
     * @return boolean return true if this fragment ion is double charged.
     */
    public boolean isDoubleCharged();

    /**
     * Return the mass error between the theoretical fragment ion and the
     * matched peak in dalton.
     * <b>The fragment ion must be a match in the mass spectrum before this
     * method can be used!</b>
     *
     * @return the mass error between the theoretical fragment ion and the
     * matched peak in dalton.
     */
    public double getMassError();
}
