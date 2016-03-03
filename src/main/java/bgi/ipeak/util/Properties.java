/**
 *
 */
package bgi.ipeak.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


/**
 * @author Administrator
 *
 */
public class Properties {

    private String modFile_path = "";
    private String usermodFile_path = "";
    private String percolator_path = "";
//	private Integer system_type;
    private static Properties the_Properties = null;

    /*
     public  Properties(String percolator_path, String mod_path)  {
     this.modFile_path=mod_path;
     this.percolator_path=percolator_path;
		
     }
     */
    private Properties() {
        String jarPath = bgi.ipeak.IPeak.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        File jarPFile = new File(jarPath);
        jarPath = jarPFile.getParent();
        
        //set the mods.xml file path
        this.modFile_path = jarPath +File.separator+ "mods.xml";
        File modFile = new File(this.modFile_path);
        //if don't find the mods.xml file in *.jar folder, then we will use the default file
        if(!modFile.isFile()){
            System.out.println("Using the default mods.xml file\n");
            this.modFile_path = "mods.xml";
            InputStream inMods = ClassLoader.getSystemClassLoader().getResourceAsStream("mods.xml");
            extractFileFromJar(inMods, this.modFile_path);
        }
        
        //set the usermods.xml file path
        this.usermodFile_path = jarPath +File.separator+ "usermods.xml";
        File usermodFile = new File(this.usermodFile_path);
        //if don't find the mods.xml file in *.jar folder, then we will use the default file
        if(!usermodFile.isFile()){
            System.out.println("Using the default usermods.xml file\n");
            this.usermodFile_path = "usermods.xml";
            InputStream inMods = ClassLoader.getSystemClassLoader().getResourceAsStream("usermods.xml");
            extractFileFromJar(inMods, this.usermodFile_path);
        }
        
        
        //set percolator path
        jarPath = jarPFile.getParent();
        String OS = System.getProperty("os.name").toLowerCase();
        //OS=OS.replaceAll("\\\\", "/");
        //for linux
        if (OS.indexOf("nux") >= 0) {
            this.percolator_path = jarPath + File.separator+ "percolator"+File.separator+ "percolator_linux.exe";
            //for windows
        } else if (OS.indexOf("win") >= 0) {
            this.percolator_path = jarPath + File.separator+ "percolator"+File.separator+ "percolator_win.exe";
            //for mac
        } else if (OS.indexOf("mac") >= 0) {
            this.percolator_path = jarPath + File.separator+ "percolator"+File.separator+ "percolator_mac.exe";
        } else {
            System.err.println("Don't support current system!");
            System.exit(0);
        }
        System.out.println("Percolator path: " + this.percolator_path);
        Check_SoftwarePath();
        
    }

    private Boolean check_file(String file_path) {
        File theFile = new File(file_path);
        if (!theFile.exists() || !theFile.isFile()) {
            System.err.println("can not find " + file_path + "\n");
        }
        return (theFile.exists() && theFile.isFile());
    }

    public String getModFile() {
        return this.modFile_path;
    }
    
    public String getUserModFile() {
        return this.usermodFile_path;
    }

    /**
     * Auto get the path of percolator
     *
     * @return Percolator path.
     */
    public String getPercolator() {
        
        return this.percolator_path;
    }

    public static String getPercolator_path() {
        if(the_Properties==null){
            the_Properties=new Properties();
        }
        return the_Properties.getPercolator();
    }

    public static String getModFile_path() {
        if(the_Properties==null){
            the_Properties=new Properties();
        }
        return the_Properties.getModFile();
    }
    
    public static String getUserModFile_path() {
        if(the_Properties==null){
            the_Properties=new Properties();
        }
        return the_Properties.getUserModFile();
    }

    private void Check_SoftwarePath() {
        if (check_file(modFile_path) && check_file(percolator_path) && check_file(usermodFile_path)) {
            System.out.println(modFile_path + "\n" + usermodFile_path+"\n"+ percolator_path + "\n");
        } else {
            System.err.println("Can not find the percolator.exe, the mod file or the usermod file\n");
        }
    }
    
     /*
     * Helper method to get the mod file out of the jar, since local file is
     * needed for OMXParser
     * This function is from Omssa2mzid.java
     */
    private void extractFileFromJar(InputStream in, String filename) {
        try {

            StringBuilder builder = new StringBuilder();
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                builder.append(line);
            }
            br.close();
            String text = builder.toString();
            FileWriter writer = new FileWriter(filename);
            writer.write(text);
            writer.close();

        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }


    }
}
