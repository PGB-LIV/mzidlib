/**
 * 
 */
package bgi.ipeak.gui;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * @author Guilin Li
 */
public class ExtFilter extends FileFilter {
	/**
	 * The extension of file.
	 */
	private String ext=null;
	/**
	 * Constructor for ExtFilter
	 * @param str The extension.
	 */
	public ExtFilter(String str){
		setExt(str);
	}
	
	@Override
	public boolean accept(File f) {
		if(f.isDirectory())
			return true;
		if(getExt()==null)
			return false;
		String extension=getExtension(f);
		if(extension==null)
			return false;
		String[] exts=getExt().split("/");
		for(String s : exts)
			if(extension.equals(s.toLowerCase()))
				return true;
		return false;
	}
	
	@Override
	public String getDescription() {
		String str="*."+getExt();
		str=str.replaceAll("/", " *.");
		return str;
	}
	
	private String getExtension(File f){
		String ext=null;
		String s=f.getName();
		int i=s.lastIndexOf(".");
		if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
	}
	/**
	 * @param ext the ext to set
	 */
	public void setExt(String ext) {
		this.ext = ext;
	}
	/**
	 * @return the ext
	 */
	public String getExt() {
		return ext;
	}

}
