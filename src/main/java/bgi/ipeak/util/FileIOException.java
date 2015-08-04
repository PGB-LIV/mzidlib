/**
 * 
 */
package bgi.ipeak.util;

import java.io.IOException;

/**
 * @author Administrator
 *
 */
public class FileIOException extends IOException {
	private static final long serialVersionUID = -1L;
	private IOException e;

	public FileIOException(){}

	public FileIOException(String message)
	{
		super(message);
	}

	public FileIOException(String message, IOException e) {
		super(message);
	    this.e = e;
	}

	public IOException getNestedIOException() {
		return this.e;
	}
}
