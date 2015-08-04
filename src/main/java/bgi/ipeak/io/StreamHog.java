package bgi.ipeak.io;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class StreamHog extends Thread {
	InputStream is;
	FileWriter fo;
	
	public StreamHog(InputStream is, FileWriter fo) {
		this.is = is;
		this.fo = fo;
		start();
	}
	public void run()
	{
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line = null;
			while ( (line = br.readLine()) != null) {
				if(this.fo!=null){
					fo.write(line+"\n");
					fo.flush();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}