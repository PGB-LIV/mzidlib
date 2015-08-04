package uk.ac.liv.mzidlib.fasta;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;


/**
 *
 * @author Fawaz Ghali
 */
public class CombineFastaFiles {
    
    private String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    public void combine(String input, String output) {
     
        try {
            if (input != null && !input.trim().equals("") && output != null && !output.trim().equals("")) {
            	PrintWriter printWriter = new PrintWriter(new FileOutputStream(output));
                String[] files =input.split(";");

                for (int i = 0; i < files.length; i++) {
                    
                    System.out.println("Reading: "+ files[i]);

                    BufferedReader bufferedReader = new BufferedReader(new FileReader(files[i]));
                    String line = bufferedReader.readLine();
                    while (line != null) {
                        if(line.startsWith(">")){
                            
                            String newLine = line.substring(0,9)+alphabet.charAt(i)+"_"+line.substring(9);
                            printWriter.println(newLine);
                        }else{
                            printWriter.println(line);
                        }

                        
                        line = bufferedReader.readLine();
                    }
                    bufferedReader.close();
                }
                printWriter.close();
                System.out.println("Done.");
            }

        } catch (IOException ex) {
            ex.printStackTrace();
        } 
    }

    public static void main(String args[]) {
        
        CombineFastaFiles combineFastaFiles = new CombineFastaFiles();
        combineFastaFiles.combine(args[0], args[1]);

    }

}
