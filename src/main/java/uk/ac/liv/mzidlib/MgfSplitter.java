package uk.ac.liv.mzidlib;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author fghali
 */
public class MgfSplitter {

    private int limit = 1000000;
    private String endTag = "END IONS";
    // C:\Users\fghali\Desktop\livepool\psi-pi\examples\1_1examples\55merge.mgf

    public static void main(String args[]) {
        MgfSplitter mgfSplitter = new MgfSplitter();
        mgfSplitter.SplitMgf(args);
    }

    public void SplitMgf(String[] args) {
        try {
            if (args == null || args.length != 2) {
                System.err.println("\n\nUsage:\n\tReverseFASTADB <input_database_file>\n");
            }
            int totalIons = 0;
            BufferedReader br = new BufferedReader(new FileReader(args[0]));
            try {
                String line = br.readLine();
                while (line != null) {
                    line = br.readLine();
                    if (line != null && line.equals("END IONS")) {
                        totalIons += 1;
                    }
                }

            } finally {
                br.close();
            }
            int numberOfFiles = Integer.valueOf(args[1]).intValue();
            int numberOfIons = totalIons / numberOfFiles;
            br = new BufferedReader(new FileReader(args[0]));
            File file = new File(args[0]);
            String filePath = file.getParent();
            String fileName = file.getName();
            for (int j = 0; j < numberOfFiles; j++) {
                BufferedWriter bw = new BufferedWriter(new FileWriter(filePath + "\\" + j + "_" + fileName));
                int k = 0;
                String line = br.readLine();
                bw.write(line);
                bw.newLine();
                while (line != null) {
                    line = br.readLine();
                    if (line != null) {
                        bw.write(line);
                        bw.newLine();
                        if (line != null && line.equals("END IONS")) {
                            k += 1;
                        }
                        if (k >= numberOfIons && j != (numberOfFiles - 1)) {
                            break;
                        }
                    }
                }
                bw.close();

            }

            br.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
