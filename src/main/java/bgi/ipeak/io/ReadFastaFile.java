package bgi.ipeak.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import bgi.ipeak.util.Protein;

public class ReadFastaFile {

    private String fasta_file;
    private String tag = ">";
    HashMap<String, Protein> dbProMap;

    public ReadFastaFile(String fasta_file, String tag) {

        this.fasta_file = fasta_file;
        this.tag = tag;
    }

    public ReadFastaFile(String fasta_file) {

        this.fasta_file = fasta_file;
    }

    public ReadFastaFile() {

    }

    public HashMap<String, Protein> get_DbProteins() throws FileNotFoundException {
        Read_FastaFile();
        return dbProMap;
    }

    private void Read_FastaFile() throws FileNotFoundException {
        dbProMap = new HashMap<String, Protein>();
        if (fasta_file != null) {
            Scanner dbPro = new Scanner(new File(fasta_file));
            dbPro.useDelimiter(tag);
            HashMap<String, Double> acidsMass = new HashMap<String, Double>();
            acidsMass.put("A", 71.037114);
            acidsMass.put("B", 114.534940);
            acidsMass.put("C", 160.030649);
            acidsMass.put("D", 115.026943);
            acidsMass.put("E", 129.042593);
            acidsMass.put("F", 147.068414);
            acidsMass.put("G", 57.021464);
            acidsMass.put("H", 137.058912);
            acidsMass.put("I", 113.084064);
            acidsMass.put("J", 0.000000);
            acidsMass.put("K", 128.094963);
            acidsMass.put("L", 113.084064);
            acidsMass.put("M", 131.040485);
            acidsMass.put("N", 114.042927);
            acidsMass.put("O", 0.000000);
            acidsMass.put("P", 97.052764);
            acidsMass.put("Q", 128.058578);
            acidsMass.put("R", 156.101111);
            acidsMass.put("S", 87.032028);
            acidsMass.put("T", 101.047679);
            acidsMass.put("U", 150.953630);
            acidsMass.put("V", 99.068414);
            acidsMass.put("W", 186.079313);
            acidsMass.put("X", 111.000000);
            acidsMass.put("Y", 163.063329);
            acidsMass.put("Z", 128.550590);
            Integer proID = 0;
            String proDesc = "";

            while (dbPro.hasNext()) {
                String proLine = dbPro.next();
                String[] proInfo = proLine.split("\n");
                if (proInfo.length == 1) {
                    proDesc += proInfo[0];
                    continue;
                } else {
                    proDesc += proInfo[0];
                }
                String[] proAccAndDesc = proDesc.split("\\s");
                proDesc = "";
                proID++;
                Double proteinMass = 0.0;
                Integer proteinLength = 0;
                String proAcc = proAccAndDesc[0];

                String proteinReference = "";
                for (int j = 1; j < proAccAndDesc.length; j++) {
                    proteinReference += proAccAndDesc[j] + " ";
                }
                String proSeq = "";
                for (int i = 1; i < proInfo.length; i++) {
                    proSeq += proInfo[i];
                    String[] acidsSeq = proInfo[i].split("");
                    for (String acids : acidsSeq) {
                        if (acidsMass.containsKey(acids)) {
                            proteinMass += acidsMass.get(acids);
                            proteinLength++;
                        }
                    }
                }
                if (proteinReference.replace(" ", "").isEmpty()) {
                    proteinReference = "-";
                }
                if (proteinReference.length() > 50) {
                    proteinReference = proteinReference.substring(0, 50) + "...";
                }
                Protein proTemp = new Protein();
                proTemp.setAccession(proAcc);
                proTemp.setDescription(proteinReference);
                proTemp.setLength(proteinLength);
                proTemp.setMass(proteinMass);
                proTemp.setProtein_sequence(proSeq);
                dbProMap.put(proAcc, proTemp);
            }
            dbPro.close();
        }
    }

}
