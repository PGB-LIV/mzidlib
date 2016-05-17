package bgi.ipeak;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import bgi.ipeak.io.PSMSummaryInformation;
import bgi.ipeak.io.ProteinSummaryInformation;
import bgi.ipeak.io.ReadFastaFile;
import bgi.ipeak.util.Protein;

public class ProteinSummary {

    @Option(name = "-i", required = true, usage = "(required),input peptideSummary file,pg file,percolator result file)")
    private String pep2pro_file;

    private Vector<PSMSummaryInformation> psm_information = new Vector<PSMSummaryInformation>();

    private HashMap<String, Vector<String>> pep_pro_lsit = new HashMap<String, Vector<String>>();
    @Option(name = "-o", required = true, usage = "out put file path")
    private String out_file;
    @Option(name = "-d", required = true, usage = "database fasta type file path")
    private String fasta_file;
    @Option(name = "-t", required = false, usage = "threshold score file path")
    private Double Max_allowFDR;
    @Option(name = "-iT", required = true, usage = "(required), file type:3:iPeak PSMSummary file, 2:pg file,1:percolator result file)")
    private Integer pep_pro_match_infor_type = 0;

    private Boolean rematch = true;
    /*HashMap<String, Protein> proMap: * key:proteins ID,value class Progein* conteins all proteins in the fasta database;	 */
    private HashMap<String, Protein> all_proteinsMap = new HashMap<String, Protein>();//�����б?��fasta�ж�ȡ��	
	/*HashMap<String, Protein> proMap: * key:proteins ID,value class Progein* conteins indentified proteins	 */
    private HashMap<String, Protein> proMap = new HashMap<String, Protein>();
    /*HashMap<String, PeptideMatchInfo> pepMap  	 * keys: ,value: 	 */
    private HashMap<String, PeptideMatchInfo> pepMap = new HashMap<String, PeptideMatchInfo>();

    private Vector<ProteinSummaryInformation> pro_summary = new Vector<ProteinSummaryInformation>();

    public static void main(String args[]) throws IOException {

        ProteinSummary pro_group = new ProteinSummary();
        CmdLineParser parser = new CmdLineParser(pro_group);
        try {
            parser.setUsageWidth(100);
            parser.parseArgument(args);
            System.err.println("\niPeak v1.0(2013-11)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
        } catch (CmdLineException e) {
            System.err.println("\niPeakPercolator v1.0(2013-07)\nWritten by Chaoqin Du (duchaoqin@genomics.cn) in the\n"
                    + "Proteomics Research Group, Department of Science & Technology, BGI-Shenzhen.\n");
            System.err.println("\nUSAGE:\n\tjava -Xmx2g -cp ipeak.jar bgi.ipeak.ProteinSummary [options...]\n");
            parser.printUsage(System.err);
            System.err.println("\nError:");
            System.err.println(e.getMessage());
            System.err.println("");
            return;
        }
        pro_group.get_proteinSummary();
        pro_group.OutPrint();
    }

    public ProteinSummary() {

    }

    public ProteinSummary(String p_pep2pro, String p_outfile, String p_fasta, Integer pp_file_type) {
        this.pep2pro_file = p_pep2pro;
        this.out_file = p_outfile;
        this.fasta_file = p_fasta;
        this.pep_pro_match_infor_type = pp_file_type;
    }

    public ProteinSummary(HashMap<String, Vector<String>> p_pep_pro_list) {
        this.pep_pro_lsit = p_pep_pro_list;
        this.pep_pro_match_infor_type = 4;
    }

    public ProteinSummary(Vector<PSMSummaryInformation> psm_inforInformation, String fasta) {

        this.psm_information = psm_inforInformation;
        this.fasta_file = fasta;
        this.pep_pro_match_infor_type = 0;
    }

    public Vector<ProteinSummaryInformation> get_proteinSummary() throws IOException {
        Read_FastaFile();
        setAnalyse_PepProInfor();
        mergeProteins();
        System.out.println("psm number: " + psm_information.size());
        System.out.println("protein db: " + all_proteinsMap.size());
        System.out.println("pep number: " + pepMap.size());
        System.out.println("protein report: " + proMap.size());
        System.out.println("accepted protein number: " + pro_summary.size());
        return pro_summary;
    }

    private void Read_FastaFile() throws FileNotFoundException {
        if (fasta_file != null) {
            ReadFastaFile read_fasta = new ReadFastaFile(fasta_file);
            all_proteinsMap = read_fasta.get_DbProteins();
        }
    }

    private void setAnalyse_PepProInfor() throws NumberFormatException, IOException {
//		System.out.println("PP_type: "+pep_pro_match_infor_type);
        if (pep_pro_match_infor_type == 0) {
            Read_PSMObject();
        } else if (pep_pro_match_infor_type == 1) {
            Read_PercolatorResult();
        } else if (pep_pro_match_infor_type == 2) {
            Read_pep_proList();
        } else if (pep_pro_match_infor_type == 3) {
            Read_iPeakPSMSummary();
        } else if (pep_pro_match_infor_type == 4) {
            Read_PPMatchHashMap();
        }
    }

    /* read the pep-pro list file; file temp is the percolator out put file;
     * "PSMid	lable	percolator_score	q_value	PEPs	peptide_sequence	proteins" */
    private void Read_PercolatorResult() throws NumberFormatException, IOException {
        BufferedReader psm = new BufferedReader(new FileReader(new File(pep2pro_file)));
        String line = null;
        while ((line = psm.readLine()) != null) {
            String[] items = line.split("\t");
            if (items[0].equals("PSMId")) {
                continue;
            }
            double qvalue = Double.parseDouble(items[2]);
            if (qvalue > Max_allowFDR) {
                continue;
            }
            String pepSeq = items[4];
            if (!pepMap.containsKey(pepSeq)) {
                pepMap.put(pepSeq, new PeptideMatchInfo());
            }
            pepMap.get(pepSeq).setSeq(pepSeq);
            String specNum = items[0];
            pepMap.get(pepSeq).getSpectra().add(specNum);

            if (rematch) {
                Vector<String> protein_accs = get_ProteinsOfthePeptide_(pepSeq);
                for (String acc : protein_accs) {

                    if (!pepMap.get(pepSeq).getPros().contains(acc)) {
                        pepMap.get(pepSeq).getPros().add(acc);
                    }

                    if (proMap.containsKey(acc)) {
                        if (!proMap.get(acc).getPeps().contains(pepSeq)) {
                            proMap.get(acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = all_proteinsMap.get(acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(acc, pro);
                    }
                }
            } else {
                for (int i = 5; i < items.length; i++) {
                    if (!pepMap.get(pepSeq).getPros().contains(items[i])) {
                        pepMap.get(pepSeq).getPros().add(items[i]);
                    }
                    if (proMap.containsKey(items[i])) {
                        if (!proMap.get(items[i]).getPeps().contains(pepSeq)) {
                            proMap.get(items[i]).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = new Protein();
                        pro.setAccession(items[i]);
                        pro.getPeps().add(pepSeq);
                        proMap.put(items[i], pro);
                    }
                }
            }
        }
        psm.close();
    }

    private void Read_pep_proList() throws NumberFormatException, IOException {
        System.out.println("file type\n\tone line one peptide-protein pair,eg:\npeptide_sequence\tprotein1;protein2;...;proteinN\n");
        BufferedReader psm = new BufferedReader(new FileReader(new File(pep2pro_file)));
        String line = null;
        Integer pep_id = 0;
        while ((line = psm.readLine()) != null) {
            pep_id++;
            String[] items = line.split("\t");

            String pepSeq = items[0];
            if (!pepMap.containsKey(pepSeq)) {
                pepMap.put(pepSeq, new PeptideMatchInfo());
            }
            pepMap.get(pepSeq).setSeq(pepSeq);

            String specNum = String.valueOf(pep_id);
            pepMap.get(pepSeq).getSpectra().add(specNum);

            if (rematch) {
                Vector<String> protein_accs = get_ProteinsOfthePeptide_(pepSeq);
                for (String acc : protein_accs) {

                    if (!pepMap.get(pepSeq).getPros().contains(acc)) {
                        pepMap.get(pepSeq).getPros().add(acc);
                    }

                    if (proMap.containsKey(acc)) {
                        if (!proMap.get(acc).getPeps().contains(pepSeq)) {
                            proMap.get(acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = all_proteinsMap.get(acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(acc, pro);
                    }
                }
            } else {
                String[] proteins = items[1].split(";");
                for (int i = 0; i < proteins.length; i++) {
                    if (!pepMap.get(pepSeq).getPros().contains(proteins[i])) {
                        pepMap.get(pepSeq).getPros().add(proteins[i]);
                    }
                    if (proMap.containsKey(proteins[i])) {
                        if (!proMap.get(proteins[i]).getPeps().contains(pepSeq)) {
                            proMap.get(proteins[i]).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = new Protein();
                        pro.setAccession(proteins[i]);
                        pro.getPeps().add(pepSeq);
                        proMap.put(proteins[i], pro);
                    }
                }
            }
        }
        psm.close();
    }

    /* read the pep-pro list file; file temp is the percolator out put file;
     * "PSMid	lable	percolator_score	q_value	PEPs	peptide_sequence	proteins" */
    private void Read_iPeakPSMSummary() throws NumberFormatException, IOException {
        BufferedReader psm = new BufferedReader(new FileReader(new File(pep2pro_file)));
        String line = psm.readLine();
        while ((line = psm.readLine()) != null) {
            String[] items = line.split("\t");
            double threshold_score = Double.parseDouble(items[3]);
            if (threshold_score > Max_allowFDR) {
                continue;
            }
            String pepSeq = items[2];
            if (!pepMap.containsKey(pepSeq)) {
                pepMap.put(pepSeq, new PeptideMatchInfo());
            }
            pepMap.get(pepSeq).setSeq(pepSeq);
            String specNum = items[1];
            pepMap.get(pepSeq).getSpectra().add(specNum);

            if (rematch) {
                Vector<String> protein_accs = get_ProteinsOfthePeptide_(pepSeq);
                for (String acc : protein_accs) {

                    if (!pepMap.get(pepSeq).getPros().contains(acc)) {
                        pepMap.get(pepSeq).getPros().add(acc);
                    }

                    if (proMap.containsKey(acc)) {
                        if (!proMap.get(acc).getPeps().contains(pepSeq)) {
                            proMap.get(acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = all_proteinsMap.get(acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(acc, pro);
                    }
                }
            } else {
                for (int i = 12; i < items.length; i++) {
                    if (!pepMap.get(pepSeq).getPros().contains(items[i])) {
                        pepMap.get(pepSeq).getPros().add(items[i]);
                    }
                    if (proMap.containsKey(items[i])) {
                        if (!proMap.get(items[i]).getPeps().contains(pepSeq)) {
                            proMap.get(items[i]).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = new Protein();
                        pro.setAccession(items[i]);
                        pro.getPeps().add(pepSeq);
                        proMap.put(items[i], pro);
                    }
                }
            }
        }
        psm.close();
    }

    private void Read_PSMObject() {
        System.out.println("read from psm project");
        Iterator<PSMSummaryInformation> iterator = psm_information.iterator();
        pep_pro_lsit.clear();
        while (iterator.hasNext()) {
            PSMSummaryInformation the_psm = iterator.next();
//			System.out.println(the_psm.getIs_decoy()+"\t"+the_psm.getPass_threshold());
            if (the_psm.getIs_decoy()) {
                continue;
            }
            if (!the_psm.getPass_threshold()) {
                continue;
            }
//			System.out.println(the_psm.getPeptide_sequence()+"\t"+the_psm.getProteins());
            String pepSeq = the_psm.getPeptide_sequence();
            if (!pepMap.containsKey(pepSeq)) {
                pepMap.put(pepSeq, new PeptideMatchInfo());
            }
            pepMap.get(pepSeq).setSeq(pepSeq);
            String specNum = the_psm.getSpectrum_index();
            pepMap.get(pepSeq).getSpectra().add(specNum);

            if (rematch) {
                Vector<String> protein_accs = new Vector<String>();
                if (pep_pro_lsit.containsKey(pepSeq)) {
                    protein_accs = pep_pro_lsit.get(pepSeq);
                } else {
                    protein_accs = get_ProteinsOfthePeptide_(pepSeq);
                    pep_pro_lsit.put(pepSeq, protein_accs);
                }
                the_psm.setProteins(protein_accs);
                for (String acc : protein_accs) {

                    if (!pepMap.get(pepSeq).getPros().contains(acc)) {
                        pepMap.get(pepSeq).getPros().add(acc);
                    }

                    if (proMap.containsKey(acc)) {
                        if (!proMap.get(acc).getPeps().contains(pepSeq)) {
                            proMap.get(acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = all_proteinsMap.get(acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(acc, pro);
                    }
                }
            } else {
                for (String protein_acc : the_psm.getProteins()) {
                    if (!pepMap.get(pepSeq).getPros().contains(protein_acc)) {
                        pepMap.get(pepSeq).getPros().add(protein_acc);
                    }
                    if (proMap.containsKey(protein_acc)) {
                        if (!proMap.get(protein_acc).getPeps().contains(pepSeq)) {
                            proMap.get(protein_acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = new Protein();
                        pro.setAccession(protein_acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(protein_acc, pro);
                    }
                }
            }
        }
    }

    private void Read_PPMatchHashMap() throws NumberFormatException, IOException {
        Iterator<String> iterator = pep_pro_lsit.keySet().iterator();
        Integer pep_id = 0;
        while (iterator.hasNext()) {
            pep_id++;
            String pepSeq = iterator.next();
            if (!pepMap.containsKey(pepSeq)) {
                pepMap.put(pepSeq, new PeptideMatchInfo());
            }
            pepMap.get(pepSeq).setSeq(pepSeq);

            String specNum = String.valueOf(pep_id);
            pepMap.get(pepSeq).getSpectra().add(specNum);

            if (rematch) {
                Vector<String> protein_accs = get_ProteinsOfthePeptide_(pepSeq);
                for (String acc : protein_accs) {

                    if (!pepMap.get(pepSeq).getPros().contains(acc)) {
                        pepMap.get(pepSeq).getPros().add(acc);
                    }

                    if (proMap.containsKey(acc)) {
                        if (!proMap.get(acc).getPeps().contains(pepSeq)) {
                            proMap.get(acc).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = all_proteinsMap.get(acc);
                        pro.getPeps().add(pepSeq);
                        proMap.put(acc, pro);
                    }
                }
            } else {
                Vector<String> proteins = pep_pro_lsit.get(pepSeq);
                for (int i = 0; i < proteins.size(); i++) {
                    if (!pepMap.get(pepSeq).getPros().contains(proteins.get(i))) {
                        pepMap.get(pepSeq).getPros().add(proteins.get(i));
                    }
                    if (proMap.containsKey(proteins.get(i))) {
                        if (!proMap.get(proteins.get(i)).getPeps().contains(pepSeq)) {
                            proMap.get(proteins.get(i)).getPeps().add(pepSeq);
                        }
                    } else {
                        Protein pro = new Protein();
                        pro.setAccession(proteins.get(i));
                        pro.getPeps().add(pepSeq);
                        proMap.put(proteins.get(i), pro);
                    }
                }
            }
        }
    }

    private Vector<String> get_ProteinsOfthePeptide_(String peptideSequence) {
        if (all_proteinsMap.isEmpty()) {
            System.err.println("protein group error: the fasta not read in. you can not rematch the peptide to proteins!\n");
            System.exit(3);
        }
        Vector<String> protein_acc = new Vector<String>();
        for (Protein protein : all_proteinsMap.values()) {
            if (protein.getProtein_sequence().contains(peptideSequence)) {
                protein_acc.add(protein.getAccession());
            }
        }
        if (protein_acc.isEmpty()) {
            System.err.println("peptide: " + peptideSequence + " can not match any protein in " + fasta_file);
        }
        return protein_acc;
    }

    /*	 * merge the proteins;	 */
    public void mergeProteins() throws IOException {

        Iterator<String> iterPro = proMap.keySet().iterator();
        int clusterID = 0;
        ArrayList<String> clusterPro = new ArrayList<String>();
        while (iterPro.hasNext()) {
            String proAcc = iterPro.next();
            Protein proTemp = proMap.get(proAcc);

            if (proTemp.getClusterID() == -1) {
                clusterID++;
                proTemp.setClusterID(clusterID);
                clusterPro.add(proAcc);

                ArrayList<String> peps = new ArrayList<String>();
                for (int i = 0; i < proTemp.getPeps().size(); i++) {
                    peps.add(proTemp.getPeps().get(i));
                }
                for (int i = 0; i < peps.size(); i++) {
                    PeptideMatchInfo pepTemp = pepMap.get(peps.get(i));
                    ArrayList<String> proList = pepTemp.getPros();
                    for (int j = 0; j < proList.size(); j++) {
                        if (proMap.get(proList.get(j)).getClusterID() == -1) {
                            proMap.get(proList.get(j)).setClusterID(clusterID);
                            proTemp.getSameset().add(proList.get(j));
                            for (int k = 0; k < proMap.get(proList.get(j)).getPeps().size(); k++) {
                                if (!peps.contains(proMap.get(proList.get(j)).getPeps().get(k))) {
                                    peps.add(proMap.get(proList.get(j)).getPeps().get(k));
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < clusterPro.size(); i++) {
            ArrayList<String> pepArray = new ArrayList<String>();
            Protein proTemp = proMap.get(clusterPro.get(i));
            ArrayList<Protein> proList = new ArrayList<Protein>();
            proList.add(proTemp);

            for (int j = 0; j < proTemp.getSameset().size(); j++) {
                boolean addFlag = false;

                for (int k = proList.size() - 1; k >= 0; k--) {
                    if (proMap.get(proTemp.getSameset().get(j)).getPeps().size()
                            > proList.get(k).getPeps().size()) {
                        proList.add(k + 1, proMap.get(proTemp.getSameset().get(j)));
                        addFlag = true;
                        break;
                    } else if (proMap.get(proTemp.getSameset().get(j)).getPeps().size()
                            == proList.get(k).getPeps().size()) {
                        if (proMap.get(proTemp.getSameset().get(j)).getLength() > proList.get(k).getLength()) {
                            proList.add(k + 1, proMap.get(proTemp.getSameset().get(j)));
                            addFlag = true;
                            break;
                        }
                    }
                }
                if (addFlag == false) {
                    proList.add(0, proMap.get(proTemp.getSameset().get(j)));
                }
            }

            proTemp.setSameset(new ArrayList<String>());
            for (int j = 0; j < proList.size(); j++) {
                for (int k = j + 1; k < proList.size(); k++) {
                    if (proList.get(j).getPeps().containsAll(proList.get(k).getPeps())
                            && proList.get(k).getPeps().containsAll(proList.get(j).getPeps())) {
                        proList.get(j).getSameset().add(proList.get(k).getAccession());
                        proList.remove(k);
                        k--;
                    }
                }
            }

            for (int j = 0; j < proList.size(); j++) {
                for (int k = 0; k < proList.get(j).getPeps().size(); k++) {
                    if (!pepArray.contains(proList.get(j).getPeps().get(k))) {
                        pepArray.add(proList.get(j).getPeps().get(k));
                    }
                }
            }

            Integer pro_group_id = 0;
            for (int j = 0; j < proList.size(); j++) {

                Boolean accepted_protein = false;
                for (int k = 0; k < proList.get(j).getPeps().size(); k++) {
                    if (pepMap.get(proList.get(j).getPeps().get(k)).getPros().size() == 1) {
                        accepted_protein = true;
                    }
                }

                if (!accepted_protein) {
                    ArrayList<String> pepComp = new ArrayList<String>();
                    for (int k = 0; k < proList.size(); k++) {
                        if (k != j) {
                            for (int l = 0; l < proList.get(k).getPeps().size(); l++) {
                                if (!pepComp.contains(proList.get(k).getPeps().get(l))) {
                                    pepComp.add(proList.get(k).getPeps().get(l));
                                }
                            }
                        }
                    }
                    if (pepArray.containsAll(pepComp) && pepComp.containsAll(pepArray)) {
                        proList.remove(j);
                        j--;
                    } else {
                        accepted_protein = true;
                    }
                }
                if (accepted_protein) {
                    proList.get(j).setNonsense(false);
                }
            }

            for (int j = 0; j < proList.size(); j++) {
                ArrayList<String> pepComp = new ArrayList<String>();
                for (int k = 0; k < proList.size(); k++) {
                    if (k != j) {
                        for (int l = 0; l < proList.get(k).getPeps().size(); l++) {
                            if (!pepComp.contains(proList.get(k).getPeps().get(l))) {
                                pepComp.add(proList.get(k).getPeps().get(l));
                            }
                        }
                    }
                }
                for (int k = 0; k < pepArray.size(); k++) {
                    if (!pepComp.contains(pepArray.get(k))) {
                        pepMap.get(pepArray.get(k)).setUniq(true);
                    }
                }
            }

            for (int j = 0; j < proList.size(); j++) {
                pro_group_id++;
                proList.get(j).setGroupID("PG_" + proList.get(j).getClusterID() + "_" + pro_group_id);
                ArrayList<String> same_setproteis = proList.get(j).getSameset();
                String behalf_protein = proList.get(j).getAccession();
                Integer protein_length = proList.get(j).getLength();

                if (!same_setproteis.isEmpty()) {
                    for (String pro_acc : same_setproteis) {
                        if (proMap.get(pro_acc).getLength() < protein_length) {
                            behalf_protein = proMap.get(pro_acc).getAccession();
                        }
                    }
                }

                ProteinSummaryInformation the_proSum = new ProteinSummaryInformation();
                String protein_group_id = proList.get(j).getGroupID();
                Vector<String> same_set = new Vector<String>();
                for (String sameset_pro : proList.get(j).getSameset()) {
                    if (behalf_protein.equals(sameset_pro)) {
                        same_set.add(proList.get(j).getAccession());
                    } else {
                        same_set.add(sameset_pro);
                    }
                }
                ArrayList<String> peptides = proMap.get(behalf_protein).getPeps();
                Vector<String> uniq_peps = new Vector<String>();
                Vector<String> razor_pep = new Vector<String>();
                Integer uniq_psm_number = 0;
                Integer razor_psm_number = 0;
                for (String pep : peptides) {
                    //						System.out.println(pepMap.get(pep).isUniq()+"\t"+pepMap.get(pep).getSpectra().size());
                    if (pepMap.get(pep).isUniq()) {
                        uniq_peps.add(pep);
                        uniq_psm_number += pepMap.get(pep).getSpectra().size();
                    } else {
                        razor_pep.add(pep);
                        razor_psm_number += pepMap.get(pep).getSpectra().size();
                    }
                }

                Integer psm_number = uniq_psm_number + razor_psm_number;
                //					System.out.println(uniq_peps.size()+"\t"+razor_pep.size()+"\t"+psm_number+"\t"+uniq_psm_number+"\t"+razor_psm_number);
                String reference = proMap.get(behalf_protein).getDescription();

                the_proSum.setProtein_acc(behalf_protein);
                the_proSum.setMass(proMap.get(behalf_protein).getMass());
                the_proSum.setPeptide_number(peptides.size());
                the_proSum.setProtein_group_id(protein_group_id);
                the_proSum.setPsm_number(psm_number);
                the_proSum.setRazor_peptides(razor_pep);
                the_proSum.setReference(reference);
                the_proSum.setSameset_proteins(same_set);
                the_proSum.setUniq_pep_number(uniq_peps.size());
                the_proSum.setUniq_peptides(uniq_peps);
                the_proSum.setUniq_psm_number(uniq_psm_number);
                pro_summary.add(the_proSum);
            }
        }

    }

    public void OutPrint() throws IOException {
        BufferedWriter output = new BufferedWriter(new FileWriter(new File(out_file)));
        output.write("Accession\tMass\tPeptideSeqs\tPepIsUniqe\tNumOfUniqPeps\tNumOfUniqSpectra\tSpectra\tSameSet\tDescription");
        output.newLine();
        output.flush();
        for (ProteinSummaryInformation proteinSummaryInformation : pro_summary) {
            proteinSummaryInformation.print_protein(output);
        }
        output.close();
    }

    static class PeptideMatchInfo {

        /**
         * * PeptideInfor sequence.
         */
        private String seq = "";
        /**
         * * The spectra matched the peptide.
         */
        private ArrayList<String> spectra = new ArrayList<String>();
        /**
         * * The proteins contained the peptide.
         */
        private ArrayList<String> pros = new ArrayList<String>();
        /**
         * * Is this peptide contained by only one protein group?
         */
        private boolean uniq = false;

        public PeptideMatchInfo() {

        }

        public String getSeq() {
            return seq;
        }

        public void setSeq(String seq) {
            this.seq = seq;
        }

        public ArrayList<String> getSpectra() {
            return spectra;
        }

        public void setSpectra(ArrayList<String> spectra) {
            this.spectra = spectra;
        }

        public ArrayList<String> getPros() {
            return pros;
        }

        public void setPros(ArrayList<String> pros) {
            this.pros = pros;
        }

        public boolean isUniq() {
            return uniq;
        }

        public void setUniq(boolean uniq) {
            this.uniq = uniq;
        }
    }
}
