package bgi.ipeak.percolator.paser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import bgi.ipeak.io.ReadFastaFile;
import bgi.ipeak.util.Protein;
import de.proteinms.omxparser.OmssaOmxFile;
import de.proteinms.omxparser.util.MSHitSet;
import de.proteinms.omxparser.util.MSHits;
import de.proteinms.omxparser.util.MSMZHit;
import de.proteinms.omxparser.util.MSModHit;
import de.proteinms.omxparser.util.MSPepHit;
import de.proteinms.omxparser.util.MSRequest;
import de.proteinms.omxparser.util.MSResponse;
import de.proteinms.omxparser.util.MSSearchSettings;
import de.proteinms.omxparser.util.MSSpectrum;
import de.proteinms.omxparser.util.OmssaModification;

/**
 * This class parse OMSSA search result using omssa-parser.
 *
 * @author Guilin Li
 */
public class OmssaParser {

    private String omssa_omx_file;
    private String output_feature_file;
    private String mods_file;
    private String usemod_file;
    private String fasta_file;
    private Integer rank = 1;
    private Double dm = 0.0;
    private Boolean dmppm = true;
    private String decoyRegex;
    private Boolean us_ms2 = true;
    private Vector<SpectraFeature> omssa_feature_list = new Vector<SpectraFeature>();
    private String feature_title = "";
    ArrayList<Integer> ionSearch = new ArrayList<Integer>();

    public OmssaParser(String omssa_omx_file, String output_feature_file,
            String usemod_file, String fasta, String decoyRegex) {
        this.omssa_omx_file = omssa_omx_file;
        this.output_feature_file = output_feature_file;
        this.mods_file = bgi.ipeak.util.Properties.getModFile_path();
        this.usemod_file = usemod_file;
        this.fasta_file = fasta;
        this.decoyRegex = decoyRegex;
        // TODO �Զ���ɵĹ��캯����
    }

    public void get_feature_file() throws Exception {
        feature_title = "Query\tTargetDecoy\tLog10Evalue\tMass\tCharge\tDeltaMass\t"
                + "DeltaMassPPM\tabsDM\tabsDMppm\tisoDM\tisoDMppm\tVarModRatio\t"
                + "TotalIntensity\tMatchedIonInt\trelTotMatchedIonInt\tMaxMatchedIonInt\t";
        if (us_ms2) {
            feature_title += "FragError\tFragDeltaM_Med\tFragDeltaM_MedPPM\tFragDeltaM_Iqr\tFragDeltaM_IqrPPM\t";
        }
        feature_title += "Longest\tQMatch\tEnzTryN\tEnzTryC\tEnzTryN\tPepLen\tLog10Pvalue";

        parseOmssaXML();
        print_features();
        clear();
    }

    public void parseOmssaXML() throws Exception {
        System.out.println("Parse OMSSA omx search result " + omssa_omx_file + "...");
        ReadFastaFile readFasta = new ReadFastaFile(fasta_file);
        HashMap<String, Protein> proMap = readFasta.get_DbProteins();
        HashMap<String, String> realignedPeps = new HashMap<String, String>();

        File omssa_file = new File(omssa_omx_file);
        String name = omssa_file.getName();
        name = name.substring(0, name.lastIndexOf("."));
        OmssaOmxFile omxFile = new OmssaOmxFile(omssa_omx_file, mods_file, usemod_file);
        List<MSRequest> mr = omxFile.getParserResult().MSSearch_request.MSRequest;
        List<MSResponse> mp = omxFile.getParserResult().MSSearch_response.MSResponse;

        for (MSRequest mrTemp : mr) {
            for (Integer iType : mrTemp.MSRequest_settings.MSSearchSettings.MSSearchSettings_ionstosearch.MSIonType) {
                if (!ionSearch.contains(iType)) {
                    ionSearch.add(iType);
                }
            }
        }

        if (ionSearch.size() == 0) {
            ionSearch.add(1);
            ionSearch.add(4);
            for (int i = 0; i < 4; i++) {
                feature_title += "\tfracIonSeries_" + i + "\trelIonSeriesInt_" + i;
            }
        } else {
            for (int i = 0; i < ionSearch.size() * 2; i++) {
                feature_title += "\tfracIonSeries_" + i + "\trelIonSeriesInt_" + i;
            }
        }
        feature_title += "\tPeptide\tProteins";
        //MSResponse 
        for (int i = 0; i < mr.size(); i++) {
            MSRequest msRequest = mr.get(i);
            MSSearchSettings msSearch = msRequest.MSRequest_settings.MSSearchSettings;
            Map<Integer, MSSpectrum> specMap = msRequest.MSRequest_spectra.MSSpectrumset.MSSpectrum;
            int requestMzScale = msSearch.MSSearchSettings_scale;
            double msmsTol = msSearch.MSSearchSettings_msmstol;
            //variable modifications, related to search parameter "-mv".
            List<Integer> varModList = msSearch.MSSearchSettings_variable.MSMod;
            ArrayList<Modification> varMod = new ArrayList<Modification>();
            for (Integer m : varModList) {
                OmssaModification amod = omxFile.getModifications().get(m);
                Modification myMod = new Modification();
                myMod.setFixed(false);
                myMod.setModType(m);
                //System.out.println(amod.getModName()+"\t"+amod.getModNumber()+"\t"+amod.getModMonoMass());
                myMod.setMonoMass(amod.getModMonoMass());
                myMod.setName(amod.getModName());
                if (amod.getModResidues().size() == 0) {
                    myMod.setResidues(null);
                } else {
                    myMod.setResidues(amod.getModResidues());
                    myMod.setCompactResidues(amod.getModResiduesAsCompactString().toUpperCase());
                }
                if (amod.getModName().toLowerCase().indexOf("peptide c-term") != -1 || amod.getModName().toLowerCase().indexOf("peptide cterm") != -1) {
                    myMod.setPepcterm(true);
                } else if (amod.getModName().toLowerCase().indexOf("peptide n-term") != -1 || amod.getModName().toLowerCase().indexOf("peptide nterm") != -1) {
                    myMod.setPepnterm(true);
                } else if (amod.getModName().toLowerCase().indexOf("protein c-term") != -1 || amod.getModName().toLowerCase().indexOf("protein cterm") != -1) {
                    myMod.setProcterm(true);
                } else if (amod.getModName().toLowerCase().indexOf("protein n-term") != -1 || amod.getModName().toLowerCase().indexOf("protein nterm") != -1) {
                    myMod.setPronterm(true);
                } else if (amod.getModName().toLowerCase().indexOf("nterm") != -1) {
                    myMod.setPronterm(true);
                    myMod.setPepnterm(true);
                } else if (amod.getModName().toLowerCase().indexOf("cterm") != -1) {
                    myMod.setProcterm(true);
                    myMod.setPepcterm(true);
                }
                varMod.add(myMod);
            }
            Iterator<Integer> iterRequest = specMap.keySet().iterator();

            MSResponse msResponse = mp.get(i);
            int responseMzScale = msResponse.MSResponse_scale;
            Map<Integer, MSHitSet> msHitMap = msResponse.MSResponse_hitsets.MSHitSet;

            while (iterRequest.hasNext()) {
                Integer specNum = iterRequest.next();
                MSHitSet msh = msHitMap.get(specNum);
                List<MSHits> mhl = msh.MSHitSet_hits.MSHits;
                if (mhl.size() == 0) {
                    continue;
                }
                //spectrum information.
                MSSpectrum specTemp = specMap.get(specNum);
                SpectraFeature spf = new SpectraFeature();
                List<Integer> mzList = specTemp.MSSpectrum_mz.MSSpectrum_mz_E;
                List<Integer> abundList = specTemp.MSSpectrum_abundance.MSSpectrum_abundance_E;
                HashMap<Integer, Integer> mzToIntensity = new HashMap<Integer, Integer>();
                Double abundScale = specTemp.MSSpectrum_iscale;
                double[] intensity = new double[abundList.size()];
                double[] mz = new double[abundList.size()];
                for (int j = 0; j < abundList.size(); j++) {
                    mzToIntensity.put(mzList.get(j), abundList.get(j));
                    intensity[j] = abundList.get(j) / abundScale;
                    mz[j] = mzList.get(j) / (double) requestMzScale;
                }
                spf.setMz(mz);
                spf.setIntensity(intensity);
                spf.setTotIne(spf.calLogTotInt());
                //only the front <rank> match be concerned. 
                for (int x = 0; x < (rank > mhl.size() ? mhl.size() : rank); x++) {
                    MSHits mh = mhl.get(x);
                    spf.setCharge(Byte.parseByte(String.valueOf(mh.MSHits_charge)));
                    //include the situation "P-value=0" or "E-value=0".
                    spf.setPepEvalue(mh.MSHits_evalue == 0 ? 100 : -Math.log(mh.MSHits_evalue) / Math.log(10));
                    spf.setPvalue(mh.MSHits_pvalue == 0 ? 100 : Math.log(mh.MSHits_pvalue) / Math.log(10));
                    //set a threshold for log-transformed P-value and E-value.
                    if (spf.getPepEvalue() > 100) {
                        spf.setPepEvalue(100);
                    }
                    if (spf.getPvalue() < -100) {
                        spf.setPvalue(-100);
                    }
                    spf.setPepString(mh.MSHits_pepstring);
                    spf.setPepLen(mh.MSHits_pepstring.length());
                    spf.setdM((mh.MSHits_mass - mh.MSHits_theomass) / (double) requestMzScale);
                    spf.setMrCalc(mh.MSHits_theomass / (double) requestMzScale);//experimental neutral mass of peptide in Da.				//calculated peptide mass.
                    //spf.setMrCalc(specTemp.MSSpectrum_charge.MSSpectrum_charge_E.get(0)*specTemp.MSSpectrum_precursormz/(double)requestMzScale
                    //-spf.getdM()-specTemp.MSSpectrum_charge.MSSpectrum_charge_E.get(0)*1.007276);
                    spf.setdMppm(spf.getdM() / spf.getMrCalc() * 1000000);
                    //spf.setdMppm(spf.getdM()*requestMzScale*1000000/mh.MSHits_mass);
                    spf.setAbsDM(Math.abs(spf.getdM()));
                    spf.setAbsDMppm(Math.abs(spf.getdMppm()));
                    //filter the PSM by dm.
                    if (dm != 0) {
                        if (dmppm && spf.getAbsDMppm() > dm) {
                            continue;
                        }
                        if (!dmppm && spf.getAbsDM() > dm) {
                            continue;
                        }
                    }
                    /*if(mhl.size()>1){
                     MSHits mhTemp=mhl.get(x+1);
                     double secondScore=mhTemp.MSHits_evalue==0?100:-Math.log(mhTemp.MSHits_evalue)/Math.log(10);
                     spf.setDeltaScore(spf.getPepEvalue()-secondScore);
                     }else{
                     spf.setDeltaScore(0);
                     }*/
                    //isotope mass error
                    double isoDm_0 = Math.abs(spf.getAbsDM());
                    double isoDm_1 = Math.abs(spf.getAbsDM() - 1.007825);
                    double isoDm_2 = Math.abs(spf.getAbsDM() - 2 * 1.007825);
                    double isoDm = isoDm_0;
                    if (isoDm_1 < isoDm) {
                        isoDm = isoDm_1;
                    }
                    if (isoDm_2 < isoDm) {
                        isoDm = isoDm_2;
                    }
                    spf.setIsoDM(isoDm);
                    spf.setIsoDMppm(isoDm / spf.getMrCalc() * 1000000);
                    //modifications, only the variable concerned.
                    String pepSeq = mh.MSHits_pepstring;
                    String preAA = mh.MSHits_pepstart, postAA = mh.MSHits_pepstop;
                    pepSeq = pepSeq.toUpperCase();
                    List<MSModHit> mml = mh.MSHits_mods.MSModHit;
                    int modifiable = countVarMod(varMod, pepSeq, preAA, postAA);
                    int modified = 0;
                    ArrayList<Integer> varModCount = new ArrayList<Integer>();
                    String modSites = "";
                    for (int j = 0; j < mml.size(); j++) {
                        MSModHit modHit = mml.get(j);
                        modSites += modHit.MSModHit_modtype.MSMod + "@" + modHit.MSModHit_site + ";";
                        if (varModList.contains(modHit.MSModHit_modtype.MSMod)) {
                            modified++;
                            if (!varModCount.contains(modHit.MSModHit_modtype.MSMod)) {
                                varModCount.add(modHit.MSModHit_modtype.MSMod);
                            }
                        }
                    }
                    spf.setVarmodsCount(varModCount.size());
                    spf.setVarmodsRatio(modifiable == 0 ? 0 : modified / (double) modifiable);
                    //spf.setVarmodsRatio(modifiable==0?0:modified/modifiable);//this is the original code.
                    spf.setModifiable(modifiable);
                    spf.setModified(modified);
                    spf.setMod(modSites);
                    //protein information.
                    String proteins = "";
                    int targetDecoy = 1;
                    if (decoyRegex.equals("target")) {
                        if (!realignedPeps.containsKey(pepSeq)) {
                            proteins = realignPeptideInDB(pepSeq, proMap);
                            if (proteins.isEmpty()) {
                                System.err.println("Can't find this sequence " + pepSeq + " in database!");
                                continue;
                            }
                            realignedPeps.put(pepSeq, proteins);
                        } else {
                            proteins = realignedPeps.get(pepSeq);
                        }
                    } else if (decoyRegex.equals("decoy")) {
                        targetDecoy = -1;
                        List<MSPepHit> pepProHit = mh.MSHits_pephits.MSPepHit;
                        for (int j = 0; j < pepProHit.size(); j++) {
                            MSPepHit pepProHitTemp = pepProHit.get(j);
                            String proteinID = pepProHitTemp.MSPepHit_accession;
                            proteins += proteinID + "\t";
                        }
                    } else {
                        if (!realignedPeps.containsKey(pepSeq)) {
                            proteins = realignPeptideInDB(pepSeq, proMap);
                            if (proteins.isEmpty()) {
                                System.err.println("Can't find this sequence " + pepSeq + " in database!");
                                continue;
                            }
                            realignedPeps.put(pepSeq, proteins);
                        } else {
                            proteins = realignedPeps.get(pepSeq);
                        }
                        Pattern proPattern = Pattern.compile(decoyRegex);
                        String[] accs = proteins.split("\t");
                        for (int k = 0; k < accs.length; k++) {
                            Matcher match = proPattern.matcher(proMap.get(accs[k]).getAccession());
                            if (match.find()) {
                                targetDecoy = -1;
                                break;
                            }
                        }
                    }

                    spf.setProteins(proteins);
                    spf.setTargetDecoy(targetDecoy);

                    //The enzymatic features.
                    String[] pepAA = pepSeq.split("");
                    //System.out.println(preAA+"\t"+pepSeq+"\t"+postAA);
                    boolean ctermEnz = isEnzymatic(pepAA[pepAA.length - 1], postAA);
                    spf.setcTermTryptic(ctermEnz);
                    boolean ntermEnz = isEnzymatic(preAA, pepAA[0]);
                    spf.setnTermTryptic(ntermEnz);
                    int enzN = 0;
                    for (int k = 1; k < pepAA.length; k++) {
                        if (isEnzymatic(pepAA[k - 1], pepAA[k])) {
                            enzN++;
                        }
                    }
                    spf.setEnzN(enzN);

                    //fragment ion series.
                    List<MSMZHit> mzHit = mh.MSHits_mzhits.MSMZHit;
                    HashMap<String, Double> ionSeriesInt = new HashMap<String, Double>();
                    HashMap<String, Integer> ionSeriesFrac = new HashMap<String, Integer>();
                    boolean[][] ionSeries = new boolean[4][spf.getPepLen() - 1];
                    double maxMatchedIonInt = 0, totMatchedIonInt = 0;
                    double fragError = 0, sumSquareError = 0;
                    double[] fragMassError = new double[mzHit.size()];
                    for (int j = 0; j < mzHit.size(); j++) {
                        MSMZHit mzHitTemp = mzHit.get(j);
                        double maxIonInte = Double.MIN_VALUE;
                        double massDiff = Double.MAX_VALUE;
                        int inteIndex = -1, terminate = mz.length - 1;
                        int mzCharge = mzHitTemp.MSMZHit_charge;
                        double mzTheo = mzHitTemp.MSMZHit_mz / (double) responseMzScale;
                        //System.out.println(mzHitTemp.MSMZHit_mz);
                        for (int k = 0; k < mz.length; k++) {
                            double mzTol = msSearch.MSSearchSettings_msmsppm ? (mz[k] - mzTheo) / (mzTheo - 1.007825) : (mz[k] - mzTheo) * mzCharge;
                            //System.out.println(mzTol);
                            if (mzTol < -msmsTol) {
                                continue;
                            } else if (Math.abs(mzTol) <= msmsTol) {
                                if (maxIonInte < intensity[k]) {
                                    maxIonInte = intensity[k];
                                    inteIndex = k;
                                    massDiff = mzTol;
                                }
                            } else if (mzTol > msmsTol) {
                                terminate = k;
                                break;
                            }
                        }
                        if (inteIndex == -1) {
                            for (int k = terminate; k >= 0; k--) {
                                double mzTol = msSearch.MSSearchSettings_msmsppm ? (mz[k] - mzTheo) / (mzTheo - 1.007825) : (mz[k] - mzTheo) * mzCharge;
                                if (Math.abs(massDiff) > Math.abs(mzTol)) {
                                    massDiff = mzTol;
                                    inteIndex = k;
                                } else if (mzTol < -msmsTol) {
                                    break;
                                }
                            }
                        }
                        fragMassError[j] = massDiff;
                        fragError += fragMassError[j];
                        sumSquareError += Math.pow(massDiff, 2);
                        int inte = abundList.get(inteIndex);
                        int seriesID = mzHitTemp.MSMZHit_ion.MSIonType;
                        int ionSub = mzHitTemp.MSMZHit_number;//begin with 0, end with peplength-2;
                        String key = seriesID + ";" + mzCharge;
                        if (inte > maxMatchedIonInt) {
                            maxMatchedIonInt = inte;
                        }
                        totMatchedIonInt += inte;
                        if (ionSeriesInt.containsKey(key)) {
                            ionSeriesInt.put(key, ionSeriesInt.get(key) + inte);
                        } else {
                            ionSeriesInt.put(key, (double) inte);
                        }
                        if (ionSeriesFrac.containsKey(key)) {
                            ionSeriesFrac.put(key, ionSeriesFrac.get(key) + 1);
                        } else {
                            ionSeriesFrac.put(key, 1);
                        }
                        //if(ionSub<=0 || ionSub>=mh.MSHits_pepstring.length())
                        //	continue;
                        if (seriesID == 1 && mzHitTemp.MSMZHit_charge == 1) {
                            ionSeries[0][ionSub] = true;
                        }
                        if (seriesID == 1 && mzHitTemp.MSMZHit_charge == 2) {
                            ionSeries[1][ionSub] = true;
                        }
                        if (seriesID == 4 && mzHitTemp.MSMZHit_charge == 1) {
                            ionSeries[0][mh.MSHits_pepstring.length() - ionSub - 2] = true;
                        }
                        if (seriesID == 4 && mzHitTemp.MSMZHit_charge == 2) {
                            ionSeries[1][mh.MSHits_pepstring.length() - ionSub - 2] = true;
                        }
                        if (seriesID == 2 && mzHitTemp.MSMZHit_charge == 1) {
                            ionSeries[2][ionSub] = true;
                        }
                        if (seriesID == 2 && mzHitTemp.MSMZHit_charge == 2) {
                            ionSeries[3][ionSub] = true;
                        }
                        if (seriesID == 5 && mzHitTemp.MSMZHit_charge == 1) {
                            ionSeries[2][mh.MSHits_pepstring.length() - ionSub - 2] = true;
                        }
                        if (seriesID == 5 && mzHitTemp.MSMZHit_charge == 2) {
                            ionSeries[3][mh.MSHits_pepstring.length() - ionSub - 2] = true;
                        }
                    }
                    for (int j = 0; j < ionSeries.length; j++) {
                        int longest = 1;
                        for (int k = 0; k < ionSeries[j].length; k++) {
                            longest = ionSeries[j][k] == false ? 1 : longest + 1;
                            spf.setLongestSeq(spf.getLongestSeq() < longest ? longest : spf.getLongestSeq());
                        }
                    }
                    for (Entry<String, Integer> entry : ionSeriesFrac.entrySet()) {
                        ionSeriesFrac.put(entry.getKey(), Integer.valueOf((int) Math.rint(100 * entry.getValue() / (double) (mh.MSHits_pepstring.length() - 1))));
                    }
                    spf.setIonSeriesFrac(ionSeriesFrac);
                    HashMap<String, Integer> outIonInte = new HashMap<String, Integer>();
                    for (Entry<String, Double> entry : ionSeriesInt.entrySet()) {
                        outIonInte.put(entry.getKey(), totMatchedIonInt == 0 ? 0 : Integer.valueOf((int) Math.rint(100 * (entry.getValue() / (double) totMatchedIonInt))));
                    }
                    spf.setOutIonInte(outIonInte);

                    spf.setIntMatchedTot(totMatchedIonInt == 0 ? 0 : Math.log(totMatchedIonInt / (double) abundScale));
                    spf.setRelIntMatchTot(totMatchedIonInt / (double) abundScale / Math.pow(Math.E, spf.getTotIne()));
                    spf.setMaxFragIonInt(maxMatchedIonInt / (double) abundScale);

                    DescriptiveStatistics ds = new DescriptiveStatistics(fragMassError);
                    double stdError = 0, calibCoef = 0;
                    calibCoef = mzHit.size() == 0 ? 0 : Math.pow(1 + Math.min(7, mzHit.size()), 2) / Math.pow(8, 2);
                    stdError = fragMassError.length == 0 ? 0 : calibCoef * ds.getStandardDeviation();
                    fragError = mzHit.size() == 0 ? 0 : calibCoef * fragError / mzHit.size();
                    spf.setFragError(fragError);
                    spf.setStdevFragError(stdError);
                    spf.setSumSquareFragError(sumSquareError * calibCoef);
                    spf.setFragDeltaMed(mzHit.size() == 0 ? 0 : ds.getPercentile(50));
                    spf.setFragDeltaMedPPM(spf.getFragDeltaMed() / spf.getMrCalc() * 1000000);
                    spf.setFragDeltaIqr(mzHit.size() == 0 ? 0 : (ds.getPercentile(75) - ds.getPercentile(25)));
                    spf.setFragDeltaIqrPPM(spf.getFragDeltaIqr() / spf.getMrCalc() * 1000000);
                    spf.setQmatch(omxFile.getPeptideToSpectrumMap().get(spf.getPepString()).size());
                    spf.setScan(specTemp.MSSpectrum_ids.MSSpectrum_ids_E.get(0));
                    String hitId = name + ":" + specNum + "." + (x + 1);
                    spf.setHitID(hitId);
                }
                omssa_feature_list.add(spf);
            }
        }
    }

    public void print_features() throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(new File(output_feature_file)));
        out.write(feature_title);
        out.newLine();
        out.flush();

        for (SpectraFeature spf : omssa_feature_list) {
            out.write(spf.getHitID() + "\t" + spf.getTargetDecoy() + "\t" + spf.getPepEvalue()
                    + "\t" + spf.getMrCalc() + "\t" + spf.getCharge() + "\t" + spf.getdM()
                    + "\t" + spf.getdMppm() + "\t" + spf.getAbsDM() + "\t" + spf.getAbsDMppm() + "\t" + spf.getIsoDM() + "\t" + spf.getIsoDMppm()
                    + "\t" + spf.getVarmodsRatio()
                    + "\t" + spf.getTotIne() + "\t" + spf.getIntMatchedTot() + "\t" + spf.getRelIntMatchTot() + "\t");
            if (us_ms2) {
                out.write(spf.getMaxFragIonInt() + "\t" + spf.getFragError() + "\t" + spf.getFragDeltaMed()
                        + "\t" + spf.getFragDeltaMedPPM() + "\t" + spf.getFragDeltaIqr() + "\t" + spf.getFragDeltaIqrPPM() + "\t");
            }
            out.write(spf.getLongestSeq() + "\t" + spf.getQmatch() + "\t" + spf.getEnzN() + "\t" + (spf.iscTermTryptic() ? "1" : "0")
                    + "\t" + (spf.isnTermTryptic() ? "1" : "0") + "\t" + spf.getPepLen() + "\t" + spf.getPvalue());
            for (int j = 0; j < ionSearch.size(); j++) {
                int type = ionSearch.get(j);
                for (int k = 1; k < 3; k++) {
                    HashMap<String, Integer> ionSeriesFrac = spf.getIonSeriesFrac();
                    HashMap<String, Integer> outIonInte = spf.getOutIonInte();
                    if (ionSeriesFrac.containsKey(type + ";" + k)) {
                        out.write("\t" + ionSeriesFrac.get(type + ";" + k) + "\t" + outIonInte.get(type + ";" + k));
                    } else {
                        out.write("\t0\t0");
                    }
                }
            }
            out.write("\t" + spf.getPepString() + "\t" + spf.getProteins());
            out.newLine();
        }
        out.flush();
        out.close();
        System.out.println("Parse OMSSA omx search result " + omssa_omx_file + " Done.");
    }

    private static int countVarMod(ArrayList<Modification> varMod, String pepSeq, String preAA, String postAA) {
        int modifiable = 0;
        for (int k = 0; k < varMod.size(); k++) {
            Modification amod = varMod.get(k);
            if (amod.getResidues() != null) {
                int index = 0;
                String[] residues = amod.getCompactResidues().split("");
                if (amod.isPepcterm() && !amod.isProcterm()) {
                    for (String str : residues) {
                        if (pepSeq.endsWith(str)) {
                            modifiable++;
                        }
                    }
                } else if (amod.isPepcterm() && amod.isProcterm() && postAA == null) {
                    for (String str : residues) {
                        if (pepSeq.endsWith(str)) {
                            modifiable++;
                        }
                    }
                }
                if (amod.isPepnterm() && !amod.isPronterm()) {
                    for (String str : residues) {
                        if (pepSeq.startsWith(str)) {
                            modifiable++;
                        }
                    }
                } else if (amod.isPepnterm() && amod.isPronterm() && preAA == null) {
                    for (String str : residues) {
                        if (pepSeq.startsWith(str)) {
                            modifiable++;
                        }
                    }
                }
                if (!amod.isPepcterm() && !amod.isPepnterm() && !amod.isProcterm() && !amod.isPronterm()) {
                    for (String str : residues) {
                        if (str == null || str.isEmpty()) {
                            continue;
                        }
                        while ((index = pepSeq.indexOf(str, index)) != -1) {
                            modifiable++;
                            index++;
                        }
                    }
                }
            } else {
                if (amod.isPepcterm() && !amod.isProcterm()) {
                    modifiable++;
                } else if (amod.isPepcterm() && amod.isProcterm() && postAA == null) {
                    modifiable++;
                }
                if (amod.isPepnterm() && !amod.isPronterm()) {
                    modifiable++;
                } else if (amod.isPepnterm() && amod.isPronterm() && preAA == null) {
                    modifiable++;
                }
            }
        }
        return modifiable;
    }

    /**
     * Check whether the site digested by trypsin.
     *
     * @param previous
     * @param next
     * @return
     */
    private static boolean isEnzymatic(String previous, String next) {
        if (next == null || previous == null
                || ((previous.toUpperCase().equals("K") || previous.toUpperCase().equals("R"))
                && !next.toUpperCase().equals("P"))) {
            return true;
        }
        return false;
    }

    private String realignPeptideInDB(String pepSeq,
            HashMap<String, Protein> proMap) {
        String pros = "";
        Iterator<String> iter = proMap.keySet().iterator();
        while (iter.hasNext()) {
            String proAcc = iter.next();
            String proSeq = proMap.get(proAcc).getProtein_sequence();
            if (proSeq.indexOf(pepSeq) != -1) {
                pros += proAcc + "\t";
                proMap.get(proAcc).getPeps().add(pepSeq);
            }
        }
        return pros;
    }

    private void clear() {
        omssa_feature_list.clear();
    }
}
