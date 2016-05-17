package bgi.ipeak.percolator.paser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.FragmentArray;
import uk.ac.ebi.jmzidml.model.mzidml.IonType;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import bgi.ipeak.percolator.util.PSMImpl;
import bgi.ipeak.percolator.util.Parameters;
import bgi.ipeak.percolator.util.PeptideInfor;
import java.util.Vector;
import javax.xml.bind.JAXBException;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;

/**
 * @author Administrator
 *
 */
public class MSGFParser {

    private Parameters param;
    private String input_mzid_file;
    private String output_feature_file;
    static private Double neutron = 1.0033548378;
    static private Boolean use_ms2feature = false;
    private Vector<MSGFFeatures> fList = new Vector<MSGFFeatures>();
    private String dr = "###";

    public MSGFParser(String search_result, String out_put_feature_file, String decoyRegex) {
        this.input_mzid_file = search_result;
        this.output_feature_file = out_put_feature_file;
        this.dr = decoyRegex;
    }

    public MSGFParser() {

    }

    public static void setUse_ms2feature(Boolean use_ms2feature) {
        MSGFParser.use_ms2feature = use_ms2feature;
    }

    public void get_feature_file() throws IOException, MalformedURLException, JAXBException {
        parseMSGFMzid();
        writeFeatures();
        clear();
    }

    /**
     * Method to parse MSGF+ search result in mzid.
     *
     * @param fList
     * @param srcFile
     * @throws MalformedURLException
     */
    private void parseMSGFMzid() throws MalformedURLException, JAXBException {

        System.out.println("Parse MSGF+ mzid search result: " + input_mzid_file + "...");

        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(input_mzid_file));
        this.param = new Parameters();

        SequenceCollection sc = unmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);
        //parse search parameters from <SpectrumIdentificationProtocol/>.
        //HashMap<String, Boolean> modFixed=new HashMap<String, Boolean>();
        AnalysisProtocolCollection apc = unmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
        SpectrumIdentificationProtocol sp = apc.getSpectrumIdentificationProtocol().get(0);
        for (int i = 0; i < sp.getAdditionalSearchParams().getCvParam().size(); i++) {
            CvParam cvTmp = sp.getAdditionalSearchParams().getCvParam().get(i);
            if (cvTmp.getAccession().equals("MS:1001211") || cvTmp.getName().equals("parent mass type mono")) {
                this.param.setParentMassType("Monoisotopic");
            }
            if (cvTmp.getAccession().equals("MS:1001256") || cvTmp.getName().equals("fragment mass type mono")) {
                this.param.setFragmentMassType("Monoisotopic");
            }
            if (cvTmp.getAccession().equals("MS:1001212") || cvTmp.getName().equals("parent mass type average")) {
                this.param.setParentMassType("Average");
            }
            if (cvTmp.getAccession().equals("MS:1001255") || cvTmp.getName().equals("fragment mass type average")) {
                this.param.setFragmentMassType("Average");
            }
            //System.out.println(cvTmp.getAccession()+"\t"+cvTmp.getName());
        }

        this.param.setParentTolerance(sp.getParentTolerance()
                != null ? Double.parseDouble(sp.getParentTolerance().getCvParam().get(0).getValue()) : 0.0);
        this.param.setParentTolUnitPPM(sp.getParentTolerance()
                == null ? false : sp.getParentTolerance().getCvParam().get(0).getAccession().equals("UO:0000221") ? false : true);

//		this.param.setFragmentTolUnitPPM(sp.getFragmentTolerance()
//				==null?false:sp.getFragmentTolerance().getCvParam().get(0).getAccession().equals("UO:0000221")?false:true);
//		this.param.setFragmentTolerance(sp.getFragmentTolerance()
//				==null?0:Double.parseDouble(sp.getFragmentTolerance().getCvParam().get(0).getValue()));
        //parse PSM information from <SpectrumIdentificationResult/>.
        //parse features from <SpectrumIdentificationItem/>.
        Iterator<SpectrumIdentificationResult> psmResult = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
        Pattern pattern = null;
        if (this.dr != null) {
            pattern = Pattern.compile(this.dr);
        }
        Matcher matcher = null;
        while (psmResult.hasNext()) {
            MSGFFeatures mf = new MSGFFeatures();
            SpectrumIdentificationResult sir = psmResult.next();
            PSMImpl psmTmp = new PSMImpl();
            psmTmp.setPSMId(Integer.parseInt(sir.getSpectrumID().split("=")[1]));
            for (int i = 0; i < sir.getCvParam().size(); i++) {
                if (sir.getCvParam().get(i).getAccession().equals("MS:1000796")) {
                    mf.setTitle(sir.getCvParam().get(i).getValue());
                }
            }
            for (int i = 0; i < sir.getSpectrumIdentificationItem().size(); i++) {
                SpectrumIdentificationItem siiTmp = sir.getSpectrumIdentificationItem().get(i);
                if (siiTmp.getRank() != 1) {
                    continue;
                }
                psmTmp.setRank(1);
                PeptideInfor iPep = new PeptideInfor();
                Peptide peptide = unmarshaller.unmarshal(Peptide.class, siiTmp.getPeptideRef());
                String peptideSequence = peptide.getPeptideSequence();
                iPep.setPeptideSequence(peptideSequence);
                iPep.setPeptideLength(peptideSequence.length());
                iPep.setProteins(new ArrayList<String>());
                PeptideEvidence peptideEvidence = unmarshaller.unmarshal(PeptideEvidence.class, siiTmp.getPeptideEvidenceRef().get(0).getPeptideEvidenceRef()); //用于后面酶切的判断等
                for (PeptideEvidenceRef peptideEvidenceRef : siiTmp.getPeptideEvidenceRef()) {
                    PeptideEvidence pEvidence = unmarshaller.unmarshal(PeptideEvidence.class, peptideEvidenceRef.getPeptideEvidenceRef());
                    //System.out.println(peptideEvidenceRef.getPeptideEvidence().getPre());
                    DBSequence dbSequence = unmarshaller.unmarshal(DBSequence.class, pEvidence.getDBSequenceRef());
                    String acc = dbSequence.getAccession();
                    iPep.getProteins().add(acc);
                    iPep.setPostAA(pEvidence.getPost());
                    iPep.setPreAA(pEvidence.getPre());
                    matcher = pattern.matcher(acc);
                    if (matcher.find()) {
                        iPep.setDecoy(true);
                        //break;
                    }
                }
                //System.out.println(iPep.getPeptideSequence());
                iPep.setCharge(Byte.parseByte(String.valueOf(siiTmp.getChargeState())));
                iPep.setExperimentalMass(siiTmp.getExperimentalMassToCharge());
                iPep.setTheoriticalMass(siiTmp.getCalculatedMassToCharge());

                mf.setCharge(iPep.getCharge());
                mf.setIndex(psmTmp.getPSMId());
                String fileName = new File(input_mzid_file).getName();
                mf.setFileName(fileName.substring(0, fileName.indexOf(".")));
                psmTmp.setPeptide(iPep);
                psmTmp.setTarget(iPep.isDecoy() ? false : true);
                ArrayList<bgi.ipeak.percolator.util.FragmentIon> fion = new ArrayList<bgi.ipeak.percolator.util.FragmentIon>();
                for (int j = 0; j < siiTmp.getFragmentation().getIonType().size(); j++) {
                    IonType ion = siiTmp.getFragmentation().getIonType().get(j);
                    for (int k = 0; k < ion.getIndex().size(); k++) {
                        bgi.ipeak.percolator.util.FragmentIon fragIon = new bgi.ipeak.percolator.util.FragmentIon();
                        fragIon.setDoubleCharged(ion.getCharge() > 1 ? true : false);
                        fragIon.setNumber(ion.getIndex().get(k));
                        for (int l = 0; l < ion.getFragmentArray().size(); l++) {
                            FragmentArray fay = ion.getFragmentArray().get(l);
                            if (fay.getMeasureRef().equals("m_mz")) {
                                fragIon.setMz(fay.getValues().get(k));
                            }
                            if (fay.getMeasureRef().equals("m_intensity")) {
                                fragIon.setIntensity(fay.getValues().get(k));
                            }
                            if (fay.getMeasureRef().equals("m_error")) {
                                fragIon.setMassError(fay.getValues().get(k));
                            }
                        }
                        fragIon.setType(ion.getCvParam().getName().split(":")[1]);
                        fragIon.setIDAccordeingToType();
                        fion.add(fragIon);
                    }
                }
                for (int j = 0; j < siiTmp.getCvParam().size(); j++) {
                    CvParam cv = siiTmp.getCvParam().get(j);
                    if (cv.getAccession().equals("MS:1002049")) {
                        mf.setRawScore(Double.parseDouble(cv.getValue()));
                    }
                    if (cv.getAccession().equals("MS:1002050")) {
                        mf.setDeNovoScore(Double.parseDouble(cv.getValue()));
                    }
                    if (cv.getAccession().equals("MS:1002052")) {
                        mf.setLnSpecEvalue(-Math.log(Double.parseDouble(cv.getValue())));
                    }
                    if (cv.getAccession().equals("MS:1002053")) {
                        mf.setLnEvalue(-Math.log(Double.parseDouble(cv.getValue())));
                    }
                }
                mf.setEnergy(mf.getRawScore() - mf.getDeNovoScore());
                mf.setScoreRatio(mf.getRawScore() / (mf.getDeNovoScore() + 0.0001));
                for (int j = 0; j < siiTmp.getUserParam().size(); j++) {
                    UserParam up = siiTmp.getUserParam().get(j);
                    if (up.getName().equals("IsotopeError")) {
                        mf.setIsotopeError(Double.parseDouble(up.getValue()));
                    }
                    if (up.getName().equals("ExplainedIonCurrentRatio")) {
                        mf.setLnExplainedIonCurrent(Math.log(Double.parseDouble(up.getValue()) + 0.0001));//+0.0001 ��Ϊ�˷�ֹ�� 0 ȡ log �������
                    }
                    if (up.getName().equals("NTermIonCurrentRatio")) {
                        mf.setLnNTermIonCurrentRatio(Math.log(Double.parseDouble(up.getValue()) + 0.0001));//+0.0001 ��Ϊ�˷�ֹ�� 0 ȡ log �������
                    }
                    if (up.getName().equals("CTermIonCurrentRatio")) {
                        mf.setLnCTermIonCurrentRatio(Math.log(Double.parseDouble(up.getValue()) + 0.0001));//+0.0001 ��Ϊ�˷�ֹ�� 0 ȡ log �������
                    }
                    if (up.getName().equals("MS2IonCurrent")) {
                        mf.setLnMS2IonCurrent(Math.log(Double.parseDouble(up.getValue())));
                    }
                    //MeanErrorTop7 StdevErrorTop7 ���������������������޷�����õ���
                    if (up.getName().equals("MeanErrorTop7")) {
                        mf.setMeanErrorTop7(Double.parseDouble(up.getValue()));
                    }
                    if (up.getName().equals("StdevErrorTop7")) {
                        mf.setStdevErrorTop7(Double.parseDouble(up.getValue()));
                    }
                }
                double dM = 1000000 * (siiTmp.getExperimentalMassToCharge() - siiTmp.getCalculatedMassToCharge()
                        - mf.getIsotopeError() * neutron / mf.getCharge()) / (siiTmp.getCalculatedMassToCharge()
                        + mf.getIsotopeError() * neutron / mf.getCharge());
                mf.setDeltaMass(dM);
                mf.setAbsDeltaMass(Math.abs(mf.getDeltaMass()));
                mf.setSqMeanErrorTop7(Math.sqrt(mf.getMeanErrorTop7()));

                String pepSeq = iPep.getPeptideSequence();
                if (pepSeq.startsWith("K") || pepSeq.startsWith("R")) {
                    mf.setEnzN(true);
                }
                if (pepSeq.endsWith("K") || pepSeq.endsWith("R")) {
                    mf.setEnzC(true);
                    mf.setEnzN(true);
                }
                int enzInt = 0, index = 0;
                while ((index = pepSeq.indexOf("R|K", index)) != -1) {
                    enzInt++;
                }
                mf.setEnzInt(enzInt);
                mf.setiPSM(psmTmp);
                fList.add(mf);
            }
        }
    }

    private void writeFeatures() throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(output_feature_file));
        String title = "PSMid\tTargetDecoy\tRawScore\tDeNovoScore\tScoreRatio\tEnergy\tlnEvalue\tlnSpecEvalue"
                + "\tIsotopeError\tlnExplainedIonCurrent\tlnNTermIonCurrentRatio\tlnCTermIonCurrectRatio"
                + "\tlnMS2IonCurrent\tMass\tPepLen\tDeltaMass(ppm)\tabsDeltaMass(ppm)\tMeanErrorTop7\tsqMeanErrorTop7"
                + "\tStdevErrorTop7\tCharge\tPeptide\tProteins\n";
        if (!use_ms2feature) {
            title = "PSMid\tTargetDecoy\tRawScore\tDeNovoScore\tScoreRatio\tEnergy\tlnEvalue\tlnSpecEvalue"
                    + "\tIsotopeError\tlnExplainedIonCurrent\tlnNTermIonCurrentRatio\tlnCTermIonCurrectRatio"
                    + "\tlnMS2IonCurrent\tMass\tPepLen\tDeltaMass(ppm)\tabsDeltaMass(ppm)"
                    + "\tCharge\tPeptide\tProteins\n";
        }
        out.write(title);
        out.flush();
        for (int i = 0; i < fList.size(); i++) {
            MSGFFeatures mf = fList.get(i);
            PSMImpl psm = mf.getiPSM();
            String psmid = mf.getFileName() + ":" + mf.getIndex() + "." + psm.getRank();
            String feature_infor;
            if (use_ms2feature) {
                feature_infor = psmid + "\t" + (psm.isTarget() ? 1 : -1) + "\t" + mf.getRawScore() + "\t" + mf.getDeNovoScore() + "\t"
                        + mf.getScoreRatio() + "\t" + mf.getEnergy() + "\t" + mf.getLnEvalue() + "\t" + mf.getLnSpecEvalue() + "\t"
                        + mf.getIsotopeError() + "\t" + mf.getLnExplainedIonCurrent() + "\t" + mf.getLnNTermIonCurrentRatio() + "\t"
                        + mf.getLnCTermIonCurrentRatio() + "\t" + mf.getLnMS2IonCurrent() + "\t"
                        + psm.getPeptide().getTheoriticalMass() + "\t" + psm.getPeptide().getPeptideLength() + "\t"
                        + mf.getDeltaMass() + "\t" + mf.getAbsDeltaMass() + "\t" + mf.getMeanErrorTop7() + "\t" + mf.getSqMeanErrorTop7() + "\t"
                        + mf.getStdevErrorTop7() + "\t" + mf.getCharge() + "\t" + psm.getPeptide().getPeptideSequence();
            } else {
                feature_infor = psmid + "\t" + (psm.isTarget() ? 1 : -1) + "\t" + mf.getRawScore() + "\t" + mf.getDeNovoScore() + "\t"
                        + mf.getScoreRatio() + "\t" + mf.getEnergy() + "\t" + mf.getLnEvalue() + "\t" + mf.getLnSpecEvalue() + "\t"
                        + mf.getIsotopeError() + "\t" + mf.getLnExplainedIonCurrent() + "\t" + mf.getLnNTermIonCurrentRatio() + "\t"
                        + mf.getLnCTermIonCurrentRatio() + "\t" + mf.getLnMS2IonCurrent() + "\t"
                        + psm.getPeptide().getTheoriticalMass() + "\t" + psm.getPeptide().getPeptideLength() + "\t"
                        + mf.getDeltaMass() + "\t" + mf.getAbsDeltaMass() + "\t" + mf.getCharge() + "\t" + psm.getPeptide().getPeptideSequence();
            }
            out.write(feature_infor);
            for (int j = 0; j < psm.getPeptide().getProteins().size(); j++) {
                out.write("\t" + psm.getPeptide().getProteins().get(j));
            }
            out.newLine();
            out.flush();
        }
        out.close();
    }

    private void clear() {
        fList.clear();
    }
}
