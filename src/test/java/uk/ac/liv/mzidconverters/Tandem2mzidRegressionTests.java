package uk.ac.liv.mzidconverters;

import java.io.File;
import java.util.List;

import junit.framework.TestCase;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DataCollection;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.converters.Tandem2mzid;

/**
 * In this class we use the jmzidentml parser (see
 * http://code.google.com/p/jmzidentml/) to read in the transformed mzIdentML
 * file and automatically verify its output based on a number of assertions.
 *
 * So: all automated regression tests can be placed here.
 *
 *
 * @author plukasse
 *
 */
public class Tandem2mzidRegressionTests extends TestCase {

    /**
     * Basic regression tests/checks
     *
     * @throws Exception
     */
    public void test_basic_55merge_tandem_file() throws Exception {
        String xTandemFile = "src/test/data/55merge_tandem.xml";
        String resultFile = "src/test/out/55merge_tandem.xml.mzid";
        new Tandem2mzid(xTandemFile, resultFile, "MS:1001348", "MS:1001062", false, true, MzIdentMLVersion.Version_1_1);

        //===================================================================================
        //========================= Checks /assertions section :=============================
        //===================================================================================
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(resultFile));

        DataCollection dc = unmarshaller.unmarshal(DataCollection.class);
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();
        assertEquals(1, sil.size());
        assertEquals(169, sil.get(0).getSpectrumIdentificationResult().size());

        int countSII = 0;
        int countPepEvid = 0;
        SpectrumIdentificationItem sIIWith2PepEvidence = null;
        SpectrumIdentificationItem sIIWith3PepEvidence = null;

        for (SpectrumIdentificationList sIdentList : sil) {
            for (SpectrumIdentificationResult spectrumIdentResult
                    : sIdentList.getSpectrumIdentificationResult()) {

                String spectrumID = spectrumIdentResult.getSpectrumID();

                for (SpectrumIdentificationItem spectrumIdentItem
                        : spectrumIdentResult.getSpectrumIdentificationItem()) {
                    countSII++;

                    if (spectrumID.equals("index=241")) {
                        sIIWith2PepEvidence = spectrumIdentItem;
                    }
                    if (spectrumID.equals("index=302")) {
                        sIIWith3PepEvidence = spectrumIdentItem;
                    }

                    countPepEvid += spectrumIdentItem.getPeptideEvidenceRef().size();

                    // If the auto-resolve mechanism is activated for SpectrumIdentificationItem
                    // then automatically resolve the Peptide Object; 
                    // (see  http://code.google.com/p/jmzidentml/wiki/AutoResolveExample)
                    if (MzIdentMLElement.SpectrumIdentificationItem.isAutoRefResolving()
                            && spectrumIdentItem.getPeptideRef() != null) {
                        Peptide peptide = spectrumIdentItem.getPeptide();
                        String peptideSequence = peptide.getPeptideSequence();
                        System.out.println("Pepetide Sequence = " + peptideSequence);
                    }

                } // end spectrum identification item
            } // end spectrum identification results
        }

        SequenceCollection sc = unmarshaller.unmarshal(SequenceCollection.class);

        //from the 169 identifications, we expect to find the following:
        // - one SIR that has 2 SIIs  (one mass peak matching to 2 different peptides)
        // - one SII (of another SIR) that has 2 PeptideEvidence items  (one of the peptides is found in 2 different proteins)
        //These are checked below.
        //Evidence for one of the SIRs having 2 SIIs
        assertEquals(170, countSII);

        //Evidence for one SII having 2 PeptideEvidence items, and one SII having 3 PeptideEvidence items:
        assertEquals(173, countPepEvid);
        assertEquals(2, sIIWith2PepEvidence.getPeptideEvidenceRef().size());
        assertEquals("GVGAER___", sIIWith2PepEvidence.getPeptideRef());
        assertEquals("GVGAER___", getPeptideEvidenceItem(sc, sIIWith2PepEvidence.getPeptideEvidenceRef().get(0)).getPeptideRef());
        assertEquals("GVGAER___", getPeptideEvidenceItem(sc, sIIWith2PepEvidence.getPeptideEvidenceRef().get(1)).getPeptideRef());

        //assert scores are set/translated correctly:
        //i.e. expect="4.0e+000" mh="588.310" delta="0.470" hyperscore="13.9"
        assertEquals("4.0", sIIWith2PepEvidence.getCvParam().get(0).getValue());
        assertEquals("X\\!Tandem:expect", sIIWith2PepEvidence.getCvParam().get(0).getName());
        assertEquals("MS:1001330", sIIWith2PepEvidence.getCvParam().get(0).getAccession());
        assertEquals("13.9", sIIWith2PepEvidence.getCvParam().get(1).getValue());
        assertEquals("X\\!Tandem:hyperscore", sIIWith2PepEvidence.getCvParam().get(1).getName());
        assertEquals("MS:1001331", sIIWith2PepEvidence.getCvParam().get(1).getAccession());

        //EEEEKGEEEK case (is a peptide matching 3 times on the same protein, at different places - is similar
        //to spectrum 6642 in test_modifications test case, so skipping refined tests for now):
        assertEquals("EEEEKGEEEK___", sIIWith3PepEvidence.getPeptideRef());
        assertEquals(3, sIIWith3PepEvidence.getPeptideEvidenceRef().size());

    }

    /**
     * Returns the PeptideEvidence object given the PeptideEvidenceRef item.
     *
     * @param sc
     * @param peptideEvidenceRef
     * @return
     */
    private PeptideEvidence getPeptideEvidenceItem(SequenceCollection sc, PeptideEvidenceRef peptideEvidenceRef) {
        return getPeptideEvidenceItem(sc, peptideEvidenceRef.getPeptideEvidenceRef());
    }

    private PeptideEvidence getPeptideEvidenceItem(SequenceCollection sc, String peptideEvidenceRefString) {
        List<PeptideEvidence> pepEvidenceList = sc.getPeptideEvidence();
        for (PeptideEvidence pepEvidence : pepEvidenceList) {
            if (pepEvidence.getId().equals(peptideEvidenceRefString)) {
                return pepEvidence;
            }
        }
        return null;
    }

    public void test_debug_File() throws Exception {
        String xTandemFileFromGalaxyStep = "src/test/data/DEBUG_MSMS_XTANDEM.bioml";
        String outMzid = "src/test/out/DEBUG_MSMS_XTANDEM.bioml.mzid";
        new Tandem2mzid(xTandemFileFromGalaxyStep, outMzid, "MS:1001348", "MS:1000584", false, true, MzIdentMLVersion.Version_1_1);

        //===================================================================================
        //========================= Checks /assertions section :=============================
        //===================================================================================
        //we expect TO find in the output the following items:
        // - 1 SIR
        // - 2 SII, each with 2 PeptideEvidence items. Sub-checks: 
        //    o the PeptideEvidence have the same sequence
        //    o the first SII is about evidence on  YIYEIAR
        //	  o the second SII is about evidence on YLYEIAR
        // - 4 dbsequences
        // - 2 peptides
        // - 4 peptideEvidence items
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(outMzid));

        DataCollection dc = unmarshaller.unmarshal(DataCollection.class);
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();
        assertEquals(1, sil.size());
        // - 1 SIR
        assertEquals(1, sil.get(0).getSpectrumIdentificationResult().size());

        int countSII = 0;
        int countPepEvid = 0;
        SpectrumIdentificationItem first_sII = sil.get(0).getSpectrumIdentificationResult().get(0).getSpectrumIdentificationItem().get(0);
        SpectrumIdentificationItem second_sII = sil.get(0).getSpectrumIdentificationResult().get(0).getSpectrumIdentificationItem().get(1);

        for (SpectrumIdentificationList sIdentList : sil) {
            for (SpectrumIdentificationResult spectrumIdentResult
                    : sIdentList.getSpectrumIdentificationResult()) {
                for (SpectrumIdentificationItem spectrumIdentItem
                        : spectrumIdentResult.getSpectrumIdentificationItem()) {
                    countSII++;
                    countPepEvid += spectrumIdentItem.getPeptideEvidenceRef().size();

                } // end spectrum identification item
            } // end spectrum identification results
        }

        SequenceCollection sc = unmarshaller.unmarshal(SequenceCollection.class);

        // - 2 SII, each with 2 PeptideEvidence items. 
        assertEquals(2, countSII);
        assertEquals(2, first_sII.getPeptideEvidenceRef().size());
        assertEquals(2, second_sII.getPeptideEvidenceRef().size());
        //Sub-checks: 
        //	o the PeptideEvidence items within a SII have the same sequence
        //	o the first SII is about evidence on  YIYEIAR
        //	o the second SII is about evidence on YLYEIAR

        assertEquals("YIYEIAR___", first_sII.getPeptideRef());
        assertEquals("YIYEIAR___", getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(0)).getPeptideRef());
        assertEquals("YIYEIAR___", getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(1)).getPeptideRef());

        assertEquals("YLYEIAR___", second_sII.getPeptideRef());
        assertEquals("YLYEIAR___", getPeptideEvidenceItem(sc, second_sII.getPeptideEvidenceRef().get(0)).getPeptideRef());
        assertEquals("YLYEIAR___", getPeptideEvidenceItem(sc, second_sII.getPeptideEvidenceRef().get(1)).getPeptideRef());

        // - 4 peptideEvidence items
        assertEquals(4, countPepEvid);

        // - 2 peptides
        assertEquals(2, sc.getPeptide().size());

        // - 4 dbsequences 
        assertEquals(4, sc.getDBSequence().size());

    }

    /**
     * Test case with special scenario that includes a peptide with
     * modifications found in 2 different proteins.
     *
     * This test also includes a check on the parsing of the ":reversed" flag
     * added by X!Tandem to the decoy matches.
     *
     * @throws Exception
     */
    public void test_modifications_and_decoyMatches() throws Exception {
        String xTandemFile = "src/test/data/S2_depl_spiked_6.mzML_xtandemOut.bioml";
        String resultFile = "src/test/out/S2_depl_spiked_6.mzML_xtandemOut.bioml.mzid";
        boolean isMs2SpectrumIdStartingAtZero = true;
        new Tandem2mzid(xTandemFile, resultFile, "MS:1001348", "MS:1000584", isMs2SpectrumIdStartingAtZero, true, MzIdentMLVersion.Version_1_1);

        //===================================================================================
        //========================= Checks /assertions section :=============================
        //===================================================================================
        //we expect TO find in the output the following items:
        // - 685 SIR
        // - when zooming in on spectrum 6642, we expect: 
        //		o 1 SII, with 2 PeptideEvidence items. Sub-checks: 
        //			> both the PeptideEvidence have the same sequence YGAEDMSGMEGSSGLGDRFGAK
        //			> the first PeptideEvidence matches region start="1531" end="1552"
        //			> the second PeptideEvidence matches region start="2035" end="2056"
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(resultFile));

        DataCollection dc = unmarshaller.unmarshal(DataCollection.class);
        AnalysisData ad = dc.getAnalysisData();

        // Get the list of SpectrumIdentification elements
        List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();
        assertEquals(1, sil.size());
        // - 685 SIR
        assertEquals(685, sil.get(0).getSpectrumIdentificationResult().size());

        SpectrumIdentificationResult sIR6642 = null;
        SpectrumIdentificationResult sIR3099 = null;

        for (SpectrumIdentificationList sIdentList : sil) {
            for (SpectrumIdentificationResult spectrumIdentResult
                    : sIdentList.getSpectrumIdentificationResult()) {
                String spectrumID = spectrumIdentResult.getSpectrumID();
                if (spectrumID.equals("index=6642")) {
                    sIR6642 = spectrumIdentResult;
                }
                if (spectrumID.equals("index=3099")) {
                    sIR3099 = spectrumIdentResult;
                }

            } // end spectrum identification results
        }

        // - when zooming in on spectrum 6642, we expect: 
        //		o 1 SII, with 2 PeptideEvidence items. Sub-checks: 
        //			> both the PeptideEvidence have the same sequence YGAEDMSGMEGSSGLGDRFGAK
        //			> the first PeptideEvidence matches region start="1531" end="1552"
        //			> the second PeptideEvidence matches region start="2035" end="2056"
        assertEquals(1, sIR6642.getSpectrumIdentificationItem().size());
        SpectrumIdentificationItem first_sII = sIR6642.getSpectrumIdentificationItem().get(0);
        assertEquals(2, first_sII.getPeptideEvidenceRef().size());

        //Sub-checks: 
        //			> both the PeptideEvidence have the same sequence YGAEDMSGMEGSSGLGDRFGAK, 
        //			  with the same modifications annotations at positions 6, 9, 12, 13, 1
        //			  that is: 15.9949@M$6;15.9949@M$9;79.9663@S$12;79.9663@S$13;79.9663@S$1;
        SequenceCollection sc = unmarshaller.unmarshal(SequenceCollection.class);
        String peptideSignature = "YGAEDMSGMEGSSGLGDRFGAK_15.9949@M$6;15.9949@M$9;79.9663@S$12;79.9663@S$13;79.9663@Y$1;__";
        assertEquals(peptideSignature, first_sII.getPeptideRef());
        assertEquals(peptideSignature, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(0)).getPeptideRef());
        assertEquals(peptideSignature, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(1)).getPeptideRef());
        //			> the first PeptideEvidence matches region start="1531" end="1552"
        assertEquals(1531, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(0)).getStart().intValue());
        assertEquals(1552, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(0)).getEnd().intValue());
        //			> the second PeptideEvidence matches region start="2035" end="2056"
        assertEquals(2035, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(1)).getStart().intValue());
        assertEquals(2056, getPeptideEvidenceItem(sc, first_sII.getPeptideEvidenceRef().get(1)).getEnd().intValue());

        //Here the case of AVPSVSSNVTAMSLQWNLIR in spectrum 3099, where it is matched 
        //with a modification on M and on its second <OR> third S !
        //What we expect is:
        // - there are 2 SII, one for each set of modifications (@M and @2ndS, @M and @3rdS)
        assertEquals(2, sIR3099.getSpectrumIdentificationItem().size());
        first_sII = sIR3099.getSpectrumIdentificationItem().get(0);
        // - each SII has 8 PeptideEvidence items (on 8 different proteins)
        assertEquals(8, first_sII.getPeptideEvidenceRef().size());
        SpectrumIdentificationItem second_sII = sIR3099.getSpectrumIdentificationItem().get(1);
        assertEquals(8, second_sII.getPeptideEvidenceRef().size());

        //Check decoy match. The DBSequence below should be linked to a PeptideEvidence that has isDecoy= true:
        //<DBSequence 
        //accession="IPI:IPI00008714.5|SWISS-PROT:O76027|ENSEMBL:ENSP00000357943|REFSEQ:NP_003559|H-INV:HIT000032464|
        //VEGA:OTTHUMP00000033027..." 
        //searchDatabase_ref="SearchDB_1" length="345" 
        //name="IPI:IPI00008714.5|SWISS-PROT:O76027|ENSEMBL:ENSP00000357943|REFSEQ:NP_003559|H-INV:HIT000032464|
        //VEGA:OTTHUMP00000033027 Tax_Id=9606 Gene_Symbol=ANXA9 Annexin A9:reversed" 
        //id="dbseq_IPI:IPI00008714.5|SWISS-PROT:O76027|ENSEMBL:ENSP00000357943|REFSEQ:NP_003559|H-INV:HIT000032464|
        //VEGA:OTTHUMP00000033027...">
        //<PeptideEvidence isDecoy="true" post="Q" pre="R" end="121" start="108" 
        //peptide_ref="QYQDFVRILHEPNR_-17.0265@Q$1;__" 
        //dBSequence_ref="dbseq_IPI:IPI00008714.5|SWISS-PROT:O76027|ENSEMBL:ENSP00000357943|REFSEQ:NP_003559|H-INV:HIT000032464|
        //VEGA:OTTHUMP00000033027..." id="PE653_2_737"/>
        PeptideEvidence decoyPepEvidence = getPeptideEvidenceItem(sc, "PE653_2_737");
        assertEquals(decoyPepEvidence.isIsDecoy(), true);

        //<PeptideEvidence isDecoy="false" post="R" pre="-" end="9" start="1" peptide_ref="MMIKETSLR___" 
        //dBSequence_ref="dbseq_IPI:IPI00218930.1|SWISS-PROT:Q06190-2|TREMBL:B4DS28|ENSEMBL:ENSP00000334748|
        //REFSEQ:NP_871626|VEGA:OTTHUMP00000216421..." id="PE679_2_800"/>
        decoyPepEvidence = getPeptideEvidenceItem(sc, "PE679_2_800");
        assertEquals(decoyPepEvidence.isIsDecoy(), false);

    }

    public void test_constructors() throws Exception {
        //the simple constructor first, which forces databaseFileFormatID and massSpecFileFormatID inference
        //based on extension names found in file:
        String xTandemFile = "src/test/data/55merge_tandem.xml";
        String resultFile = "src/test/out/55merge_tandem.xml.mzid";
        new Tandem2mzid(xTandemFile, resultFile, MzIdentMLVersion.Version_1_1);

        //===================================================================================
        //========================= Checks /assertions section :=============================
        //===================================================================================
        //the following lines are found in the xtandem file:
        //
        //<note type="input" label="spectrum, path">D:/TestSpace/NeoTestMarch2011/55merge.mgf</note>
        //<note label="list path, sequence source #1">D:/Software/Databases/Neospora_3rndTryp/Neo_rndTryp_3times.fasta.pro</note>
        //
        //so we expect the algorithm to set the following in the mzid file:
        //
        //"MS:1001348";
        //"FASTA format"
        //in ../DataCollection/Inputs/SearchDatabase/FileFormat/cvParam accession="MS:1001348" cvRef="PSI-MS" name="FASTA format"
        //
        //"MS:1001062";
        //"Mascot MGF file";
        //in ../DataCollection/Inputs/SpectraData/FileFormat/cvParam accession="MS:1001062" cvRef="PSI-MS" name="Mascot MGF file"
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(resultFile));

        DataCollection dc = unmarshaller.unmarshal(DataCollection.class);
        //check Fasta format:
        CvParam searchDBCvParam = dc.getInputs().getSearchDatabase().get(0).getFileFormat().getCvParam();
        assertEquals(searchDBCvParam.getAccession(), "MS:1001348");
        //check mascot format:
        CvParam spectraDataCvParam = dc.getInputs().getSpectraData().get(0).getFileFormat().getCvParam();
        assertEquals(spectraDataCvParam.getAccession(), "MS:1001062");

        //Overrule with custom codes and check:
        resultFile = resultFile + "_dummy.mzid";
        new Tandem2mzid(xTandemFile, resultFile, "MS:1001349", "MS:1000566", false, true, MzIdentMLVersion.Version_1_1);

        unmarshaller = new MzIdentMLUnmarshaller(new File(resultFile));

        dc = unmarshaller.unmarshal(DataCollection.class);
        //check Fasta format:
        searchDBCvParam = dc.getInputs().getSearchDatabase().get(0).getFileFormat().getCvParam();
        assertEquals(searchDBCvParam.getAccession(), "MS:1001349");
        //check mascot format:
        spectraDataCvParam = dc.getInputs().getSpectraData().get(0).getFileFormat().getCvParam();
        assertEquals(spectraDataCvParam.getAccession(), "MS:1000566");

    }

    public void test_NTERM_modifications() throws Exception {
        String xTandemFile = "src/test/data/55merge_tandem_NTERM_MODIF.xml";
        String resultFile = "src/test/out/55merge_tandem_NTERM_MODIF.xml.mzid";
        new Tandem2mzid(xTandemFile, resultFile, "MS:1001348", "MS:1001062", false, true, MzIdentMLVersion.Version_1_1);

        //===================================================================================
        //========================= Checks /assertions section :=============================
        //===================================================================================
        //we expect TO find in the output the following items:
        // - a Peptide item with PeptideSequence = LCYIALDFDEEMKAAEDSSDIEK and modification name = Propionyl:13C(3)
        //   at location 0. 
        MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(new File(resultFile));

        // Get the list of Peptides
        SequenceCollection sc = unmarshaller.unmarshal(SequenceCollection.class);
        List<Peptide> peptideList = sc.getPeptide();

        //the peptide we are looking for is the first one:
        Peptide firstPeptide = peptideList.get(0);
        assertEquals(firstPeptide.getPeptideSequence(), "LCYIALDFDEEMKAAEDSSDIEK");
        Modification firstModification = firstPeptide.getModification().get(0);
        assertEquals(firstModification.getLocation(), new Integer(0));
        assertEquals(firstModification.getMonoisotopicMassDelta(), 59.036279);
        assertEquals(firstModification.getCvParam().size(), 1);
        assertEquals(firstModification.getCvParam().get(0).getName(), "Propionyl:13C(3)");

    }

    /*	
     public void test_galaxyFile() throws Exception
     {
     String xTandemFileFromGalaxyStep = "src/test/data/Galaxy-[X_Tandem_on_data_N].bioml";
     new Tandem2mzid(xTandemFileFromGalaxyStep, xTandemFileFromGalaxyStep + ".mzid", true);
     }

	
     public void test_multiple_sii_case() throws Exception
     {
     String xTandemFileFromGalaxyStep = "src/test/data/20091015_spiking_0_05fmol_CytC_MSMS_XTANDEM.out";
     String outMzid = xTandemFileFromGalaxyStep + ".mzid";
     new Tandem2mzid(xTandemFileFromGalaxyStep, outMzid, true);
     }
	
	
	
     */
}
