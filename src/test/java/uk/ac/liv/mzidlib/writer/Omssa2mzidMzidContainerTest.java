/*
 * Date: 06-Sep-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.Omssa2mzidMzidContainerTest.java
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package uk.ac.liv.mzidlib.writer;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;
import static junit.framework.TestCase.fail;

import java.io.File;
import java.util.List;

import de.proteinms.omxparser.util.MSSpectrum;

import org.apache.commons.io.FileUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidenceRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabase;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SourceFile;
import uk.ac.ebi.jmzidml.model.mzidml.SpectraData;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.model.mzidml.UserParam;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;

/**
 *
 * @author Da Qi
 */
public class Omssa2mzidMzidContainerTest {

    private static Omssa2mzidMzidContainer testContainer;
    private final double diff = 0.000001;
    private static File testFile;

    public Omssa2mzidMzidContainerTest() {

    }

    @BeforeClass
    public static void setUpClass() {
        testFile = FileUtils.getFile("src", "test", "data",
                                     "55merge_omssa.omx");

        testContainer = new Omssa2mzidMzidContainer(
                testFile.getAbsolutePath(), Boolean.FALSE, "REV_", null, null,
                MzIdentMLVersion.Version_1_1);
    }

    @AfterClass
    public static void tearDownClass() {
    }

    /**
     * Test of getAnalysisCollection method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetAnalysisCollection() {
        System.out.println("getAnalysisCollection");
        AnalysisCollection analysisCollection = testContainer
                .getAnalysisCollection();
        assertNotNull("AnalysisCollection is null.", analysisCollection);

        assertNull(analysisCollection.getProteinDetection());
        List<SpectrumIdentification> specIdentList = analysisCollection
                .getSpectrumIdentification();
        assertEquals(specIdentList.size(), 1);

        SpectrumIdentification specIdent = specIdentList.get(0);
        assertEquals(specIdent.getSpectrumIdentificationProtocolRef(),
                     "SearchProtocol_1");
        assertEquals(specIdent.getSpectrumIdentificationListRef(), "SI_List_1");
        assertEquals(specIdent.getId(), "SpecIdent_1");

        List<SearchDatabaseRef> searchDbRefList = specIdent
                .getSearchDatabaseRef();
        assertEquals(searchDbRefList.size(), 1);
        SearchDatabaseRef searchDbRef = searchDbRefList.get(0);
        assertEquals(searchDbRef.getSearchDatabaseRef(), "SearchDB_1");
    }

    /**
     * Test of getAnalysisProtocolCollection method, of class
     * Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetAnalysisProtocolCollection() {
        System.out.println("getAnalysisProtocolCollection");

        AnalysisProtocolCollection apc = testContainer
                .getAnalysisProtocolCollection();

        assertNull(apc.getProteinDetectionProtocol());
        List<SpectrumIdentificationProtocol> sipList = apc
                .getSpectrumIdentificationProtocol();
        assertEquals(sipList.size(), 1);

        SpectrumIdentificationProtocol sip = sipList.get(0);
        assertEquals(sip.getAnalysisSoftwareRef(), "ID_software");
        assertEquals(sip.getId(), "SearchProtocol_1");

        //test SearchType
        assertEquals(sip.getSearchType().getCvParam(), CvConstants.MS_MS_SEARCH);

        //test AdditionalSearchParams
        assertEquals(sip.getAdditionalSearchParams().getCvParam().size(), 2);
        assertTrue(sip.getAdditionalSearchParams().getCvParam().contains(
                CvConstants.FRAGMENT_MASS_TYPE_MONO));
        assertTrue(sip.getMassTable().isEmpty());

        //test ModificationParams
        assertEquals(sip.getModificationParams().getSearchModification().size(),
                     2);

        //test SearchModification
        SearchModification searchMod = sip.getModificationParams()
                .getSearchModification().get(1);
        assertFalse(searchMod.isFixedMod());
        assertEquals(15.994915, searchMod.getMassDelta(), diff);
        assertEquals(searchMod.getResidues().size(), 1);
        assertEquals(searchMod.getResidues().get(0), "M");
        CvParam cp = searchMod.getCvParam().get(0);
        assertEquals(cp.getCvRef(), "UNIMOD");
        assertEquals(cp.getCv(), CvConstants.UNIMOD_CV);
        assertEquals(cp.getAccession(), "UNIMOD:35");
        assertEquals(cp.getName(), "Oxidation");

        //test Enzymes
        assertEquals(sip.getEnzymes().getEnzyme().size(), 1);
        assertFalse(sip.getEnzymes().isIndependent());

        //test Enzyme
        Enzyme en = sip.getEnzymes().getEnzyme().get(0);
        assertEquals(en.getCTermGain(), "OH");
        assertEquals(en.getNTermGain(), "H");
        assertFalse(en.isSemiSpecific());
        assertEquals(en.getMissedCleavages().intValue(), 1);
        assertEquals(en.getEnzymeName().getCvParam().get(0).getAccession(),
                     "MS:1001251");
        assertEquals(en.getEnzymeName().getCvParam().get(0).getName(), "Trypsin");

        //test FragmentTolerance
        assertEquals(sip.getFragmentTolerance().getCvParam().size(), 2);
        cp = sip.getFragmentTolerance().getCvParam().get(0);
        assertEquals(cp.getUnitCvRef(), "UO");
        assertEquals(cp.getUnitAccession(), "UO:0000221");
        assertEquals(cp.getUnitName(), "dalton");
        assertEquals(cp.getAccession(), "MS:1001412");
        assertEquals(cp.getValue(), "0.8");

        //test ParentTolerance
        assertEquals(sip.getParentTolerance().getCvParam().size(), 2);
        assertEquals(sip.getParentTolerance().getCvParam().get(0).getName(),
                     "search tolerance plus value");

        //test Threshold
        assertEquals(sip.getThreshold().getCvParam().get(0).getAccession(),
                     "MS:1001494");
    }

    /**
     * Test of getAnalysisSampleCollection method, of class
     * Omssa2mzidMzidContainer.
     */
//    @Test
//    public void testGetAnalysisSampleCollection() {
//        System.out.println("getAnalysisSampleCollection");
//        Omssa2mzidMzidContainer instance = null;
//        AnalysisSampleCollection expResult = null;
//        AnalysisSampleCollection result = instance.getAnalysisSampleCollection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    /**
     * Test of getAnalysisSoftwareList method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetAnalysisSoftwareList() {
        System.out.println("getAnalysisSoftwareList");

        AnalysisSoftwareList softwareList = testContainer
                .getAnalysisSoftwareList();
        //TODO        
        //assertEquals(softwareList.getAnalysisSoftware().size(), 2); 
        AnalysisSoftware software = softwareList.getAnalysisSoftware().get(0);
        assertEquals(software.getName(), "OMSSA");
        assertEquals(software.getSoftwareName().getCvParam(), CvConstants.OMSSA);
    }

    /**
     * Test of getAuditCollection method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetAuditCollection() {
        System.out.println("getAuditCollection");

        AuditCollection auditCol = testContainer.getAuditCollection();
        assertEquals(auditCol.getPerson().size(), 1);
        assertEquals(auditCol.getOrganization().size(), 1);
        assertEquals(auditCol.getPerson().get(0).getLastName(), "secondName");
        assertEquals(auditCol.getPerson().get(0).getAffiliation().get(0)
                .getOrganization(), auditCol.getOrganization().get(0));
    }

    /**
     * Test of getBibliographicReference method, of class
     * Omssa2mzidMzidContainer.
     */
//    @Test
//    public void testGetBibliographicReference() {
//        System.out.println("getBibliographicReference");
//        BibliographicReference bib = testContainer.getBibliographicReference();
//        assertNull(bib);
//    }
    /**
     * Test of getCvList method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetCvList() {
        System.out.println("getCvList");

        CvList cvList = testContainer.getCvList();
        assertTrue(cvList.getCv().contains(CvConstants.PSI_CV));
        assertTrue(cvList.getCv().contains(CvConstants.UNIT_CV));
        assertTrue(cvList.getCv().contains(CvConstants.UNIMOD_CV));
    }

    /**
     * Test of getInputs method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetInputs() {
        System.out.println("getInputs");

        Inputs inputs = testContainer.getInputs();

        //test SourceFile
        List<SourceFile> srcFileList = inputs.getSourceFile();
        assertEquals(srcFileList.size(), 1);
        SourceFile srcFile = srcFileList.get(0);
        assertEquals(srcFile.getLocation(), testFile.getAbsolutePath());
        assertEquals(srcFile.getFileFormat().getCvParam().getAccession(),
                     "MS:1001400");
        assertEquals(srcFile.getFileFormat().getCvParam().getName(),
                     "OMSSA xml file");

        //test SearchDatabase
        List<SearchDatabase> schDbList = inputs.getSearchDatabase();
        assertEquals(schDbList.size(), 1);
        SearchDatabase schDb = schDbList.get(0);
        assertTrue(schDb.getNumDatabaseSequences() == 22348);
        assertEquals(schDb.getLocation(),
                     "D:/Software/Databases/Neospora_3rndTryp/Neo_rndTryp_3times.fasta");
        CvParam cpDbName = schDb.getDatabaseName().getCvParam();
        assertNull(cpDbName);
        UserParam upDbName = schDb.getDatabaseName().getUserParam();
        assertEquals(upDbName.getName(),
                     "D:/Software/Databases/Neospora_3rndTryp/Neo_rndTryp_3times.fasta");

        //test SpectraData
        List<SpectraData> spcDataList = inputs.getSpectraData();
        assertEquals(spcDataList.size(), 1);
        SpectraData spcData = spcDataList.get(0);
        assertEquals(spcData.getLocation(),
                     "D:/TestSpace/NeoTestMarch2011/55merge.mgf");
        assertEquals(spcData.getFileFormat().getCvParam().getAccession(),
                     CvConstants.MASCOT_MGF_FORMAT.getAccession());
        assertEquals(spcData.getSpectrumIDFormat().getCvParam().getName(),
                     "multiple peak list nativeID format");
        assertEquals(spcData.getSpectrumIDFormat().getCvParam().getAccession(),
                     "MS:1000774");
    }

    /**
     * Test of getProteinDetectionList method, of class Omssa2mzidMzidContainer.
     */
//    @Test
//    public void testGetProteinDetectionList() {
//        System.out.println("getProteinDetectionList");
//        Omssa2mzidMzidContainer instance = null;
//        ProteinDetectionList expResult = null;
//        ProteinDetectionList result = instance.getProteinDetectionList();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    /**
     * Test of getProvider method, of class Omssa2mzidMzidContainer.
     */
//    @Test
//    public void testGetProvider() {
//        System.out.println("getProvider");
//        Omssa2mzidMzidContainer instance = null;
//        Provider expResult = null;
//        Provider result = instance.getProvider();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    /**
     * Test of getSequenceCollection method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetSequenceCollection() {
        System.out.println("getSequenceCollection");

        SequenceCollection seqCol = testContainer.getSequenceCollection();

        List<DBSequence> dbSeqList = seqCol.getDBSequence();
        assertEquals(dbSeqList.size(), 66);
        for (DBSequence dbSeq : dbSeqList) {
            if (dbSeq.getId().equals("dbseq_Rnd3psu|NC_LIV_135580")) {
                assertEquals(dbSeq.getLength(), new Integer(1893));
                assertEquals(dbSeq.getAccession(), "Rnd3psu|NC_LIV_135580");
                assertEquals(dbSeq.getCvParam().get(0).getValue(),
                             "Rnd3psu|NC_LIV_135580 Rnd3psu|NC_LIV_135580 Decoy sequence, "
                             + "was | organism=Neospora_caninum | product=hypothetical "
                             + "protein | location=Neo_chrXI:5803226-5810048(+) | length=1893");
                break;
            }
        }

        List<Peptide> pepList = seqCol.getPeptide();
        assertEquals(pepList.size(), 69);
        for (Peptide pep : pepList) {
            if (pep.getId().equals("GGECTMGASSGSVDTCGETQR_1@5")) {
                assertEquals(pep.getPeptideSequence(), "GGECTMGASSGSVDTCGETQR");
                List<Modification> modList = pep.getModification();
                assertEquals(modList.size(), 3);
                Modification mod = modList.get(0);
                assertEquals(mod.getLocation(), new Integer(6));
                assertTrue(mod.getResidues().contains("M"));
                assertTrue(mod.getMonoisotopicMassDelta() - 15.994915 < diff);
                assertEquals(mod.getCvParam().get(0).getName(), "Oxidation");
                assertEquals(mod.getCvParam().get(0).getAccession(), "UNIMOD:35");
                break;
            }
        }

        List<PeptideEvidence> pepEvdList = seqCol.getPeptideEvidence();
        assertEquals(pepEvdList.size(), 71);
        for (PeptideEvidence pepEvd : pepEvdList) {
            if (pepEvd.getId().equals("PE1_2_0")) {
                assertEquals(pepEvd.getPeptideRef(), "FVGGEVLRQPPLLSR");
                assertEquals(pepEvd.getDBSequenceRef(),
                             "dbseq_Rnd2psu|NC_LIV_112900");
                assertEquals(pepEvd.getStart(), new Integer(4096));
                assertEquals(pepEvd.getEnd(), new Integer(4110));
                assertEquals(pepEvd.getPre(), "R");
                assertEquals(pepEvd.getPost(), "G");
                assertFalse(pepEvd.isIsDecoy());
            }

            if (pepEvd.getId().equals("PE39_2_70")) {
                assertEquals(pepEvd.getPeptideRef(), "LSAQRGTSSLEPPVAPR");
                assertEquals(pepEvd.getDBSequenceRef(),
                             "dbseq_Rnd3psu|NC_LIV_135580");
                assertEquals(pepEvd.getStart(), new Integer(1095));
                assertEquals(pepEvd.getEnd(), new Integer(1111));
                assertEquals(pepEvd.getPre(), "R");
                assertEquals(pepEvd.getPost(), "A");
                assertFalse(pepEvd.isIsDecoy());
            }
        }

    }

    /**
     * Test of getSpectrumIdentificationList method, of class
     * Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetSpectrumIdentificationList() {
        System.out.println("getSpectrumIdentificationList");

        SpectrumIdentificationList specIdentList
                = testContainer.getSpectrumIdentificationList();

        List<SpectrumIdentificationResult> specIdentResList = specIdentList
                .getSpectrumIdentificationResult();

        assertEquals(specIdentResList.size(), 39);
        for (SpectrumIdentificationResult sir : specIdentResList) {
            if (sir.getSpectrumID().equals("index=272")) {
                assertEquals(sir.getSpectraDataRef(), "SID_1");
                List<SpectrumIdentificationItem> specIdentItemList = sir
                        .getSpectrumIdentificationItem();
                assertEquals(specIdentItemList.size(), 6);
                for (SpectrumIdentificationItem sii : specIdentItemList) {
                    if (sii.getPeptideRef().equals("YSYGESPVNLR")) {
                        assertEquals(sii.getRank(), 3);
                        assertTrue(sii.isPassThreshold());
                        assertEquals(sii.getChargeState(), 3);
                        List<CvParam> cpList = sii.getCvParam();
                        for (CvParam cp : cpList) {
                            if (cp.getName().equals("OMSSA:evalue")) {
                                assertEquals(cp.getAccession(), "MS:1001328");
                                assertEquals(Double.parseDouble(cp.getValue()),
                                             0.671623654735357, diff);
                            }
                            if (cp.getName().equals("OMSSA:pvalue")) {
                                assertEquals(Double.parseDouble(
                                        "2.90363207020441E-6"),
                                             2.90363207020441E-6, diff);
                            }
                        }
                        assertEquals(427.871, sii.getCalculatedMassToCharge(),
                                     diff);
                        PeptideEvidence pepEvd = sii.getPeptideEvidenceRef()
                                .get(0).getPeptideEvidence();
                        assertEquals(pepEvd.getStart(), new Integer(298));
                        assertEquals(pepEvd.getEnd(), new Integer(308));
                        assertEquals(pepEvd.getPre(), "K");
                        assertEquals(pepEvd.getPost(), "S");
                        assertEquals(pepEvd.getPeptideRef(), "YSYGESPVNLR");

                        break;
                    }
                }
                break;
            }
        }

    }

}
