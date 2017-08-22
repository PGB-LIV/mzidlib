/*
 * Date: 11-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.fdr.FalseDiscoveryRateGlobalTest.java
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

package uk.ac.liv.mzidlib.fdr;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import javax.xml.bind.JAXBException;

import org.apache.commons.io.FileUtils;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisData;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationResult;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.liv.mzidlib.constants.CvConstants;

/**
 *
 * @author Da Qi
 */
public class FalseDiscoveryRateGlobalTest {

    private final FalseDiscoveryRateGlobal fdrPSMSeq;
    private final FalseDiscoveryRateGlobal fdrPSMPar;
    private final FalseDiscoveryRateGlobal fdrPep;
    private final double DIFF = 0.000001;

    public FalseDiscoveryRateGlobalTest() {

        File testFile = FileUtils.getFile("src","test","data","Adult_Adrenalgland_Gel_Velos_2_f26.t_tandem.mzid");
        
        
                
        fdrPSMSeq
                = new FalseDiscoveryRateGlobal(testFile.getAbsolutePath(),
                                               "0.01", "REVERSED",
                                               CvConstants.XTANDEM_EXPECT.
                                               getAccession(), true, "PSM",
        
                        "PAG", "1.2");
        final long startFdrPSMSeq = System.currentTimeMillis();
        
        fdrPSMSeq.computeFDRusingJonesMethod();
        
        final long endFdrPSMSeq = System.currentTimeMillis();
        
        fdrPSMPar
                = new FalseDiscoveryRateGlobal(testFile.getAbsolutePath(),
                                               "0.01", "REVERSED",
                                               CvConstants.XTANDEM_EXPECT.
                                               getAccession(), true, "PSM",
        
                        "PAG", "1.2");
        
        final long startFdrPSMPar = System.currentTimeMillis();
        
        fdrPSMPar.computeFDRusingJonesMethodPar();
        
        final long endFdrPSMPar = System.currentTimeMillis();
        
        System.out.println("FalseDiscoveryRateGlobal.computeFDRusingJonesMethod() takes " + (endFdrPSMSeq - startFdrPSMSeq) + " Millis.\n");
        
        System.out.println("FalseDiscoveryRateGlobal.computeFDRusingJonesMethodPar() takes " + (endFdrPSMPar - startFdrPSMPar) + " Millis.\n");
        
        System.out.println("Speedup is: " + (double) (endFdrPSMSeq - startFdrPSMSeq) / (double) (endFdrPSMPar - startFdrPSMPar) + ".\n");
        
        
        fdrPep
                = new FalseDiscoveryRateGlobal(testFile.getAbsolutePath(),
                                               "0.01", "REVERSED",
                                               CvConstants.XTANDEM_EXPECT.
                                               getAccession(), true, "Peptide",
                                               "PAG", "1.2");
        
        fdrPep.computeFDRusingJonesMethod();
    }

    @BeforeClass
    public static void setUpClass() {
        List<Double> evalues = Arrays.asList(new Double[]{0.2, 0.002, 0.023,
            0.43, 0.0056, 0.0021, 0.0082, 0.00045, 0.0027, 0.065, 0.000032,
            0.045});

    }

    @AfterClass
    public static void tearDownClass() {
    }

//    /**
//     * Test of computeFDRusingJonesMethod method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testComputeFDRusingJonesMethod() {
//        System.out.println("computeFDRusingJonesMethod");
//        FalseDiscoveryRateGlobal instance = null;
//        instance.computeFDRusingJonesMethod();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    /**
     * Test of writeToMzIdentMLFile method, of class FalseDiscoveryRateGlobal.
     *
     * @throws java.io.IOException
     * @throws javax.xml.bind.JAXBException
     */
    @Test
    public void testWriteToMzIdentMLFile()
            throws IOException, JAXBException {
        System.out.println("writeToMzIdentMLFile");
        System.out.println("Test peptide level");
        File tempPepFile = File.createTempFile("fdr-peptide-test", ".mzid");
        tempPepFile.deleteOnExit();
        fdrPep.writeToMzIdentMLFile(tempPepFile.getAbsolutePath());

        MzIdentMLUnmarshaller pepUm = new MzIdentMLUnmarshaller(tempPepFile);

        AnalysisData ad = pepUm.unmarshal(MzIdentMLElement.AnalysisData);
        List<SpectrumIdentificationResult> sirList = ad.
                getSpectrumIdentificationList().get(0).
                getSpectrumIdentificationResult();
        for (SpectrumIdentificationResult sir : sirList) {
            if (sir.getId().equals("SIR_650")) {

                // test SpectrumIdentificationResult attributes
                assertEquals(sir.getSpectrumID(), "index=723");
                assertEquals(sir.getSpectraDataRef(), "SID_1");
                assertEquals(sir.getSpectrumIdentificationItem().size(), 1);

                // test SpectrumIdentificationItem attributes
                SpectrumIdentificationItem sii = sir.
                        getSpectrumIdentificationItem().
                        get(0);
                assertEquals(436.897851, sii.getCalculatedMassToCharge(), DIFF);
                assertEquals(437.237152, sii.getExperimentalMassToCharge(), DIFF);
                assertEquals(sii.getPeptideRef(), "GLGSIFGSVGGETK___");
                assertEquals(sii.getRank(), 1);
                assertEquals(sii.getChargeState(), 3);
                assertTrue(sii.isPassThreshold());
                assertEquals(sii.getId(), "SII_650_1");

                // test PeptideEvidenceRef
                assertEquals(sii.getPeptideEvidenceRef().size(), 6);
                assertEquals(sii.getPeptideEvidenceRef().get(3).
                        getPeptideEvidenceRef(),
                             "PE650_2_4098");

                int num = 0;
                for (CvParam cp : sii.getCvParam()) {

                    if (cp.getAccession().equals(
                            CvConstants.DISTINCT_PEPTIDE_LEVEL_LOCAL_FDR.
                            getAccession())) {
                        num++;
                        assertEquals(53.92986698911729, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    } else if (cp.getAccession().equals(
                            CvConstants.DISTINCT_PEPTIDE_LEVEL_Q_VALUE.
                            getAccession())) {
                        num++;
                        assertEquals(53.89728096676737, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    } else if (cp.getAccession().equals(
                            CvConstants.DISTINCT_PEPTIDE_LEVEL_FDRSCORE.
                            getAccession())) {
                        num++;
                        assertEquals(53.89728096676737, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    }
                }
                assertEquals("Doesn't contain all peptide level score:", num, 3);
            }
        }
        
        System.out.println("Test PSM level");
        File tempPSMFile = File.createTempFile("fdr-psm-test", ".mzid");
        tempPSMFile.deleteOnExit();
        fdrPSMSeq.writeToMzIdentMLFile(tempPSMFile.getAbsolutePath());
        MzIdentMLUnmarshaller psmUm = new MzIdentMLUnmarshaller(tempPSMFile);

        ad = psmUm.unmarshal(MzIdentMLElement.AnalysisData);
        sirList = ad.
                getSpectrumIdentificationList().get(0).
                getSpectrumIdentificationResult();
        for (SpectrumIdentificationResult sir : sirList) {
            if (sir.getId().equals("SIR_1")) {

                // test SpectrumIdentificationResult attributes
                assertEquals(sir.getSpectrumID(), "index=14");
                assertEquals(sir.getSpectraDataRef(), "SID_1");
                assertEquals(sir.getSpectrumIdentificationItem().size(), 1);

                // test SpectrumIdentificationItem attributes
                SpectrumIdentificationItem sii = sir.
                        getSpectrumIdentificationItem().
                        get(0);
                assertEquals(958.6044, sii.getCalculatedMassToCharge(), DIFF);
                assertEquals(958.597168, sii.getExperimentalMassToCharge(), DIFF);
                assertEquals(sii.getPeptideRef(), "KVSVLKER___");
                assertEquals(sii.getRank(), 1);
                assertEquals(sii.getChargeState(), 1);
                assertTrue(sii.isPassThreshold());
                assertEquals(sii.getId(), "SII_1_1");

                assertEquals(sir.getCvParam().size(), 1);
                CvParam cvp = sir.getCvParam().get(0);
                assertEquals(cvp.getValue(), "Adult_Adrenalgland_Gel_Velos_2_f26.499.499.1 RTINSECONDS=835.9027");
                
                // test PeptideEvidenceRef
                assertEquals(sii.getPeptideEvidenceRef().size(), 2);
                assertEquals(sii.getPeptideEvidenceRef().get(0).
                        getPeptideEvidenceRef(),
                             "PE1_2_0");

                int num = 0;
                for (CvParam cp : sii.getCvParam()) {

                    if (cp.getAccession().equals(
                            CvConstants.PSM_LEVEL_LOCAL_FDR.
                            getAccession())) {
                        num++;
                        assertEquals(53.943044906900326, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    } else if (cp.getAccession().equals(
                            CvConstants.PSM_LEVEL_Q_VALUE.
                            getAccession())) {
                        num++;
                        assertEquals(53.92827812756639, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    } else if (cp.getAccession().equals(
                            CvConstants.PSM_LEVEL_FDRSCORE.
                            getAccession())) {
                        num++;
                        assertEquals(53.92827812756639, Double.
                                     parseDouble(cp.getValue()), DIFF);
                    }
                }
                assertEquals("Doesn't contain all PSM level score:", num, 3);
            }
        }
    }

    /**
     * Test of getSorted_spectrumResult method, of class
     * FalseDiscoveryRateGlobal.
     */
//    @Test
//    public void testGetSorted_spectrumResult() {
//        System.out.println("getSorted_spectrumResult");
//        FalseDiscoveryRateGlobal instance = null;
//        List<String> expResult = null;
//        List<String> result = instance.getSorted_spectrumResult();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_peptideNames method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_peptideNames() {
//        System.out.println("getSorted_peptideNames");
//        FalseDiscoveryRateGlobal instance = null;
//        List<String> expResult = null;
//        List<String> result = instance.getSorted_peptideNames();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_evalues method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_evalues() {
//        System.out.println("getSorted_evalues");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getSorted_evalues();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_decoyOrNot method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_decoyOrNot() {
//        System.out.println("getSorted_decoyOrNot");
//        FalseDiscoveryRateGlobal instance = null;
//        List<String> expResult = null;
//        List<String> result = instance.getSorted_decoyOrNot();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_simpleFDR method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_simpleFDR() {
//        System.out.println("getSorted_simpleFDR");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getSorted_simpleFDR();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_qValues method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_qValues() {
//        System.out.println("getSorted_qValues");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getSorted_qValues();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSorted_estimatedFDR method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSorted_estimatedFDR() {
//        System.out.println("getSorted_estimatedFDR");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getSorted_estimatedFDR();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getTP method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetTP() {
//        System.out.println("getTP");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getTP();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getFP method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetFP() {
//        System.out.println("getFP");
//        FalseDiscoveryRateGlobal instance = null;
//        List<Double> expResult = null;
//        List<Double> result = instance.getFP();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of clearAllData method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testClearAllData() {
//        System.out.println("clearAllData");
//        FalseDiscoveryRateGlobal instance = null;
//        instance.clearAllData();
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getdBSequenceHashMap method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetdBSequenceHashMap() {
//        System.out.println("getdBSequenceHashMap");
//        FalseDiscoveryRateGlobal instance = null;
//        Map<String, DBSequence> expResult = null;
//        Map<String, DBSequence> result = instance.getdBSequenceHashMap();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setdBSequenceHashMap method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testSetdBSequenceHashMap() {
//        System.out.println("setdBSequenceHashMap");
//        Map<String, DBSequence> dBSequenceMap = null;
//        FalseDiscoveryRateGlobal instance = null;
//        instance.setdBSequenceHashMap(dBSequenceMap);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPeptideEvidenceHashMap method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetPeptideEvidenceHashMap() {
//        System.out.println("getPeptideEvidenceHashMap");
//        FalseDiscoveryRateGlobal instance = null;
//        Map<String, PeptideEvidence> expResult = null;
//        Map<String, PeptideEvidence> result
//                = instance.getPeptideEvidenceHashMap();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setPeptideEvidenceHashMap method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testSetPeptideEvidenceHashMap() {
//        System.out.println("setPeptideEvidenceHashMap");
//        Map<String, PeptideEvidence> peptideEvidenceMap = null;
//        FalseDiscoveryRateGlobal instance = null;
//        instance.setPeptideEvidenceHashMap(peptideEvidenceMap);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPeptideHashMap method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetPeptideHashMap() {
//        System.out.println("getPeptideHashMap");
//        FalseDiscoveryRateGlobal instance = null;
//        Map<String, Peptide> expResult = null;
//        Map<String, Peptide> result = instance.getPeptideHashMap();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setPeptideHashMap method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testSetPeptideHashMap() {
//        System.out.println("setPeptideHashMap");
//        Map<String, Peptide> peptideMap = null;
//        FalseDiscoveryRateGlobal instance = null;
//        instance.setPeptideHashMap(peptideMap);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAnalysisSoftwareList method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetAnalysisSoftwareList() {
//        System.out.println("getAnalysisSoftwareList");
//        FalseDiscoveryRateGlobal instance = null;
//        AnalysisSoftwareList expResult = null;
//        AnalysisSoftwareList result = instance.getAnalysisSoftwareList();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAuditCollection method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetAuditCollection() {
//        System.out.println("getAuditCollection");
//        FalseDiscoveryRateGlobal instance = null;
//        AuditCollection expResult = null;
//        AuditCollection result = instance.getAuditCollection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getProvider method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetProvider() {
//        System.out.println("getProvider");
//        FalseDiscoveryRateGlobal instance = null;
//        Provider expResult = null;
//        Provider result = instance.getProvider();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAnalysisProtocolCollection method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetAnalysisProtocolCollection() {
//        System.out.println("getAnalysisProtocolCollection");
//        FalseDiscoveryRateGlobal instance = null;
//        AnalysisProtocolCollection expResult = null;
//        AnalysisProtocolCollection result
//                = instance.getAnalysisProtocolCollection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getCvList method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetCvList() {
//        System.out.println("getCvList");
//        FalseDiscoveryRateGlobal instance = null;
//        CvList expResult = null;
//        CvList result = instance.getCvList();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getAnalysisCollection method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetAnalysisCollection() {
//        System.out.println("getAnalysisCollection");
//        FalseDiscoveryRateGlobal instance = null;
//        AnalysisCollection expResult = null;
//        AnalysisCollection result = instance.getAnalysisCollection();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getInputs method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetInputs() {
//        System.out.println("getInputs");
//        FalseDiscoveryRateGlobal instance = null;
//        Inputs expResult = null;
//        Inputs result = instance.getInputs();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getFromXMLPeptideSequenceHash method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetFromXMLPeptideSequenceHash() {
//        System.out.println("getFromXMLPeptideSequenceHash");
//        FalseDiscoveryRateGlobal instance = null;
//        Map<String, String> expResult = null;
//        Map<String, String> result = instance.getFromXMLPeptideSequenceHash();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSearchDatabase_Ref method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testGetSearchDatabase_Ref() {
//        System.out.println("getSearchDatabase_Ref");
//        FalseDiscoveryRateGlobal instance = null;
//        String expResult = "";
//        String result = instance.getSearchDatabase_Ref();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of writeToTsvFile method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testWriteToTsvFile()
//            throws Exception {
//        System.out.println("writeToTsvFile");
//        String fileName = "";
//        FalseDiscoveryRateGlobal instance = null;
//        instance.writeToTsvFile(fileName);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of writeToCsvFile method, of class FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testWriteToCsvFile()
//            throws Exception {
//        System.out.println("writeToCsvFile");
//        String fileName = "";
//        FalseDiscoveryRateGlobal instance = null;
//        instance.writeToCsvFile(fileName);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of writeTheSortedDataToFile method, of class
//     * FalseDiscoveryRateGlobal.
//     */
//    @Test
//    public void testWriteTheSortedDataToFile()
//            throws Exception {
//        System.out.println("writeTheSortedDataToFile");
//        String fileName = "";
//        FalseDiscoveryRateGlobal instance = null;
//        instance.writeTheSortedDataToFile(fileName);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
}
