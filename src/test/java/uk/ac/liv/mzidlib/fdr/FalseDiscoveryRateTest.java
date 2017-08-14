/*
 * Date: 26-Jan-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.fdr.FalseDiscoveryRateTest.java
 *
 * jmzquantml is Copyright 2017 University of Liverpool.
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

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.model.mzidml.AnalysisCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisProtocolCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.DBSequence;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.Peptide;
import uk.ac.ebi.jmzidml.model.mzidml.PeptideEvidence;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationItem;

/**
 *
 * @author Da Qi
 */
public class FalseDiscoveryRateTest {

    static FalseDiscoveryRate fdr;
    
    public FalseDiscoveryRateTest() {
        
    }

    @BeforeClass
    public static void setUpClass() {
        String fileName = "Adult_Adrenalgland_Gel_Velos_2_f26.t_tandem.mzid";
        String searchEngine = "Search Engine";
        String decoyRatio = "10";
        String decoy = "REVERSE";
        String mzidFile = FalseDiscoveryRateTest.class.getClassLoader().
                getResource(
                        fileName).getFile();
        fdr = new FalseDiscoveryRate(mzidFile, searchEngine,
                                                        decoyRatio, decoy);
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of main method, of class FalseDiscoveryRate.
     */
    @Test
    public void testMain() {
        
        System.out.println("main");
        String[] args = null;
        FalseDiscoveryRate.main(args);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFromXMLPeptideModificationHash method, of class
     * FalseDiscoveryRate.
     */
    @Test
    public void testGetFromXMLPeptideModificationHash() {
        System.out.println("getFromXMLPeptideModificationHash");
        FalseDiscoveryRate instance = null;
        Map<String, List<List<String>>> expResult = null;
        Map<String, List<List<String>>> result
                = instance.getFromXMLPeptideModificationHash();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFromXMLPeptideSequenceHash method, of class
     * FalseDiscoveryRate.
     */
    @Test
    public void testGetFromXMLPeptideSequenceHash() {
        System.out.println("getFromXMLPeptideSequenceHash");
        FalseDiscoveryRate instance = null;
        Map<String, String> expResult = null;
        Map<String, String> result = instance.getFromXMLPeptideSequenceHash();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFromXMLSpectrumInfoHash method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetFromXMLSpectrumInfoHash() {
        System.out.println("getFromXMLSpectrumInfoHash");
        FalseDiscoveryRate instance = null;
        Map<String, List<List<Object>>> expResult = null;
        Map<String, List<List<Object>>> result
                = instance.getFromXMLSpectrumInfoHash();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeTheSortedDataToFile method, of class FalseDiscoveryRate.
     */
    @Test
    public void testWriteTheSortedDataToFile()
            throws Exception {
        System.out.println("writeTheSortedDataToFile");
        String fileName = "";
        FalseDiscoveryRate instance = null;
        instance.writeTheSortedDataToFile(fileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeMzidFile method, of class FalseDiscoveryRate.
     */
    @Test
    public void testWriteMzidFile() {
        System.out.println("writeMzidFile");
        String csvFileName = "";
        FalseDiscoveryRate instance = null;
        instance.writeMzidFile(csvFileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToMzIdentMLFile method, of class FalseDiscoveryRate.
     */
    @Test
    public void testWriteToMzIdentMLFile() {
        System.out.println("writeToMzIdentMLFile");
        String fileName = "";
        FalseDiscoveryRate instance = null;
        instance.writeToMzIdentMLFile(fileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToCsvFile method, of class FalseDiscoveryRate.
     */
    @Test
    public void testWriteToCsvFile()
            throws Exception {
        System.out.println("writeToCsvFile");
        String fileName = "";
        FalseDiscoveryRate instance = null;
        instance.writeToCsvFile(fileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToTsvFile method, of class FalseDiscoveryRate.
     */
    @Test
    public void testWriteToTsvFile()
            throws Exception {
        System.out.println("writeToTsvFile");
        String fileName = "";
        FalseDiscoveryRate instance = null;
        instance.writeToTsvFile(fileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of computeFDRusingJonesMethod method, of class FalseDiscoveryRate.
     */
    @Test
    public void testComputeFDRusingJonesMethod() {
        System.out.println("computeFDRusingJonesMethod");
        FalseDiscoveryRate instance = null;
        instance.computeFDRusingJonesMethod();
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_spectrumResult method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_spectrumResult() {
        System.out.println("getSorted_spectrumResult");
        FalseDiscoveryRate instance = null;
        List<String> expResult = null;
        List<String> result = instance.getSorted_spectrumResult();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_peptideNames method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_peptideNames() {
        System.out.println("getSorted_peptideNames");
        FalseDiscoveryRate instance = null;
        List<String> expResult = null;
        List<String> result = instance.getSorted_peptideNames();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_evalues method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_evalues() {
        System.out.println("getSorted_evalues");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getSorted_evalues();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_scores method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_scores() {
        System.out.println("getSorted_scores");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getSorted_scores();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_decoyOrNot method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_decoyOrNot() {
        System.out.println("getSorted_decoyOrNot");
        FalseDiscoveryRate instance = null;
        List<String> expResult = null;
        List<String> result = instance.getSorted_decoyOrNot();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_simpleFDR method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_simpleFDR() {
        System.out.println("getSorted_simpleFDR");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getSorted_simpleFDR();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_qValues method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_qValues() {
        System.out.println("getSorted_qValues");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getSorted_qValues();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSorted_estimatedFDR method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSorted_estimatedFDR() {
        System.out.println("getSorted_estimatedFDR");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getSorted_estimatedFDR();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getTP method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetTP() {
        System.out.println("getTP");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getTP();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFP method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetFP() {
        System.out.println("getFP");
        FalseDiscoveryRate instance = null;
        List<Double> expResult = null;
        List<Double> result = instance.getFP();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of clearAllData method, of class FalseDiscoveryRate.
     */
    @Test
    public void testClearAllData() {
        System.out.println("clearAllData");
        FalseDiscoveryRate instance = null;
        instance.clearAllData();
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getdBSequenceHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetdBSequenceHashMap() {
        System.out.println("getdBSequenceHashMap");
        FalseDiscoveryRate instance = null;
        Map<String, DBSequence> expResult = null;
        Map<String, DBSequence> result = instance.getdBSequenceHashMap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setdBSequenceHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testSetdBSequenceHashMap() {
        System.out.println("setdBSequenceHashMap");
        Map<String, DBSequence> dBSequenceHashMap = null;
        FalseDiscoveryRate instance = null;
        instance.setdBSequenceHashMap(dBSequenceHashMap);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getPeptideEvidenceHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetPeptideEvidenceHashMap() {
        System.out.println("getPeptideEvidenceHashMap");
        FalseDiscoveryRate instance = null;
        Map<String, PeptideEvidence> expResult = null;
        Map<String, PeptideEvidence> result
                = instance.getPeptideEvidenceHashMap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getPeptideEvidenceMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetPeptideEvidenceMap() {
        System.out.println("getPeptideEvidenceMap");
        FalseDiscoveryRate instance = null;
        Map<String, PeptideEvidence> expResult = null;
        Map<String, PeptideEvidence> result = instance.getPeptideEvidenceMap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setPeptideEvidenceHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testSetPeptideEvidenceHashMap() {
        System.out.println("setPeptideEvidenceHashMap");
        Map<String, PeptideEvidence> peptideEvidenceHashMap = null;
        FalseDiscoveryRate instance = null;
        instance.setPeptideEvidenceHashMap(peptideEvidenceHashMap);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getPeptideHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetPeptideHashMap() {
        System.out.println("getPeptideHashMap");
        FalseDiscoveryRate instance = null;
        Map<String, Peptide> expResult = null;
        Map<String, Peptide> result = instance.getPeptideHashMap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setPeptideHashMap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testSetPeptideHashMap() {
        System.out.println("setPeptideHashMap");
        Map<String, Peptide> peptideHashMap = null;
        FalseDiscoveryRate instance = null;
        instance.setPeptideHashMap(peptideHashMap);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSpectrumIdentificationItemHashMap method, of class
     * FalseDiscoveryRate.
     */
    @Test
    public void testGetSpectrumIdentificationItemHashMap() {
        System.out.println("getSpectrumIdentificationItemHashMap");
        FalseDiscoveryRate instance = null;
        Map<String, SpectrumIdentificationItem> expResult = null;
        Map<String, SpectrumIdentificationItem> result
                = instance.getSpectrumIdentificationItemHashMap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setSpectrumIdentificationItemHashMap method, of class
     * FalseDiscoveryRate.
     */
    @Test
    public void testSetSpectrumIdentificationItemHashMap() {
        System.out.println("setSpectrumIdentificationItemHashMap");
        Map<String, SpectrumIdentificationItem> spectrumIdentificationItemHashMap
                = null;
        FalseDiscoveryRate instance = null;
        instance.setSpectrumIdentificationItemHashMap(
                spectrumIdentificationItemHashMap);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAnalysisSoftwareList method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetAnalysisSoftwareList() {
        System.out.println("getAnalysisSoftwareList");
        FalseDiscoveryRate instance = null;
        AnalysisSoftwareList expResult = null;
        AnalysisSoftwareList result = instance.getAnalysisSoftwareList();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAuditCollection method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetAuditCollection() {
        System.out.println("getAuditCollection");
        FalseDiscoveryRate instance = null;
        AuditCollection expResult = null;
        AuditCollection result = instance.getAuditCollection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getProvider method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetProvider() {
        System.out.println("getProvider");
        FalseDiscoveryRate instance = null;
        Provider expResult = null;
        Provider result = instance.getProvider();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAnalysisProtocolCollection method, of class
     * FalseDiscoveryRate.
     */
    @Test
    public void testGetAnalysisProtocolCollection() {
        System.out.println("getAnalysisProtocolCollection");
        FalseDiscoveryRate instance = null;
        AnalysisProtocolCollection expResult = null;
        AnalysisProtocolCollection result
                = instance.getAnalysisProtocolCollection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getCvList method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetCvList() {
        System.out.println("getCvList");
        FalseDiscoveryRate instance = null;
        CvList expResult = null;
        CvList result = instance.getCvList();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAnalysisCollection method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetAnalysisCollection() {
        System.out.println("getAnalysisCollection");
        FalseDiscoveryRate instance = null;
        AnalysisCollection expResult = null;
        AnalysisCollection result = instance.getAnalysisCollection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getInputs method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetInputs() {
        System.out.println("getInputs");
        FalseDiscoveryRate instance = null;
        Inputs expResult = null;
        Inputs result = instance.getInputs();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSearchDatabase_Ref method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSearchDatabase_Ref() {
        System.out.println("getSearchDatabase_Ref");
        FalseDiscoveryRate instance = null;
        String expResult = "";
        String result = instance.getSearchDatabase_Ref();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getUnimodHashmap method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetUnimodHashmap() {
        System.out.println("getUnimodHashmap");
        FalseDiscoveryRate instance = null;
        Map<String, String> expResult = null;
        Map<String, String> result = instance.getUnimodHashmap();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSpectraData_ref method, of class FalseDiscoveryRate.
     */
    @Test
    public void testGetSpectraData_ref() {
        System.out.println("getSpectraData_ref");
        FalseDiscoveryRate instance = null;
        String expResult = "";
        String result = instance.getSpectraData_ref();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}
