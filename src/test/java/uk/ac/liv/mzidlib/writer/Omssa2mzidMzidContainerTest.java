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
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSampleCollection;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftware;
import uk.ac.ebi.jmzidml.model.mzidml.AnalysisSoftwareList;
import uk.ac.ebi.jmzidml.model.mzidml.AuditCollection;
import uk.ac.ebi.jmzidml.model.mzidml.BibliographicReference;
import uk.ac.ebi.jmzidml.model.mzidml.CvList;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;
import uk.ac.ebi.jmzidml.model.mzidml.Enzyme;
import uk.ac.ebi.jmzidml.model.mzidml.Inputs;
import uk.ac.ebi.jmzidml.model.mzidml.ProteinDetectionList;
import uk.ac.ebi.jmzidml.model.mzidml.Provider;
import uk.ac.ebi.jmzidml.model.mzidml.SearchDatabaseRef;
import uk.ac.ebi.jmzidml.model.mzidml.SearchModification;
import uk.ac.ebi.jmzidml.model.mzidml.SequenceCollection;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentification;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationList;
import uk.ac.ebi.jmzidml.model.mzidml.SpectrumIdentificationProtocol;
import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.constants.CvConstants;

/**
 *
 * @author Da Qi
 */
public class Omssa2mzidMzidContainerTest {

    private static Omssa2mzidMzidContainer testContainer;
    private final double diff = 0.000001;

    public Omssa2mzidMzidContainerTest() {

    }

    @BeforeClass
    public static void setUpClass() {
        File testFile = FileUtils.getFile("src", "test", "data",
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
    @Test
    public void testGetAnalysisSampleCollection() {
        System.out.println("getAnalysisSampleCollection");
        Omssa2mzidMzidContainer instance = null;
        AnalysisSampleCollection expResult = null;
        AnalysisSampleCollection result = instance.getAnalysisSampleCollection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAnalysisSoftwareList method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetAnalysisSoftwareList() {
        System.out.println("getAnalysisSoftwareList");

        AnalysisSoftwareList softwareList = testContainer
                .getAnalysisSoftwareList();
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
    @Test
    public void testGetBibliographicReference() {
        System.out.println("getBibliographicReference");
        BibliographicReference bib = testContainer.getBibliographicReference();
        assertNull(bib);
    }

    /**
     * Test of getCvList method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetCvList() {
        System.out.println("getCvList");
        Omssa2mzidMzidContainer instance = null;
        CvList expResult = null;
        CvList result = instance.getCvList();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getInputs method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetInputs() {
        System.out.println("getInputs");
        Omssa2mzidMzidContainer instance = null;
        Inputs expResult = null;
        Inputs result = instance.getInputs();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getProteinDetectionList method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetProteinDetectionList() {
        System.out.println("getProteinDetectionList");
        Omssa2mzidMzidContainer instance = null;
        ProteinDetectionList expResult = null;
        ProteinDetectionList result = instance.getProteinDetectionList();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getProvider method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetProvider() {
        System.out.println("getProvider");
        Omssa2mzidMzidContainer instance = null;
        Provider expResult = null;
        Provider result = instance.getProvider();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSequenceCollection method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetSequenceCollection() {
        System.out.println("getSequenceCollection");
        Omssa2mzidMzidContainer instance = null;
        SequenceCollection expResult = null;
        SequenceCollection result = instance.getSequenceCollection();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getSpectrumIdentificationList method, of class
     * Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetSpectrumIdentificationList() {
        System.out.println("getSpectrumIdentificationList");
        Omssa2mzidMzidContainer instance = null;
        SpectrumIdentificationList expResult = null;
        SpectrumIdentificationList result
                = instance.getSpectrumIdentificationList();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getMatchedIon method, of class Omssa2mzidMzidContainer.
     */
    @Test
    public void testGetMatchedIon() {
        System.out.println("getMatchedIon");
        MSSpectrum spectrum = null;
        double ionMz = 0.0;
        Omssa2mzidMzidContainer instance = null;
        double[] expResult = null;
        double[] result = instance.getMatchedIon(spectrum, ionMz);
        //assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}
