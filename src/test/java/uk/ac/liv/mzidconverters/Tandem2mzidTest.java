package uk.ac.liv.mzidconverters;

import junit.framework.TestCase;
import uk.ac.liv.mzidlib.MzIdentMLLib;
import uk.ac.liv.mzidlib.converters.Tandem2mzid;



/**
 * This class can be used to test the results of Tandem2mzid transformation tool.
 * 
 *      
 * @author plukasse
 *
 */
public class Tandem2mzidTest extends TestCase 
{
	
	public void test_basic2() throws Exception
	{
		String xTandemFile = "src/test/data/OutputBrueba.2014_09_23_16_44_27.t.xml";
		String resultFile = "src/test/out/OutputBrueba.2014_09_23_16_44_27.t.xml.mzid";
		new Tandem2mzid(xTandemFile, resultFile);
		
		
	}
	
	/**
	 * Testing the basic functionality
	 * @throws Exception 
	 */
	public void test_basic() throws Exception
	{
		String xTandemFile = "src/test/data/55merge_tandem.xml";
		String resultFile = "src/test/data/55merge_tandem.xml.mzid";
		new Tandem2mzid(xTandemFile, resultFile, "MS:1001348","MS:1001062",false, true);
		
		
	}
	
	
	public void test_basic_commandLine() throws Exception
	{
		String xTandemFile = "src/test/data/55merge_tandem.xml";
		String resultFile = "src/test/data/55merge_tandem.xml.mzid";
		String args[] = {"Tandem2mzid", xTandemFile, resultFile, "-decoyRegex", "a regex", "-outputFragmentation", "true"};
		
		MzIdentMLLib.main(args);
		
		
	}
	

	public void test_galaxyFile() throws Exception
	{
		String xTandemFileFromGalaxyStep = "src/test/data/Galaxy-[X_Tandem_on_data_N].bioml";
		new Tandem2mzid(xTandemFileFromGalaxyStep, xTandemFileFromGalaxyStep + ".mzid", "MS:1001348","MS:1000584",true, true);
	}
		
	public void test_debug_File() throws Exception
	{
		String xTandemFileFromGalaxyStep = "src/test/data/DEBUG_MSMS_XTANDEM.bioml";
		String outMzid = xTandemFileFromGalaxyStep + ".mzid";
		new Tandem2mzid(xTandemFileFromGalaxyStep, outMzid, "MS:1001348","MS:1000584",false, true);
		
		//Possible check if code below would work...:
		/*
		uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller unmarshaller = new uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller(new File(outMzid));
		SpectrumIdentificationItem firstIdentification = (SpectrumIdentificationItem)unmarshaller.unmarshal(".//SpectrumIdentificationItem");
		//Get 1st identification:
		assertTrue(firstIdentification.getExperimentalMassToCharge() == 464.285614);
		*/
	}
	
	
	public void test_debug_File_with_protRegex() throws Exception
	{
		String xTandemFile = "src/test/data/DEBUG_MSMS_XTANDEM.bioml";
		String outMzid = xTandemFile + "_REGEX_TEST.mzid";
		new Tandem2mzid(
				xTandemFile, 
				outMzid, 
				"MS:1001348","MS:1000584", //databaseFileFormatID,  massSpecFileFormatID,
				false, //isMs2SpectrumIdStartingAtZero 
				null, //decoyRegularExpression
				"\\S+", //proteinCodeRegex
				true); //outputFragmentation
	}
	
	public void test_proteinCodeRegex_commandLine() throws Exception
	{
		String xTandemFile = "src/test/data/DEBUG_MSMS_XTANDEM.bioml";
		String outMzid = xTandemFile + "_REGEX_TEST_CMDLINE.mzid";
		String args[] = {"Tandem2mzid", xTandemFile, outMzid, 
				"-decoyRegex", "a regex",
				"-proteinCodeRegex", "\\S+",
				"-outputFragmentation", "true"};
		
		MzIdentMLLib.main(args);
		
		
	}
	
	public void test_multiple_sii_case() throws Exception
	{
		String xTandemFileFromGalaxyStep = "src/test/data/20091015_spiking_0_05fmol_CytC_MSMS_XTANDEM.out";
		String outMzid = xTandemFileFromGalaxyStep + ".mzid";
		new Tandem2mzid(xTandemFileFromGalaxyStep, outMzid, "MS:1001348","MS:1000584",true, true);
	}
	
	
	public void test_modifications() throws Exception
	{
		String xTandemFile = "src/test/data/S2_depl_spiked_6.mzML_xtandemOut.bioml";
		String resultFile = "src/test/data/S2_depl_spiked_6.mzML_xtandemOut.bioml.mzid";
		new Tandem2mzid(xTandemFile, resultFile, "MS:1001348","MS:1000584",true, true);
		
		
	}
	
	
	public void test_substitutions() throws Exception
	{
		String xTandemFile = "src/test/data/55merge_tandem_WITH_SUBSTITUTIONS.xml";
		String resultFile = "src/test/data/55merge_tandem_WITH_SUBSTITUTIONS.xml.mzid";
		
		//Also testing different databaseformat and spectrum file formats:
		new Tandem2mzid(xTandemFile, resultFile, "MS:1001349","MS:1000584",false, true);
		
		
	}
	
	
	/**
	 * Test for X!Tandem file generated with the option "output, parameters=no" 
	 * (see http://thegpm.org/tandem/api/opara.html )
	 * 
	 * @throws Exception
	 */
	public void test_ABSENT_PARAMETERS_OUTPUT() throws Exception
	{
		String xTandemFile = "src/test/data/DEBUG_NO_PARAMETERS_IN_OUTPUT_XTANDEM.bioml";
		String resultFile = "src/test/data/DEBUG_NO_PARAMETERS_IN_OUTPUT_XTANDEM.bioml.mzid";
		try
		{
			new Tandem2mzid(xTandemFile, resultFile, "MS:1001348","MS:1000584",false, true);
		}
		catch(Exception e)
		{
			assertTrue(e.getMessage().contains("Expected parameter not found in X!Tandem file"));
		}
	}
	
			
}
