package uk.ac.liv.mzidlib.util;

/**
 * This class contains some utility methods related to the Controlled Vocabulary
 * (CV) at
 *
 * http://psidev.cvs.sourceforge.net/viewvc/\*checkout*\/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
 *
 * @author Pieter
 *
 */
//CORRECT URL : see http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo
public class CVUtils {

    /**
     * Tries to infer the correct CV term id and name based on the extension of
     * the given databaseFileName. It falls back to ["MS:1001348","FASTA
     * format"] if it cannot infer the format based on the extension.
     *
     * @param databaseFileName : the file reported as the protein database where
     * the ms/ms search was done.
     *
     * @return : returns an array with 2 strings, the first one being the CV
     * term id and the second one the name.
     */
    public static String[] getDatabaseFileFormat(String databaseFileName) {
        //The default:
        //The sequence database was stored in the FASTA format
        String databaseFileFormatID = "MS:1001348";
        String databaseFileFormatName = "FASTA format";

        if (databaseFileName.toLowerCase().contains(".asn1")) {
            //The sequence database was stored in the Abstract Syntax Notation 1 format.
            databaseFileFormatID = "MS:1001349";
            databaseFileFormatName = "ASN.1";
        } else if (databaseFileName.toLowerCase().contains("formatdb")) {//TODO: not sure what the extension would be here...
            //The sequence database was stored in the NCBI formatdb (*.p*) format
            databaseFileFormatID = "MS:1001350";
            databaseFileFormatName = "NCBI *.p*";
        } else if (databaseFileName.toLowerCase().contains(".aln")) {
            //ClustalW ALN (multiple alignment) format
            databaseFileFormatID = "MS:1001351";
            databaseFileFormatName = "clustal aln";
        } else if (databaseFileName.toLowerCase().contains(".embl")) {
            //ClustalW ALN (multiple alignment) format
            databaseFileFormatID = "MS:1001352";
            databaseFileFormatName = "embl em";
        } else if (databaseFileName.toLowerCase().contains(".pir")) {
            //The NBRF PIR was used as format.
            databaseFileFormatID = "MS:1001353";
            databaseFileFormatName = "NBRF PIR";
        }
        String[] result = {databaseFileFormatID, databaseFileFormatName};
        return result;
    }

    /**
     * Tries to infer the correct CV term id and name based on the extension of
     * the given spectrumFileName. It falls back to the generic
     * ["MS:1000560","mass spectrometer file format"] if it cannot infer the
     * format based on the extension.
     *
     * @param spectrumFileName : the file reported as the ms spectra file that
     * contained the ms/ms spectra.
     *
     * @return : returns an array with 2 strings, the first one being the CV
     * term id and the second one the name.
     */
    public static String[] getMassSpecFileFormatID(String spectrumFileName) {
        //the default:
        //TODO check if: MS:1000560  --> would this be Ok as standard? Are there rules about what can be in the tags? Or should 
        //we leave them empty when we don't have the information?
        String massSpecFileFormatID = "MS:1000560";
        String massSpecFileFormatName = "mass spectrometer file format";

        if (spectrumFileName.toLowerCase().contains(".mzml")) {
            //"Proteomics Standards Inititative mzML file format."
            massSpecFileFormatID = "MS:1000584";
            massSpecFileFormatName = "mzML file";
        } else if (spectrumFileName.toLowerCase().contains(".dta")) {//TODO not sure about the extension
            //Sequest DTA file format
            massSpecFileFormatID = "MS:1000613";
            massSpecFileFormatName = "DTA file";
        } else if (spectrumFileName.toLowerCase().contains(".mzxml")) {
            //Institute of Systems Biology mzXML file format.
            massSpecFileFormatID = "MS:1000566";
            massSpecFileFormatName = "ISB mzXML file";
        } else if (spectrumFileName.toLowerCase().contains(".mzdata")) {
            //Proteomics Standards Inititative mzData file format
            massSpecFileFormatID = "MS:1000564";
            massSpecFileFormatName = "PSI mzData file";
        } else if (spectrumFileName.toLowerCase().contains(".mgf")) {
            //Mascot MGF file
            massSpecFileFormatID = "MS:1001062";
            massSpecFileFormatName = "Mascot MGF file";
        } else if (spectrumFileName.toLowerCase().contains(".mz5")) {
            //mz5 file format, modeled after mzML
            massSpecFileFormatID = "MS:1001881";
            massSpecFileFormatName = "mz5 file";
        }

        //TODO - other important formats (but we have to figure out what are their extensions first: 
        /*    		
         id: MS:1001463
         name: Phenyx XML format
         is_a: MS:1000560 ! mass spectrometer file format
         is_a: MS:1001040 ! intermediate analysis format
		
         id: MS:1001509
         name: Agilent MassHunter file
         def: "A data file found in an Agilent MassHunter directory which contains raw data acquired by an Agilent mass spectrometer." [PSI:PI]
         is_a: MS:1000560 ! mass spectrometer file format
		
         id: MS:1001527
         name: Proteinscape spectra
         def: "Spectra from Bruker/Protagen Proteinscape database." [PSI:MS]
         is_a: MS:1000560 ! mass spectrometer file format
		

         [Term]
         id: MS:1000563
         name: Thermo RAW file
         def: "Thermo Scientific RAW file format." [PSI:MS]
         is_a: MS:1000560 ! mass spectrometer file format

         [Term]
         id: MS:1000565
         name: Micromass PKL file
         def: "Micromass PKL file format." [PSI:MS]
         is_a: MS:1000560 ! mass spectrometer file format
         */
        String[] result = {massSpecFileFormatID, massSpecFileFormatName};
        return result;
    }
}
