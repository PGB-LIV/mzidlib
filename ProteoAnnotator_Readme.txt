=========================================================================================================================
=                                                                     													=
=                                                                     													=
=                ProteoAnnotator – Open Source Proteogenomics Annotation Software Supporting PSI Standards              =
=                                                                     													=
=                                           http://www.proteoannotator.org/   											=
=                                                                     													=
=========================================================================================================================


Description
	ProteoAnnotator – Open Source Proteogenomics Annotation Software Supporting PSI Standards        

How to test ProteoAnnotator

	1- Download the latest version from: http://www.proteoannotator.org/?q=installation
	2- Extract the zip file. Example, extract ProteoAnnotator-1.0.zip
	3- Download the test files from: http://www.proteoannotator.org/datasets/releases/ProteoAnnotator_test_files.zip
	4- Extract the zip file. Example, extract the ProteoAnnotator_test_files.zip
	5- Change the directory to the extracted ProteoAnnotator-1.0.zip. Example change the directory to mzidentml-lib-1.6.10-SNAPSHOT
	6- Run the following command:

			java -Xmx10G -jar mzidentml-lib-1.6.10-SNAPSHOT.jar ProteoAnnotator -inputGFF PATH_TO_ProteoAnnotator_test_files\canonical_gene_model\ToxoDB-10.0_TgondiiME49.gff -spectrum_files PATH_TO_ProteoAnnotator_test_files\mgf\Orbitap-1DE\combined.mgf -outputFolder PATH_TO_ProteoAnnotator_test_files\output_results -inputPredicted "PATH_TO_ProteoAnnotator_test_files\non_canonical_gene_model\augustus\augustus.gff;PATH_TO_ProteoAnnotator_test_files\non_canonical_gene_model\augustus\augustus.aa##PATH_TO_ProteoAnnotator_test_files\non_canonical_gene_model\glimmer\Toxo-Glimmer-ME49.gff" -prefix orbitap-1DE_toxo10_ -compress false -searchParameters "PATH_TO_ProteoAnnotator_test_files\search.txt"

IMMPORTANT:
	Replace PATH_TO_ProteoAnnotator_test_files with the actual path to the unzipped ProteoAnnotator_test_files folder.

Example (change it to match your folder structures):

			java -Xmx10G -jar mzidentml-lib-1.6.10-SNAPSHOT.jar ProteoAnnotator -inputGFF C:\ProteoAnnotator_test_files\canonical_gene_model\ToxoDB-10.0_TgondiiME49.gff -spectrum_files C:\ProteoAnnotator_test_files\mgf\Orbitap-1DE\combined.mgf -outputFolder C:\ProteoAnnotator_test_files\output_results -inputPredicted "C:\ProteoAnnotator_test_files\non_canonical_gene_model\augustus\augustus.gff;C:\ProteoAnnotator_test_files\non_canonical_gene_model\augustus\augustus.aa##C:\ProteoAnnotator_test_files\non_canonical_gene_model\glimmer\Toxo-Glimmer-ME49.gff" -prefix orbitap-1DE_toxo10_ -compress false -searchParameters "C:\ProteoAnnotator_test_files\search.txt"

			In this example the "PATH_TO_ProteoAnnotator_test_files" is "C:\ProteoAnnotator_test_files"
	
Command line parameters:
	-prefix (optional): a prefix to be attached to the output file names. 
	-inputGFF_A (mandatory): The canonical GFF file.
	-inputFasta_A (optional if the canonical GFF contains the FASTA): The protein database.
	-outputFolder (mandatory): The output folder for the analysis.
	-spectrum_files (mandatory): The MGF file.
	-searchParameters (mandatory): The search parameters file to be used for the search.
	-inputPredicted (optional): The non-canonical gene models, these are a set of GFF/FASTA files. The GFF and FASTA are separated by ';' and the pairs are separated by '##'

search.txt content example (see here for more info https://code.google.com/p/compomics-utilities/wiki/IdentificationParametersCLI):
	-prec_ppm 1 -prec_tol 5 -fixed_mods "carbamidomethyl c" -variable_mods "oxidation of m"
	

