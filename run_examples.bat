java -Xms1024m -jar "mzidentml-lib.jar" Csv2mzid example_files\iprg_omssa.csv test_outputs\iprg_omssa_from_csv.mzid -paramsFile  example_files\iprg_omssa_params.csv -cvAccessionForPSMOrdering MS:1001328 -decoyRegex RRRRR

java -Xms1024m -jar "mzidentml-lib.jar" Omssa2mzid example_files\iprg_omssa.omx test_outputs\iprg_omssa_include_fragmentation.mzid -outputFragmentation true -decoyRegex RRRRR

java -Xms1024m -jar "mzidentml-lib.jar" Tandem2mzid example_files\iprg2008_tandem.xml test_outputs\iprg2008_tandem_no_frags.mzid -outputFragmentation false -decoyRegex RRRRR

java -Xms1024m -jar "mzidentml-lib.jar" FalseDiscoveryRate test_outputs\iprg2008_tandem_no_frags.mzid test_outputs\iprg2008_tandem_no_frags_fdr.mzid -decoyRegex Rnd -decoyValue 3 -decoyValue 3 -cvTerm "MS:1001330" -betterScoresAreLower true

java -Xms1024m -jar "mzidentml-lib.jar" Threshold test_outputs\iprg2008_tandem_no_frags_fdr.mzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold.mzid -isPSMThreshold true -cvAccessionForScoreThreshold MS:1001874 -threshValue 0.01  -betterScoresAreLower true -deleteUnderThreshold true

java -Xms1024m -jar "mzidentml-lib.jar" InsertMetaDataFromFasta test_outputs\iprg2008_tandem_no_frags_fdr_threshold.mzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated.mzid -fastaFile example_files/MGI_UniProtMouse_20071203_iPRG2008_fixed__tdr.fasta -accessionSplitRegex "/ /"

java -Xms1024m -jar "mzidentml-lib.jar" ProteoGrouper test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated.mzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated_grouped.mzid -cvAccForSIIScore MS:1001874 -logTransScore true -requireSIIsToPassThreshold true  -verboseOutput false

java -Xms1024m -jar "mzidentml-lib.jar" AddEmpaiToMzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated_grouped.mzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated_grouped_empai.mzid  -fastaFile example_files/MGI_UniProtMouse_20071203_iPRG2008_fixed__tdr.fasta -accessionSplitRegex "/ /" -verboseOutput false

java -Xms1024m -jar "mzidentml-lib.jar" Mzid2Csv test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated_grouped_empai.mzid test_outputs\iprg2008_tandem_no_frags_fdr_threshold_annotated_grouped_empai.csv  -exportType exportProteinGroups

