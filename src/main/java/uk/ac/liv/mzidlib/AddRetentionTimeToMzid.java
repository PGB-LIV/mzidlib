package uk.ac.liv.mzidlib;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.Scan;
import uk.ac.ebi.jmzml.model.mzml.ScanList;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;
import uk.ac.liv.mzidlib.util.Java7Mapper;
import uk.ac.liv.mzidlib.util.Java7Optional;
import uk.ac.liv.mzidlib.util.Java7Predicate;
import static uk.ac.liv.mzidlib.util.Java7Stream.stream;
import uk.ac.liv.mzidlib.util.MzidLibUtils;

/**
 *
 *
 * @author Fawaz Ghali
 */
public class AddRetentionTimeToMzid {

    private String inputMzid;
    private String inputSourceFile;
    private String outputFile;
    private File mgfFile;
    private MgfFile inputMGFParser;
    private MzMLUnmarshaller mzMLUnmarshallerunmarshaller;
    private Map<String, String> spectrumRT = new HashMap<>();

    private Ms2Query ms2Query;
    private MzIdentMLUnmarshaller mzIdentMLUnmarshaller;
    private Cv psiCV;
    private Cv unitCV;
    private List<SpectrumIdentificationList> spectrumIdentificationList = new ArrayList<>();
    //metadata
    private AnalysisSoftwareList analysisSoftwareList;
    private AuditCollection auditCollection;
    private Provider provider;
    private AnalysisProtocolCollection analysisProtocolCollection;
    private CvList cvList;
    private AnalysisCollection analysisCollection;
    private Inputs inputs;
    private MzidLibUtils mzidLibUtils;

    public static void main(String args[]) {
        new AddRetentionTimeToMzid(args[0], args[1], args[2]);
    }

    public AddRetentionTimeToMzid(String mzidIn, String inputSourceFile, String mzidOut) {
        this.inputMzid = mzidIn;
        this.outputFile = mzidOut;
        this.inputSourceFile = inputSourceFile;
        mgfFile = new File(inputSourceFile);
        mzidLibUtils = new MzidLibUtils();
        init();
    }

    private void init() {

        insertDetailsIntoMzid();
        writeMzidFile(outputFile);
    }

    private void insertDetailsIntoMzid() {
        try {
            mzIdentMLUnmarshaller = new MzIdentMLUnmarshaller(new File(inputMzid));
            if (inputSourceFile.endsWith(".mgf")) {
                inputMGFParser = new MgfFile(mgfFile);
            } else {
                File xmlFile = new File(inputSourceFile);
                mzMLUnmarshallerunmarshaller = new MzMLUnmarshaller(xmlFile);
                MzMLObjectIterator<Spectrum> spectrumIterator = mzMLUnmarshallerunmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);
                while (spectrumIterator.hasNext()) {
                    Spectrum spectrum = spectrumIterator.next();
                    ScanList scanList = spectrum.getScanList();

                    List<Scan> scans = scanList.getScan();
                    for (Scan scan : scans) {
                        List<CVParam> cvParamList = scan.getCvParam();
                        for (CVParam cvParam : scan.getCvParam()) {
                            if (cvParam.getAccession().equals("MS:1000016")) {
                                spectrumRT.put(spectrum.getId(), cvParam.getValue());
                                break;
                            }
                        }
                    }
                }
            }

            //cvList = mzIdentML.getCvList();
            cvList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.CvList);
            //analysisSoftwareList = mzIdentML.getAnalysisSoftwareList();
            analysisSoftwareList = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisSoftwareList);
            //auditCollection = mzIdentML.getAuditCollection();
            auditCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AuditCollection);
            //provider = mzIdentML.getProvider();
            provider = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Provider);
            // analysisProtocolCollection = mzIdentML.getAnalysisProtocolCollection();
            analysisProtocolCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisProtocolCollection);
            //analysisCollection = mzIdentML.getAnalysisCollection();
            analysisCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.AnalysisCollection);
            //inputs = mzIdentML.getDataCollection().getInputs();
            inputs = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.Inputs);
            // searchDatabase_Ref = inputs.getSearchDatabase().get(0).getId();

            Iterator<Cv> iterCv = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.CV);
            while (iterCv.hasNext()) {
                Cv cv = iterCv.next();
                if (cv.getUri().toLowerCase().contains("psi")) {
                    psiCV = cv;
                } else if (cv.getUri().toLowerCase().contains("unit")) {
                    unitCV = cv;
                }
            }

        } catch (JMzReaderException ex) {
            ex.printStackTrace();
        } catch (OutOfMemoryError error) {

            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + error.getMessage() + "."
                    + "\nPlease see the reference guide at 05 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);

        }

    }

    // Write the new data into a file
    public void writeMzidFile(String csvFileName) {
        try {
            String outFile = csvFileName;
            Writer writer = new FileWriter(outFile);

            MzIdentMLMarshaller marshaller;
            marshaller = new MzIdentMLMarshaller();

            writer.write(marshaller.createXmlHeader() + "\n");

            String mzID = mzIdentMLUnmarshaller.getMzIdentMLId();
            if (mzID != null) {
                writer.write(marshaller.createMzIdentMLStartTag(mzID) + "\n");
            } else {
                writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");
            }

            if (cvList != null) {
                marshaller.marshal(cvList, writer);
            }
            writer.write("\n");
            if (analysisSoftwareList != null) {
                AnalysisSoftware analysisSoftware = new AnalysisSoftware();
                Date date = new Date();
                SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
                analysisSoftware.setName(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
                analysisSoftware.setId(this.getClass().getSimpleName() + "_" + dateFormat.format(date));
                Param param = new Param();
                param.setParam(mzidLibUtils.makeCvParam("MS:1002237", "mzidLib", psiCV));
                analysisSoftware.setSoftwareName(param);
                analysisSoftwareList.getAnalysisSoftware().add(analysisSoftware);
                marshaller.marshal(analysisSoftwareList, writer);
            }
            writer.write("\n");

            if (provider != null) {
                marshaller.marshal(provider, writer);
            }
            writer.write("\n");

            if (auditCollection != null) {
                marshaller.marshal(auditCollection, writer);
            }
            writer.write("\n");

            SequenceCollection sequenceCollection = mzIdentMLUnmarshaller.unmarshal(MzIdentMLElement.SequenceCollection);

            if (sequenceCollection != null) {
                marshaller.marshal(sequenceCollection, writer);
            }
            writer.write("\n");

            if (analysisCollection != null) {
                marshaller.marshal(analysisCollection, writer);
            }
            writer.write("\n");

            if (analysisProtocolCollection != null) {
                marshaller.marshal(analysisProtocolCollection, writer);
            }
            writer.write("\n");

            writer.write(marshaller.createDataCollectionStartTag() + "\n");

            writer.write("\n");

            if (inputs != null) {
                marshaller.marshal(inputs, writer);
            }
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");

            String spectrumIdentificationListRef = "";
            if (analysisCollection.getSpectrumIdentification().size() > 0) {
                spectrumIdentificationListRef = analysisCollection.getSpectrumIdentification().get(0).getSpectrumIdentificationListRef();
            }
            SpectrumIdentificationList siList;
            siList = new SpectrumIdentificationList();

            siList.setId(spectrumIdentificationListRef);

            Iterator<FragmentationTable> iterFragmentationTable = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.FragmentationTable);
            while (iterFragmentationTable.hasNext()) {

                FragmentationTable fr = iterFragmentationTable.next();
                siList.setFragmentationTable(fr);
            }

            Iterator<SpectrumIdentificationResult> iterSpectrumIdentificationResult = mzIdentMLUnmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.SpectrumIdentificationResult);
            int srProcessed = 0;
            int srNoMatch = 0;            
            while (iterSpectrumIdentificationResult.hasNext()) {
                SpectrumIdentificationResult sr = iterSpectrumIdentificationResult.next();
                String spectrumID = sr.getSpectrumID();

                if (inputSourceFile.endsWith(".mgf")) {
                    String spectrumIndex = spectrumID.substring(6);
                    Integer index1 = Integer.valueOf(spectrumIndex);
                    ms2Query = inputMGFParser.getMs2Query(index1.intValue());
                    Java7Optional<SpectrumIdentificationItem> optional = stream(sr.getSpectrumIdentificationItem()).filter(new MzMatchingPredicate(ms2Query, 0.1)).any();
                    srProcessed++;
                    if (optional.isEmpty()) {
                        // If we are here, then there was no match within 0.1 m/z in the SIR
                        // Let's list the values within 2 m/z units, or throw an exception if there are not any..
                        List<SpectrumIdentificationItem> items = stream(sr.getSpectrumIdentificationItem()).filter(new MzMatchingPredicate(ms2Query, 2)).toList();                        
                        if (items.isEmpty()) {                            
                            throw new RuntimeException("No match to " + ms2Query.getPrecursorMZ() + " within 0.1 m/z for spectrum with index " + index1.intValue());
                        } else {
                            List<String> mzs = stream(items).map(mzMapperToString).toList();                            
                            System.out.println("Warning: No match to " + ms2Query.getPrecursorMZ() 
                                    + " within 0.1 m/z for spectrum with index: " 
                                    + index1.intValue() + ". Within 2 m/z: " 
                                    + String.join(" ", mzs));
                            srNoMatch++;
                        }                        
                    }                  

                    CvParam rtParam = mzidLibUtils.makeCvParam("MS:1000016", "scan start time", psiCV, "" + ms2Query.getRetentionTime());
                    rtParam.setUnitAccession("UO:0000010");
                    rtParam.setUnitName("second");
                    sr.getCvParam().add(rtParam);
                    sr.getCvParam().add(mzidLibUtils.makeCvParam("MS:1000796", "spectrum title", psiCV, "" + ms2Query.getTitle()));
                } else {
                    String rt = spectrumRT.get(spectrumID);
                    if (rt != null) {
                        sr.getCvParam().add(mzidLibUtils.makeCvParam("MS:1000016", "scan start time", psiCV, "" + rt));
                    }
                }

                siList.getSpectrumIdentificationResult().add(sr);
            }
            
            System.out.println("Processed: " + srProcessed + ". No match within 0.1, match within 2.0: " + srNoMatch);

            marshaller.marshal(siList, writer);
            writer.write("\n");

            writer.write(marshaller.createProteinDetectionListStartTag("PDL_1", null) + "\n");
            writer.write(marshaller.createProteinDetectionListClosingTag() + "\n");
            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");
            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

            writer.close();

            System.out.println("Output written to " + outFile);

        } catch (JMzReaderException ex) {
            ex.printStackTrace();
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = this.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

    private static class MzMatchingPredicate implements Java7Predicate<SpectrumIdentificationItem> {
        private final Ms2Query query;
        private final double tolerance;
        public MzMatchingPredicate(Ms2Query query, double tolerance) {
            this.query = query;
            this.tolerance = tolerance;
        }

        @Override
        public boolean test(SpectrumIdentificationItem item) {
            double delta = Math.abs(item.getExperimentalMassToCharge() - query.getPrecursorMZ());
            return delta <= tolerance;
        }
    }   
   
    private final Java7Mapper<SpectrumIdentificationItem, String> mzMapperToString = new Java7Mapper<SpectrumIdentificationItem, String>() {    
        @Override
        public String map(SpectrumIdentificationItem item) {
            return String.valueOf(item.getExperimentalMassToCharge());
        }
    };    
}