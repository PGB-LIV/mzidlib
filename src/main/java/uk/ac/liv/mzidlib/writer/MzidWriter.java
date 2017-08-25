/*
 * Date: 23-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.writer.MzidWriter.java
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

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;

import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

/**
 * MzIdentML file writer class.
 * @author Da Qi
 *     University of Liverpool
 * @time 23-Aug-2017 14:24:47
 */
public class MzidWriter {

    /**
     * Static generic MzIdentML writer.
     * @param fileName output MzIdentML file name
     * @param mc MzIdentML elements container
     * @throws FileNotFoundException file not found exception
     * @throws IOException IO exception
     */
    public static void write(String fileName, MzidContainer mc)
            throws FileNotFoundException, IOException {

        OutputStream fos = new FileOutputStream(fileName);

        try (Writer writer = new OutputStreamWriter((fos),
                                                    StandardCharsets.UTF_8)) {
            MzIdentMLMarshaller marshaller = new MzIdentMLMarshaller();

            writer.write(marshaller.createXmlHeader() + "\n");

            // mzIdentML start tag
            writer.write(marshaller.createMzIdentMLStartTag("12345") + "\n");

            marshaller.marshal(mc.getCvList(), writer);
            writer.write("\n");

            marshaller.marshal(mc.getAnalysisSoftwareList(), writer);
            writer.write("\n");

            marshaller.marshal(mc.getProvider(), writer);
            writer.write("\n");

            marshaller.marshal(mc.getAuditCollection(), writer);
            writer.write("\n");

            //m.marshall(analysisSampleCollection, writer);     //TODO - complete this part
            //writer.write("\n");
            
            marshaller.marshal(mc.getSequenceCollection(), writer);
            writer.write("\n");

            marshaller.marshal(mc.getAnalysisCollection(), writer);
            writer.write("\n");

            marshaller.marshal(mc.getAnalysisProtocolCollection(), writer);
            writer.write("\n");

            writer.write(marshaller.createDataCollectionStartTag() + "\n");
            marshaller.marshal(mc.getInputs(), writer);
            writer.write("\n");

            writer.write(marshaller.createAnalysisDataStartTag() + "\n");

            marshaller.marshal(mc.getSpectrumIdentificationList(), writer);

            marshaller.marshal(mc.getProteinDetectionList(), writer);

            writer.write(marshaller.createAnalysisDataClosingTag() + "\n");

            writer.write(marshaller.createDataCollectionClosingTag() + "\n");

            writer.write(marshaller.createMzIdentMLClosingTag());

        } catch (IOException ex) {
            String methodName = Thread.currentThread()
                    .getStackTrace()[1].getMethodName();
            String className = MzidWriter.class.getClass().getName();
            String message = "The task \"" + methodName + "\" in the class \""
                    + className + "\" was not completed because of "
                    + ex.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more "
                    + "information on this error. "
                    + "https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);
        }
    }

}
