/*
 * Date: 31-Aug-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.search.SearchSingleMgfAction.java
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

package uk.ac.liv.mzidlib.search;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.RecursiveAction;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.lang.ArrayUtils;
import org.xml.sax.SAXException;

import uk.ac.ebi.jmzidml.model.utils.MzIdentMLVersion;
import uk.ac.liv.mzidlib.SearchGUICLI;
import uk.ac.liv.mzidlib.constants.CvConstants;
import uk.ac.liv.mzidlib.converters.Convert2MzidTask;
import uk.ac.liv.mzidlib.util.Utils;
import uk.ac.liv.mzidlib.writer.Omssa2mzidMzidContainer;
import uk.ac.liv.mzidlib.writer.Tandem2mzidMzidContainer;

/**
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 31-Aug-2017 09:06:03
 */
public class SearchSingleMgfAction extends RecursiveAction {

    private SearchGUICLI searchGuiCli;
    private String[] searchParams;

    public SearchSingleMgfAction(SearchGUICLI sgc, String[] sp) {
        this.searchGuiCli = sgc;
        System.arraycopy(sp, 0, this.searchParams, 0, sp.length);
    }

    @Override
    protected void compute() {

        runSearchGuiCli(this.searchParams);

        String newOutput = Utils.getCmdParameter(searchParams, "output_folder",
                                                 true);
        File outputDir = new File(newOutput);
        String[] outputListFiles = outputDir.list();
        String omssaFileName = "";
        String tandemFileName = "";
        String msgfFileName = "";

        if (outputListFiles != null) {
            for (String fn : outputListFiles) {
                if (fn.endsWith(".omx")) {
                    omssaFileName = newOutput + File.separator + fn;
                } else if (fn.endsWith(".t.xml")) {
                    tandemFileName = newOutput + File.separator + fn;
                } else if (fn.endsWith(".msgf.mzid")) {
                    msgfFileName = newOutput + File.separator + fn;
                }
            }
        }

        File tandemMzidFile;
        File omssaMzidFile;
        try {
            if (!tandemFileName.isEmpty()) {
                Tandem2mzidMzidContainer tandemContainer
                        = new Tandem2mzidMzidContainer(tandemFileName,
                                                       CvConstants.FASTA_FORMAT
                                                       .getAccession(),
                                                       CvConstants.MASCOT_MGF_FORMAT
                                                       .getAccession(),
                                                       Boolean.FALSE,
                                                       null, "\\S+",
                                                       MzIdentMLVersion.Version_1_2);
                String tandemOut = tandemFileName.substring(0, tandemFileName
                                                            .lastIndexOf("."))
                        + "_tandem.mzid";
                Convert2MzidTask convertTandemTask = new Convert2MzidTask(tandemOut,
                                                         tandemContainer);
                convertTandemTask.fork();
                tandemMzidFile = convertTandemTask.join();
            }

            if (!omssaFileName.isEmpty()) {
                Omssa2mzidMzidContainer omssaContainer
                        = new Omssa2mzidMzidContainer(omssaFileName,
                                                      Boolean.FALSE,
                                                      "REVERSED",
                                                      null,
                                                      null,
                                                      MzIdentMLVersion.Version_1_2);
                String omssaOut = omssaFileName.substring(0, omssaFileName
                                                          .lastIndexOf("."))
                        + "_omssa.mzid";
                Convert2MzidTask convertOmssaTask
                        = new Convert2MzidTask(omssaOut, omssaContainer);
                convertOmssaTask.fork();
                omssaMzidFile = convertOmssaTask.join();
            }
            
            
        } catch (SAXException ex) {
            Logger.getLogger(SearchSingleMgfAction.class.getName())
                    .log(Level.SEVERE, null, ex);
        } catch (ParserConfigurationException ex) {
            Logger.getLogger(SearchSingleMgfAction.class.getName())
                    .log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(SearchSingleMgfAction.class.getName())
                    .log(Level.SEVERE, null, ex);
        }

    }

    public void runSearchGuiCli(String[] searchParamters) {
        String searchGuiPath = getSearchGuiFile();

        try {
            System.out.println(
                    "Running SearchGUI search function using the following command:");

            String[] javaCli = {"java", "-Djava.awt.headless=true", "-cp",
                searchGuiPath, "eu.isas.searchgui.cmd.SearchCLI"};
            String[] args = (String[]) ArrayUtils.addAll(javaCli,
                                                         searchParamters);

            for (int i = 0; i < args.length; i++) {
                System.out.print(args[i] + " ");

            }
            System.out.println();

            ProcessBuilder pb = new ProcessBuilder(args);
            pb.redirectErrorStream(true);
            final Process process = pb.start();
            InputStream is = process.getInputStream();
            InputStreamReader isr = new InputStreamReader(is,
                                                          StandardCharsets.UTF_8);
            try (BufferedReader br = new BufferedReader(isr)) {
                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);

                }
                try {
                    process.waitFor();

                } catch (InterruptedException ex) {
                    process.destroy();
                    is.close();
                }

            }
            System.out.println("SearchGUI search function finished.");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Returns the path to the jar file's parent folder.
     *
     * @return the path to the jar file's parent folder
     */
    public String getSearchGuiPath() {
        return new File(SearchSingleMgfAction.class.getProtectionDomain()
                .getCodeSource()
                .getLocation().getPath()).getParent() + File.separator
                + "searchgui";
    }

    /**
     * Returns SearchGUI jar file.
     *
     * @return SearchGUI jar file.
     */
    public String getSearchGuiFile() {
        String path = null;
        try {
            File dir = new File(URLDecoder.decode(getSearchGuiPath(), "UTF-8"));
            System.out.println("SearchGUI dir: " + dir);

            File[] matches = dir.listFiles(new FilenameFilter() {

                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".jar");
                }

            });
            path = matches[0].getAbsolutePath();

        } catch (UnsupportedEncodingException ex) {
            ex.printStackTrace();
        }
        return path;
    }

}
