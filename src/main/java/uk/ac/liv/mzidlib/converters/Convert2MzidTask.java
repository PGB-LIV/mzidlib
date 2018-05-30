/*
 * Date: 01-Sep-2017
 * Author: Da Qi
 * File: uk.ac.liv.mzidlib.converters.Convert2MzidTask.java
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

package uk.ac.liv.mzidlib.converters;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.RecursiveTask;
import java.util.logging.Level;
import java.util.logging.Logger;

import uk.ac.liv.mzidlib.writer.MzidContainer;
import uk.ac.liv.mzidlib.writer.MzidWriter;

/**
 * Convert file to mzid task.
 *
 * @author Da Qi
 * @institute University of Liverpool
 * @time 01-Sep-2017 10:01:04
 */
public class Convert2MzidTask extends RecursiveTask<File> {

    private final String outputFileName;
    private final MzidContainer mzidContainer;

    public Convert2MzidTask(String out, MzidContainer mc) {
        this.outputFileName = out;
        this.mzidContainer = mc;
    }

    @Override
    protected File compute() {
        try {
            MzidWriter.write(outputFileName, mzidContainer);

        } catch (IOException ex) {
            Logger.getLogger(Convert2MzidTask.class.getName())
                    .log(Level.SEVERE, null, ex);
        }
        return new File(this.outputFileName);
    }

}
