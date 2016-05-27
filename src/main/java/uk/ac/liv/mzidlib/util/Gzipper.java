package uk.ac.liv.mzidlib.util;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author fghali
 */
public class Gzipper {

    public static File extractFile(File zipped_file) {
        GZIPInputStream gin = null;
        File outFile = null;
        try {

            OutputStream fos = null;
            gin = new GZIPInputStream(new FileInputStream(zipped_file));
            outFile = new File(zipped_file.getParent(), "tmp_" + zipped_file.getName().replaceAll("\\.gz$", ""));
            fos = new FileOutputStream(outFile);
            byte[] buf = new byte[100000];
            int len;
            while ((len = gin.read(buf)) > 0) {
                fos.write(buf, 0, len);
            }
            fos.close();

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                gin.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        return outFile;

    }

//    public static void deleteFile(File file) {
//
//
//        boolean success = file.delete();
//        if (!success) {
//            System.out.println("Deletion failed.");
//        } else {
//            System.out.println("File deleted.");
//
//        }
//
//    }
    public static void compressFile(File file) throws InterruptedException {
        try {
            DateFormat format = new SimpleDateFormat("yyyy_MM_dd_hh_mm_ss");
            String timeStamp = format.format(new Date());
            OutputStream fos = new FileOutputStream(file + "_" + timeStamp + ".gz");
            GZIPOutputStream gzos = new GZIPOutputStream(fos);

            InputStream fin = new FileInputStream(file);
            InputStream in = new BufferedInputStream(fin);
            byte[] buffer = new byte[1024];
            int i;
            while ((i = in.read(buffer)) >= 0) {
                gzos.write(buffer, 0, i);
            }
            System.out.println("the file is in now gzip format");

            in.close();
            in = null;

            fin.close();
            fin = null;

            gzos.flush();
            gzos.close();
            gzos = null;

            fos.flush();
            fos.close();
            fos = null;

            System.gc();
            // delete the input file.

            //Thread.sleep(1000);
            file.deleteOnExit();
        } catch (IOException e) {
            String methodName = Thread.currentThread().getStackTrace()[1].getMethodName();
            String className = Gzipper.class.getName();
            String message = "The task \"" + methodName + "\" in the class \"" + className + "\" was not completed because of " + e.getMessage() + "."
                    + "\nPlease see the reference guide at 02 for more information on this error. https://code.google.com/p/mzidentml-lib/wiki/CommonErrors ";
            System.out.println(message);

        }

    }
}
