package uk.ac.liv.mzidlib.util;

import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.FileStore;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

// The idea is to see whether the underlying file system is one of a
// network kind (NFS, SMB, CIFS, OpenAFS, GlusterFS) and then to copy
// the file to the local scratch space for faster processing. 

// This implementation is primarily NFS (any version), as this is the
// most common version of networked file system, and it is flagged as
// the one even in the presence of other network enabled file system
// through masking.

public class FileHandler {
    
    private static final String NFSFileSystemPrefix = "nfs";

    private FileHandler()
    {

    }

    /**
     * This method copies the given file to a local drive if the given file is
     * in a network file system. The calling method should check for null object.
     *
     * @param filePath     as {@link String} the path of the file to be processed
     * @param processLocal decides whether the file should be copied to local drive
     * @return Locally available {@link File} for the given path
     */
    public static File handleFile(final String filePath, final boolean processLocal)
    {
        try {
            File originalFile = new File(filePath);
            if (processLocal) {
                URI rootURI = new URI("file:///");
                Path rootPath = Paths.get(rootURI);
                Path dirPath = rootPath.resolve(filePath);
                FileStore dirFileStore = Files.getFileStore(dirPath);
                File mzidFile = new File(originalFile.getName());
                if (dirFileStore.type().contains(NFSFileSystemPrefix)){
                    Files.copy(originalFile.toPath(), mzidFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                } else {
                    mzidFile = originalFile;
                }
                return mzidFile;
            } else {
                return originalFile;
            }
        } catch (URISyntaxException e) {
            // TODO: Not sure mzid has a builtin  logger. If so, needs logging. 
            return null;
        } catch (IOException e) {
            // TODO: Not sure mzid has a builtin  logger. If so, needs logging. 
            return null;
        }
    }

}
