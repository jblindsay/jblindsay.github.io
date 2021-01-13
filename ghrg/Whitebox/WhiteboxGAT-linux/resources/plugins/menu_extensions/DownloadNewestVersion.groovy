import java.io.File
import java.io.FileInputStream
import java.io.FileOutputStream
import java.io.IOException
import java.util.List
import java.util.zip.*
import java.net.HttpURLConnection;
import java.net.URL;
import javax.swing.*;
import whiteboxgis.WhiteboxGui
import org.apache.commons.io.FilenameUtils
import groovy.time.*

// The following lines are necessary for the script 
// to be recognized as a menu extension.
parentMenu = "Help"
menuLabel = "Check For Update"

try {
	def timeStart = new Date()
	WhiteboxGui wg = (WhiteboxGui)(pluginHost)

	if (wg.currentVersionNumber == null || wg.currentVersionNumber.isEmpty()) {
		if (wg.isVersionUpToDate()) {
			wg.showFeedback("Your version of Whitebox GAT is current.")
			return
		}
	} else if (wg.currentVersionNumber <= wg.versionNumber) {
		wg.showFeedback("Your version of Whitebox GAT is current.")
		return
	} else {
		if (wg.showFeedback("There is a newer version. Would you like to update?", 
		    JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE) == 1) {
			return
		}
	}
	
	String downloadArtifact = wg.updateDownloadArtifact
	
	// ask the user for the directory into which the file should be saved.
	JFileChooser fc = new JFileChooser()
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY)
    String outFile = FilenameUtils.getBaseName(downloadArtifact) + ".zip"
    fc.setSelectedFile(new File(outFile));
    //appDirectory = wg.getApplicationDirectory() + File.separator
    //fc.setCurrentDirectory(new File(appDirectory));
    
    int result = fc.showSaveDialog(wg);
    File file = null;
    if (result == JFileChooser.APPROVE_OPTION) {
        String zipFileName = fc.getSelectedFile(); 
		File zipFile = new File(zipFileName)
		if (zipFile.exists()) {
			zipFile.delete()
			println "File already exists. Deleting old version..."
		}
	
		Downloader d = new Downloader(wg, zipFileName)
		if (d.doSomething()) {
			wg.updateProgress("Unzipping file...", 0)
			
			extractFolder(zipFileName)
		
			// once extracted, delete the original zip file
			if (zipFile.exists()) {
				zipFile.delete()
			}
			
			TimeDuration duration = TimeCategory.minus(new Date(), timeStart)
			wg.showFeedback("Operation completed in ${duration} seconds.")
	
			wg.updateProgress("Progress:", 0)
		} else {
			wg.showFeedback("Download was unsuccessful.")
		}
	} else {
    	wg.showFeedback("Operation Cancelled")
    }
} catch (IOException ioe) {
	println ioe.getMessage()
} catch (Exception e) {
	println e.getMessage()
}

// This method is from NeilMonday, http://stackoverflow.com/questions/981578/how-to-unzip-files-recursively-in-java
// Modified by J. Lindsay
static public void extractFolder(String zipFile) throws ZipException, IOException {
    int BUFFER = 2048;
    File file = new File(zipFile);

    ZipFile zip = new ZipFile(file);
    String newPath = zipFile.substring(0, zipFile.length() - 4);

    new File(newPath).mkdir();
    Enumeration zipFileEntries = zip.entries();

    // Process each entry
    while (zipFileEntries.hasMoreElements()) {
        // grab a zip file entry
        ZipEntry entry = (ZipEntry) zipFileEntries.nextElement();
        String currentEntry = entry.getName();
        File destFile = new File(newPath, currentEntry);
        File destinationParent = destFile.getParentFile();

        // create the parent directory structure if needed
        if (!destinationParent.exists()) {
        	destinationParent.mkdirs();
        }
		
        if (!entry.isDirectory()) {
            BufferedInputStream bis = new BufferedInputStream(zip.getInputStream(entry));
            int currentByte;
            // establish buffer for writing file
            byte[] data = new byte[BUFFER];

            // write the current file to disk
            FileOutputStream fos = new FileOutputStream(destFile);
            BufferedOutputStream dest = new BufferedOutputStream(fos,
            BUFFER);

            // read and write until last byte is encountered
            while ((currentByte = bis.read(data, 0, BUFFER)) != -1) {
                dest.write(data, 0, currentByte);
            }
            dest.flush();
            dest.close();
            bis.close();
        } else {
        	destFile.mkdirs();
        }

        if (currentEntry.endsWith(".zip")) {
            // found a zip file, try to open
            extractFolder(destFile.getAbsolutePath());
        }
    }
}

public class Downloader {
	String site = ""
	String filename = ""
	WhiteboxGui pluginHost;
	
	public Downloader(WhiteboxGui pluginHost, String filename) {
		this.pluginHost = pluginHost
		site = pluginHost.updateDownloadArtifact
		this.filename = filename
	}
	
    public boolean doSomething() throws Exception {
    
	    try {
	        URL url=new URL(site);
	        HttpURLConnection connection = (HttpURLConnection) url.openConnection();
	        int filesize = connection.getContentLength();
	        float totalDataRead = 0;
            java.io.BufferedInputStream inStream = new java.io.BufferedInputStream(connection.getInputStream());
            java.io.FileOutputStream fos = new java.io.FileOutputStream(filename);
            java.io.BufferedOutputStream bout = new BufferedOutputStream(fos,1024);
            byte[] data = new byte[1024];
            int i = 0;
            while((i = inStream.read(data, 0, 1024)) >= 0) {
	            totalDataRead=totalDataRead+i;
	            bout.write(data,0,i);
	            float percent = (totalDataRead * 100f) / filesize;
	            pluginHost.updateProgress("Downloading file...", (int) percent)
	            // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return false
				}
            }
            bout.close();
            inStream.close();
            pluginHost.updateProgress("Progress:", 0)
            return true
	    } catch(Exception e) {
	    	pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
		    pluginHost.logException("Error in DownloadNewestVersion", e)
		    return false
	    }
	}
}
