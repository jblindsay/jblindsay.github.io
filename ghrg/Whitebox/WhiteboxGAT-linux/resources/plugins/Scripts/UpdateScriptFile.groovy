/*
 * Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
import java.io.File
import java.net.URL
import javax.swing.*
import groovy.io.FileType
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities
import whiteboxgis.WhiteboxGui

/* The following lines are necessary for the script 
   to be recognized as a menu extension. */
def parentMenu = "Tools"
def menuLabel = "Update Scripts From Repository"

public class Downloader {
	String site = ""
	String filename = ""
	WhiteboxGui pluginHost;
	
	public Downloader(WhiteboxGui pluginHost, String filename, String site) {
		this.pluginHost = pluginHost
		this.site = site
		this.filename = filename
	}

	public void setSite(String site) {
		this.site = site
	}

	public void setFilename(String filename) {
		this.filename = filename
	}
	
    public int download(String fn, long localLastModified) throws Exception {
    
	    try {
	    	
	        URL url=new URL(site);
	        HttpURLConnection connection = (HttpURLConnection) url.openConnection();
	        connection.setReadTimeout(100000)

			long lastModified = connection.getLastModified()

			if (localLastModified > lastModified) {
				if (pluginHost.showFeedback("The locally modifed file is newer than the file in the code repository." +
				"\nAre you sure you want to replace the local file?", 
				JOptionPane.YES_NO_OPTION, 
				JOptionPane.QUESTION_MESSAGE) == JOptionPane.NO_OPTION) {
                	return -1
				}
			}
	        
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
	            pluginHost.updateProgress("Downloading ${fn}...", (int) percent)
	            // check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return -1
				}
            }
            bout.close();
            inStream.close();
            pluginHost.updateProgress("Progress:", 0)
            return 0
	    } catch (FileNotFoundException e) {
		    return 1
	    } catch(SocketTimeoutException e) {
		    return 2
	    } catch(IOException e) {
		    return 3
	    } catch(Exception e) {
	    	pluginHost.logException("Error in " + descriptiveName, e)
		    return 4
	    }
	}
}

def updateScript = { scriptName ->
	
	WhiteboxGui wg = (WhiteboxGui)(pluginHost)

	// First get the existing file. In case something goes wrong 
	// with the download, you want to be able to bring back the 
	// old file
	String pathSep = File.separator
	String scriptDir = "${pluginHost.resourcesDirectory}plugins${pathSep}Scripts${pathSep}"
	String localFileName = "${scriptDir}${scriptName}"
	def localFile = new File(localFileName)
	long localFileModifiedTime = localFile.lastModified()
	String fileContents = localFile.text
	
	// Figure out the name of the remote file
	//String remoteDir = "https://whitebox-geospatial-analysis-tools.googlecode.com/svn/trunk/WhiteboxGIS/resources/plugins/Scripts/"
	String remoteDir = "https://raw.githubusercontent.com/jblindsay/whitebox-geospatial-analysis-tools/master/WhiteboxGIS/resources/plugins/Scripts"
	String downloadArtifact = "${remoteDir}/${scriptName}"

	// Download the remote file
	//String outFileName = "/Users/johnlindsay/Documents/${scriptName}" 
	Downloader d = new Downloader(wg, localFileName, downloadArtifact) 
	pluginHost.showFeedback(localFileName)
	pluginHost.showFeedback(downloadArtifact)
	int ret = d.download(scriptName, localFileModifiedTime)
	//pluginHost.showFeedback(downloadArtifact)
	if (ret == 0) {
		// see if there is a help file for this script that you can download
		String localHelpDir = "${pluginHost.resourcesDirectory}Help${pathSep}"
		//String remoteHelpDir = "https://whitebox-geospatial-analysis-tools.googlecode.com/svn/trunk/WhiteboxGIS/resources/Help/"
		String remoteHelpDir = "https://raw.githubusercontent.com/jblindsay/whitebox-geospatial-analysis-tools/master/WhiteboxGIS/resources/Help/"
	
		String helpFileName = scriptName.substring(0, scriptName.indexOf(".")) + ".html"
		downloadArtifact = "${remoteHelpDir}${helpFileName}"
		localFileName = "${localHelpDir}${helpFileName}"
		d.setFilename(localFileName)
		d.setSite(downloadArtifact)
		if (d.download(scriptName, -1) == 0) {
			
		}
				
		pluginHost.showFeedback("Successful update")
	} else {
		pluginHost.showFeedback("The update was unsuccessful")
		pluginHost.returnData(fileContents)
		localFile.delete()
		localFile.withWriter { out ->
		    out.write(fileContents)
		}
	}
}

def localScriptList = {
	def list = []
	String pathSep = File.separator
	String scriptDir = "${pluginHost.resourcesDirectory}plugins${pathSep}Scripts${pathSep}"
	def localDir = new File(scriptDir)
	localDir.eachFileRecurse(FileType.FILES) { file ->
	  list << file.getName()
	}
	list
}

def updateAllScripts = {

	def localScripts = localScriptList()
	WhiteboxGui wg = (WhiteboxGui)(pluginHost)
	String pathSep = File.separator
	String scriptDir = "${pluginHost.resourcesDirectory}plugins${pathSep}Scripts${pathSep}"
	//String dir = "https://whitebox-geospatial-analysis-tools.googlecode.com/svn/trunk/WhiteboxGIS/resources/plugins/Scripts/"
	String dir = "https://raw.githubusercontent.com/jblindsay/whitebox-geospatial-analysis-tools/master/WhiteboxGIS/resources/plugins/Scripts/"

	String localHelpDir = "${pluginHost.resourcesDirectory}Help${pathSep}"
	//String remoteHelpDir = "https://whitebox-geospatial-analysis-tools.googlecode.com/svn/trunk/WhiteboxGIS/resources/Help/"
	String remoteHelpDir = "https://raw.githubusercontent.com/jblindsay/whitebox-geospatial-analysis-tools/master/WhiteboxGIS/resources/Help/"
	
	URL url = new URL(dir)
	HttpURLConnection connection = (HttpURLConnection) url.openConnection();
	long lastModified = connection.getLastModified()
	connection.setReadTimeout(100000)
	int filesize = connection.getContentLength();
	//java.io.BufferedInputStream inStream = new java.io.BufferedInputStream(connection.getInputStream());
	BufferedReader inStream = new BufferedReader(new InputStreamReader(url.openStream()));
	
	def scriptNames = new ArrayList<String>()
	String line;
	while ((line = inStream.readLine()) != null) {
	    if (line.contains("a href=")) {
	    	int i = line.indexOf(">", line.indexOf(">") + 1)
	    	int j = line.indexOf("<", i)
	    	def content = line.substring(i + 1, j).trim()
	    	if (!content.equals("..") && !content.isEmpty()) {
	    		scriptNames.add(content)
	    	}
	    }
	}
	inStream.close()
	
	def listOfUpdatedScripts = []
	def listOfNewScripts = []
	scriptNames.each() { n -> 
		def downloadArtifact = "${dir}${n}"
		def localFileName = "${scriptDir}${n}"
		if (localScripts.find { j -> j.equals(n) } == null) {
			Downloader d = new Downloader(wg, localFileName, downloadArtifact) 
			int ret = d.download(n, -1)
			if (ret == 0) {
				listOfNewScripts.add(n)

				// see if there is a help file for this script that you can download
				String helpFileName = n.substring(0, n.indexOf(".")) + ".html"
				downloadArtifact = "${remoteHelpDir}${helpFileName}"
				localFileName = "${localHelpDir}${helpFileName}"
				d.setFilename(localFileName)
				d.setSite(downloadArtifact)
				if (d.download(n, -1) == 0) {
					
				}
			}
		} else {
			// the script is in both the repo and locally. Is the repo newer?
			url = new URL(downloadArtifact)
			connection = (HttpURLConnection) url.openConnection();
			lastModified = connection.getLastModified()
			def localFile = new File(localFileName)
			long localFileModifiedTime = localFile.lastModified()
			if (localFileModifiedTime < lastModified) {
				//println "${n} ${localFileModifiedTime} ${lastModified}"
				
				Downloader d = new Downloader(wg, localFileName, downloadArtifact) 
				int ret = d.download(n, -1)
				if (ret == 0) {
					listOfUpdatedScripts.add(n)
	
					// see if there is a help file for this script that you can download
					String helpFileName = n.substring(0, n.indexOf(".")) + ".html"
					downloadArtifact = "${remoteHelpDir}${helpFileName}"
					localFileName = "${localHelpDir}${helpFileName}"
					d.setFilename(localFileName)
					d.setSite(downloadArtifact)
					if (d.download(n, -1) == 0) {
						
					}
				}
			}
		}
	}
	if (listOfUpdatedScripts.size() > 0 || listOfNewScripts.size() > 0) {
		def output = new StringBuilder()
		output.append("<!DOCTYPE html>")
		output.append('<html lang="en">')
		output.append("<head>")
        output.append("<title>Script Update Summary</title>").append("\n")
        output.append("<style>")
		output.append("table {margin-left: 15px;} ")
		output.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
		output.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
		output.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
		output.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
		output.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
		output.append("caption {font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
		output.append(".numberCell { text-align: right; }") 
        output.append("</style></head>").append("\n")
        output.append("<body><h1>Script Update Summary</h1>").append("\n")

		if (listOfUpdatedScripts.size() > 0) {
			output.append("<p>The following scripts have newer version in the code repository than \n" + 
			"your local script versions and have been updated. Note that if any of the updated \n" + 
			"scripts require changes made to the Whitebox GAT program, the updates may result \n" +
			"in apparent errors in the scripts. If this is the case, you should check for a newer \n" +
			"version of Whitebox GAT, if there is one.</p>")
			output.append("<ul>")
			listOfUpdatedScripts.each() { s -> output.append("<li>${s}</li>") }
			output.append("</ul>")
		}
		if (listOfNewScripts.size() > 0) {
			output.append("<p>New Scripts:</p>")
			output.append("<ul>")
			listOfNewScripts.each() { s -> output.append("<li>${s}</li>") }
			output.append("</ul>")
			output.append("<p>You will need to restart Whitebox GAT to use these new scripts.</p>")
			pluginHost.returnData(output.toString())
		}

		ret.append("</body></html>")
	} else {
		pluginHost.showFeedback("No new scripts were found in the code repository.")
	}
}

if (args.length >= 1) {
	updateScript(args[0])
} else {
	updateAllScripts()
}
						
