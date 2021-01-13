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
 
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.text.DecimalFormat
import java.util.Date
import java.util.GregorianCalendar
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours
import whitebox.structures.KdTree
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.LASReader.VariableLengthRecord
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "GetLasFileSummary"
def descriptiveName = "Get LAS File Summary"
def description = "Returns a summary of statistics and properties for a LAS file."
def toolboxes = ["LidarTools"]

public class GetLasFileSummary implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public GetLasFileSummary(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, this)	
		
			// Specifying the help file will display the html help
			// file in the help pane. This file should be be located 
			// in the help directory and have the same name as the 
			// class, with an html extension.
			sd.setHelpFile(name)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogMultiFile("Select the input LAS files", "Input LAS Files:", "LAS Files (*.las), LAS")
			
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	// The CompileStatic annotation can be used to significantly
	// improve the performance of a Groovy script to nearly 
	// that of native Java code.
	@CompileStatic
	private void execute(String[] args) {
	  	try {
		  	if (args.length != 1) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			DecimalFormat df = new DecimalFormat("###,###,###,###.###")
			DecimalFormat df2 = new DecimalFormat("##0.0##%")
			// read the input parameters
			final String inputFileString = args[0]
			
			String[] inputFiles = inputFileString.split(";")
	
			// check for empty entries in the inputFiles array
			int i = 0;
			for (String inputFile : inputFiles) {
				if (!inputFile.isEmpty()) {
					i++
				}
			}
	
			if (i != inputFiles.length) {
				// there are empty entries in the inputFiles array
				// remove them.
				ArrayList<String> temp = new ArrayList<>();
				for (String inputFile : inputFiles) {
					if (!inputFile.isEmpty()) {
						temp.add(inputFile)
					}
				}
				inputFiles = new String[i];
	    		temp.toArray(inputFiles);
			}
			
			
			final int numFiles = inputFiles.length
			i = 0
			for (String inputFile : inputFiles) {


				StringBuilder ret = new StringBuilder()
				ret.append("<!DOCTYPE html>")
				ret.append('<html lang="en">')
	
				ret.append("<head>")
	//            ret.append("<meta content=\"text/html; charset=iso-8859-1\" http-equiv=\"content-type\">")
	            ret.append("<title>LAS File Summary Statistics</title>").append("\n")
	            
	            ret.append("<style>")
				ret.append("table {margin-left: 15px;} ")
				ret.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
				ret.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
				ret.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
				ret.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
				ret.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
				ret.append("caption {font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
				ret.append(".numberCell { text-align: right; }") 
	            ret.append("</style></head>").append("\n")
	            ret.append("<body><h1>LAS File Summary Statistics</h1>").append("\n")

            
				ret.append("<p><b>File Name:</b> &nbsp ").append(new File(inputFile).getName()).append("</p>\n")
				
				LASReader las = new LASReader(inputFile)
				long totalPoints = las.getNumPointRecords()
				
				// get the day an month
				int year = las.getFileCreationYear()
				int dayOfYear = las.getFileCreationDay()
				GregorianCalendar gc = new GregorianCalendar()
				gc.set(GregorianCalendar.YEAR, year)
				gc.set(GregorianCalendar.DAY_OF_YEAR, dayOfYear)
				int dayOfMonth = gc.get(GregorianCalendar.DAY_OF_MONTH)
				int monthOfYear = gc.get(GregorianCalendar.MONTH) + 1
				String fileDate = dayOfMonth + "/" + monthOfYear + "/" + year

				ret.append("<p><b>File Creation Date:</b> &nbsp ").append(fileDate).append("<br>")
				ret.append("<b>Generating Software:</b> &nbsp ").append(las.getGeneratingSoftware()).append("<br>\n")
				ret.append("<b>LAS Version:</b> &nbsp ").append(las.getVersionMajor()).append(".").append(las.getVersionMinor()).append("</p>\n")
				
				int[] classHisto = new int[64]
				
				double minScanAngle = Double.POSITIVE_INFINITY
				double maxScanAngle = Double.NEGATIVE_INFINITY
				double minIntensity = Double.POSITIVE_INFINITY
				double maxIntensity = Double.NEGATIVE_INFINITY
				int numFirstRet = 0
				int numLastRet = 0
				int intermediateRet = 0
				int onlyRet = 0
				
				PointRecord point
				int progress
				int oldProgress = -1
				for (int n = 0; n < totalPoints; n++) {
					point = las.getPointRecord(n)
					byte scanAngle = point.getScanAngle()
					if (scanAngle < minScanAngle) { minScanAngle = scanAngle }
					if (scanAngle > maxScanAngle) { maxScanAngle = scanAngle }
					
					int intensity = point.getIntensity()
					if (intensity < minIntensity) { minIntensity = intensity }
					if (intensity > maxIntensity) { maxIntensity = intensity }
					
					int classification = point.getClassification()
					classHisto[classification]++
					
					byte retNum = point.getReturnNumber()
					byte numRets = point.getNumberOfReturns()

					if (numRets == 1) {
						onlyRet++
					} else if (retNum == 1) {
						numFirstRet++
					} else if (retNum == numRets) {
						numLastRet++
					} else {
						intermediateRet++
					}

					progress = (int)(100f * n / (totalPoints - 1))
            		if (progress != oldProgress) {
						pluginHost.updateProgress(progress)
            			oldProgress = progress
            			// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
            		}
				}

				double rangeX = las.getMaxX() - las.getMinX()
				double rangeY = las.getMaxY() - las.getMinY()
				double area = rangeX * rangeY
				double avgDensity = totalPoints / area

				ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
				ret.append("<caption>File Summary</caption>")
				
				ret.append("<tr><th>Statistic</th><th>Value</th></tr>")

				ret.append("<tr><td>Total Num. Points</td><td class=\"numberCell\">")
				ret.append(df.format(totalPoints)).append("</td></tr>").append("\n")

				ret.append("<tr><td>Avg. Density</td><td class=\"numberCell\">")
				ret.append(df.format(avgDensity)).append(" pts. m<sup>-2</sup></td></tr>").append("\n")

				ret.append("<tr><td>Min. X</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMinX())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Max. X</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMaxX())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Min. Y</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMinY())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Max. Y</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMaxY())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Min. Z</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMinZ())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Max. Z</td><td class=\"numberCell\">")
				ret.append(df.format(las.getMaxZ())).append("</td></tr>").append("\n")

				ret.append("<tr><td>Min. Intensity</td><td class=\"numberCell\">")
				ret.append(df.format(minIntensity)).append("</td></tr>").append("\n")

				ret.append("<tr><td>Max. Intensity</td><td class=\"numberCell\">")
				ret.append(df.format(maxIntensity)).append("</td></tr>").append("\n")

				ret.append("<tr><td>Min. Scan Angle</td><td class=\"numberCell\">")
				ret.append(df.format(minScanAngle)).append("</td></tr>").append("\n")

				ret.append("<tr><td>Max. Scan Angle</td><td class=\"numberCell\">")
				ret.append(df.format(maxScanAngle)).append("</td></tr>").append("\n")

				ret.append("</table>")


				ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
				ret.append("<caption>Point Return Table</caption>")
				
				ret.append("<tr><th>Return Value</th><th>Number</th></th><th>Percentage</th></tr>")

				long[] pointsByRet = las.getNumPointsByReturn()
				for (int n = 0; n < pointsByRet.length; n++) {
					if (pointsByRet[n] > 0) {
						ret.append("<tr><td>").append(n + 1)
						ret.append("</td><td class=\"numberCell\">")
						ret.append(df.format(pointsByRet[n]))
						ret.append("<td class=\"numberCell\">")
						ret.append(df2.format(pointsByRet[n] * 1.0d / totalPoints))
						ret.append("</td></tr>").append("\n")
					}
				}
				
				ret.append("</table>")


				ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
				ret.append("<caption>Point Position Table</caption>")
				
				ret.append("<tr><th>Return Position</th><th>Number</th></th><th>Percentage</th></tr>")

				ret.append("<tr><td>Only</td><td class=\"numberCell\">")
				ret.append(df.format(onlyRet))
				ret.append("<td class=\"numberCell\">")
				ret.append(df2.format(1f * onlyRet / totalPoints))
				ret.append("</td></tr>").append("\n")

				ret.append("<tr><td>First</td><td class=\"numberCell\">")
				ret.append(df.format(numFirstRet))
				ret.append("<td class=\"numberCell\">")
				ret.append(df2.format(1f * numFirstRet / totalPoints))
				ret.append("</td></tr>").append("\n")

				ret.append("<tr><td>Intermediate</td><td class=\"numberCell\">")
				ret.append(df.format(intermediateRet))
				ret.append("<td class=\"numberCell\">")
				ret.append(df2.format(1f * intermediateRet / totalPoints))
				ret.append("</td></tr>").append("\n")

				ret.append("<tr><td>Last</td><td class=\"numberCell\">")
				ret.append(df.format(numLastRet))
				ret.append("<td class=\"numberCell\">")
				ret.append(df2.format(1f * numLastRet / totalPoints))
				ret.append("</td></tr>").append("\n")

				ret.append("</table>")


				ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
				ret.append("<caption>Point Classification Table</caption>")
				
				ret.append("<tr><th>Classification</th><th>Number</th></th><th>Percentage</th></tr>")

				String[] classificationVals = ["Never classified", 
				 "Unclassified", "Ground", "Low Vegetation",
				 "Medium Vegetation", "High Vegetation", "Building",
				 "Low Point (noise)", "Model Key-point", "Water"]
				
				for (int n = 0; n < classHisto.length; n++) {
					if (classHisto[n] > 0) {
						if (n < classificationVals.length) {
							ret.append("<tr><td>").append(df.format(n + 1)).append(" ").append(classificationVals[n])
							
						} else {
							ret.append("<tr><td>").append(df.format(n + 1))
							
						}
						ret.append("</td><td class=\"numberCell\">")
						ret.append(df.format(classHisto[n])).append("</td>")
						ret.append("<td class=\"numberCell\">")
						ret.append(df2.format(1f * classHisto[n] / totalPoints))
						ret.append("</td></tr>").append("\n")
					}
				}
				ret.append("</table><br>")

				ret.append("<p><b>Number of Variable Length Records:</b>").append(" &nbsp ").append(las.getNumVLR()).append("</p>")

				ArrayList<VariableLengthRecord> vlrecs = las.getVariableLengthRecords()
				int r = 1
				for (VariableLengthRecord vlr : vlrecs) {
					ret.append("<p><b><i>Record ")
					ret.append(r).append("</i></b> &nbsp ")
					ret.append("<br>User ID: &nbsp ")
					ret.append(vlr.getUserID())
					ret.append("<br>Record ID: &nbsp ")
					ret.append(vlr.getRecordID())
					ret.append("<br>Description: &nbsp ")
					ret.append(vlr.getDescription())
					ret.append("<br>Length: &nbsp ")
					ret.append(vlr.getRecordLengthAfterHeader())
					ret.append("<br>Data: &nbsp ")
					ret.append(vlr.getFormatedData())
					ret.append("</p>")
					r++
				}

				ret.append("</body></html>")
				pluginHost.returnData(ret.toString());
				i++
			}
	
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
        }
	}
	
	@Override
    public void actionPerformed(ActionEvent event) {
    	if (event.getActionCommand().equals("ok")) {
    		final def args = sd.collectParameters()
			sd.dispose()
			final Runnable r = new Runnable() {
            	@Override
            	public void run() {
                	execute(args)
            	}
        	}
        	final Thread t = new Thread(r)
        	t.start()
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new GetLasFileSummary(pluginHost, args, name, descriptiveName)
}
