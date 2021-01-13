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
import java.util.Arrays
import java.util.stream.IntStream
import java.util.concurrent.atomic.AtomicInteger
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours
import whitebox.structures.KdTree
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.LASReader.VariableLengthRecord
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import groovy.transform.CompileStatic
import whitebox.structures.BoundingBox

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "FilterLasScanAngles"
def descriptiveName = "Filter LAS Scan Angles"
def description = "Removes points in a LAS file with scan angles greater than a threshold."
def toolboxes = ["LidarTools"]

public class FilterLasScanAngles implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName

	private AtomicInteger numSolved = new AtomicInteger(0)
	private AtomicInteger curProgress = new AtomicInteger(-1)
	private int numFiles
	
	public FilterLasScanAngles(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogDataInput("Enter the maximum scan angle threshold", "Scan Angle Threshold (deg.):", "", true, false)
            
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
		  	if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			final String inputFileString = args[0]
			double maxScanAngle = Double.parseDouble(args[1])

            String[] inputFiles = inputFileString.split(";")

			// check for empty entries in the inputFiles array
			numFiles = inputFiles.length
			IntStream.range(0, numFiles).parallel().forEach{ 
				performFilter(inputFiles[it], maxScanAngle) 
			};

            if (numFiles > 1) {
				pluginHost.showFeedback("Operation complete.\n${numFiles} files were created.")
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

	private performFilter(String lasFile, double threshold) {
		double scanAngle
		int progress, oldProgress = -1

		if (lasFile == null || lasFile.isEmpty()) {
			return 
		}
		
		if (pluginHost.isRequestForOperationCancelSet()) {
			return
        }
		
		LASReader las = new LASReader(lasFile)
		int totalPoints = (int)las.getNumPointRecords()

		// first scan through and count how many points have 
		// scan angles less than the threshold
		int numPoints = 0
		for (int p = 0; p < totalPoints; p++) {
			PointRecord point = las.getPointRecord(p)
			scanAngle = Math.abs(point.getScanAngle())
			if (scanAngle < threshold) { numPoints++ }
			if (numFiles == 1) {
				progress = (int) (100f * (p + 1) / totalPoints);
                if (progress != oldProgress) {
                    oldProgress = progress;
                    pluginHost.updateProgress("Reading Scan Angles:", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
                }
			}
		}

		// output the points
		DBFField[] fields = new DBFField[1];

        fields[0] = new DBFField();
        fields[0].setName("FID");
        fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
        fields[0].setFieldLength(10);
        fields[0].setDecimalCount(0);

		String outputFile = lasFile.replace(".las", ".shp");

        // see if the output files already exist, and if so, delete them.
        File outFile = new File(outputFile);
        if (outFile.exists()) {
            outFile.delete();
            (new File(outputFile.replace(".shp", ".dbf"))).delete();
            (new File(outputFile.replace(".shp", ".shx"))).delete();
        }
        
        ShapeFile output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);
		double[][] xyData = new double[numPoints][2]
		double[] zData = new double[numPoints]
		double[] mData = new double[numPoints]

		// scan it a second time to create the output shapefile
		int index = 0
		for (int p = 0; p < totalPoints; p++) {
			PointRecord point = las.getPointRecord(p)
			scanAngle = Math.abs(point.getScanAngle())
			if (scanAngle < threshold) { 
				xyData[index][0] = point.getX()
				xyData[index][1] = point.getY()
				zData[index] = point.getZ()
				mData[index] = point.getIntensity()
				index++
			}
			if (numFiles == 1) {
				progress = (int) (100f * (p + 1) / totalPoints);
                if (progress != oldProgress) {
                    oldProgress = progress;
                    pluginHost.updateProgress("Filtering Points:", progress);
                    if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
                }
			}
		}

		if (numFiles == 1) {
			pluginHost.updateProgress("Saving data:", 0);
		}
		MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
		Object[] rowData = new Object[1]
        rowData[0] = new Double(1)
        output.addRecord(wbPoint, rowData)
        output.write()

		if (numFiles == 1) {
			pluginHost.returnData(outputFile);
		} else {
			int solved = numSolved.incrementAndGet()
			progress = (int) (100f * solved / numFiles)
			if (progress > curProgress.get()) {
				curProgress.incrementAndGet()
				pluginHost.updateProgress("$solved of $numFiles filtered:", progress)
			}
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
	def myClass = new FilterLasScanAngles(pluginHost, args, name, descriptiveName)
}
