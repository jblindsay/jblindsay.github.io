/*
 * Copyright (C) 2013 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import java.text.DecimalFormat
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.PriorityQueue
import java.util.Arrays
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import whitebox.structures.BooleanBitArray1D
import java.util.Arrays
import org.apache.commons.math3.distribution.NormalDistribution
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
def name = "EliminateCoincidentPoints"
def descriptiveName = "Eliminate Coincident Points"
def description = "Eliminates the coincident points within a vector points coverage"
def toolboxes = ["VectorTools"]

public class EliminateCoincidentPoints implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
    private String basisFunctionType = ""
	
    public EliminateCoincidentPoints(WhiteboxPluginHost pluginHost, 
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
        	sd.addDialogFile("Input file", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Vector Points File:", "save", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogDataInput("Distance tolerance", "Tolerance Distance:", "0.0", true, false)
			
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
        	int progress, oldProgress, i, j
        	double x, y
        	ArrayList<Double> xList = new ArrayList<>()
            ArrayList<Double> yList = new ArrayList<>()
            
			if (args.length != 3) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
        	// read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
            double distThreshold = Double.parseDouble(args[2])
            distThreshold = distThreshold * distThreshold
            
			ShapeFile input = new ShapeFile(inputFile)
			if (input.getShapeType() != ShapeType.POINT) {
            	pluginHost.showFeedback("The input file must be of a POINT ShapeType.")
            	return
            }

			AttributeTable table = input.getAttributeTable()
			
			DBFField[] fields = table.getAllFields()

            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields)

			int numFeatures = input.getNumberOfRecords()
			double[][] point
			for (ShapeFileRecord record : input.records) {
				point = record.getGeometry().getPoints()
				for (int p = 0; p < point.length; p++) {
					xList.add(point[p][0])
					yList.add(point[p][1])
				}
				i++
                progress = (int)(100f * i / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Reading Points:", progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
			
			int numPoints = xList.size()
			BooleanBitArray1D isDuplicate = new BooleanBitArray1D(numPoints)
			double dist, xN, yN
			for (i = 0; i < (numPoints - 1); i++) {
				if (!isDuplicate.getValue(i)) {
 					x = xList.get(i)
					y = yList.get(i)
					for (j = i + 1; j < numPoints; j++) {
						xN = xList.get(j)
						dist = (xN - x) * (xN - x)
						if (dist < distThreshold) {
							yN = yList.get(j)
							dist += (yN - y) * (yN - y)
							if (dist < distThreshold) {
								isDuplicate.setValue(j, true)
							}
						}
					}
				}
				progress = (int)(100f * i / numPoints)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Searching For Duplicates:", progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			whitebox.geospatialfiles.shapefile.Point wbGeometry
			oldProgress = -1
			int numNonduplicate = 0
			for (i = 0; i < numPoints; i++) {
				if (!isDuplicate.getValue(i)) {
					numNonduplicate++
					x = xList.get(i)
					y = yList.get(i)
					wbGeometry = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
                    Object[] rowData = table.getRecord(i)
                    output.addRecord(wbGeometry, rowData);
				}
				progress = (int)(100f * i / numPoints)
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Writing New File:", progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			output.write()
			
			// display the output image
			pluginHost.returnData(outputFile)

			// output the number of duplicates found
			int numduplicates = numPoints - numNonduplicate
			String str = "There were " + numduplicates + " coincident points located in the input file."
			pluginHost.returnData(str)
			
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
    def f = new EliminateCoincidentPoints(pluginHost, args, name, descriptiveName)
}
