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
 
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger
import java.util.Date
import java.util.ArrayList
import whitebox.ui.plugin_dialog.*
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.ShapeFileInfo
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.structures.BoundingBox
import whitebox.structures.KdTree
import whitebox.structures.BooleanBitArray1D
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "MinimumInterpolation"
def descriptiveName = "Minimum Interpolation"
def description = "Interpolates a vector point file using a minimum value criterion."
def toolboxes = ["Interpolation"]

public class MinimumInterpolation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
//	private KdTree<Double> pointsTree
		
	private AtomicInteger numSolvedRows = new AtomicInteger(0)
	
	public MinimumInterpolation(WhiteboxPluginHost pluginHost, 
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
			DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and height field.", "Input Height Field:", false)
            DialogCheckBox dcb = sd.addDialogCheckBox("Use z-values", "Use z-values", false)
            dcb.setVisible(false)
            sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, true)
            sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)
            sd.addDialogDataInput("Max Search Distance", "Max Search Distance", "", true, false)
			
			def listener = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = dfs.getValue()
            		if (value != null && !value.isEmpty()) {
            			value = value.trim()
            			String[] strArray = dfs.getValue().split(";")
            			String fileName = strArray[0]
            			File file = new File(fileName)
            			if (file.exists()) {
	            			ShapeFileInfo shapefile = new ShapeFileInfo(fileName)
		            		if (shapefile.getShapeType().getDimension() == ShapeTypeDimension.Z) {
		            			dcb.setVisible(true)
		            		} else {
		            			dcb.setVisible(false)
		            		}
		            	} else {
		            		if (dcb.isVisible()) {
		            			dcb.setVisible(false)
		            		}
		            	}
            		}
            	} 
            } as PropertyChangeListener
            dfs.addPropertyChangeListener(listener)

            
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
		long start = System.currentTimeMillis()  
		  try {
		  	int progress, oldProgress, rows, cols, row, col
	    	double x, y, z
	    	double north, south, east, west
	    	double cellSize = -1.0
	    	double nodata = -32768
	    	ArrayList<Double> xList = new ArrayList<>()
	        ArrayList<Double> yList = new ArrayList<>()
	        ArrayList<Double> zList = new ArrayList<>()
		        
		  	if (args.length != 6) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters
			String[] inputData = args[0].split(";")
	        String inputFile = inputData[0]
			boolean useZValues = Boolean.parseBoolean(args[1])
	        String outputFile = args[2]
			if (!args[3].toLowerCase().contains("not specified")) {
	            cellSize = Double.parseDouble(args[3]);
	        }
	        String baseFileHeader = args[4]
	        if (baseFileHeader == null || baseFileHeader.isEmpty()) {
	        	baseFileHeader = "not specified"
	        }
	        double maxDist = Double.parseDouble(args[5]);
			
	        AttributeTable table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
			ShapeFile input = new ShapeFile(inputFile)
			if (cellSize < 0 && baseFileHeader.toLowerCase().equals("not specified")) {
				cellSize = Math.max((input.getxMax() - input.getxMin()) / 998, (input.getyMax() - input.getyMin()) / 998)
			}
			ShapeType shapeType = input.getShapeType()
	        if (shapeType.getDimension() != ShapeTypeDimension.Z && useZValues) {
	        	useZValues = false
	        }
			int heightField = -1
	        if (inputData.length == 2 && !inputData[1].trim().isEmpty()) {
	        	String heightFieldName = inputData[1].trim()
	        	String[] fieldNames = input.getAttributeTableFields()
	        	for (int i = 0; i < fieldNames.length; i++) {
	        		if (fieldNames[i].trim().equals(heightFieldName)) {
	        			heightField = i
	        			break
	        		}
	        	}
	        } else if (!useZValues) {
	        	pluginHost.showFeedback("A field within the input file's attribute table must be selected to assign point heights.")
	        	return
	        }

			double[][] point
			Object[] recData
			int i = 0
			int numFeatures = input.getNumberOfRecords()
			oldProgress = -1
			if (!useZValues) {
				for (ShapeFileRecord record : input.records) {
					recData = table.getRecord(i)
					z = (Double)(recData[heightField])
					point = record.getGeometry().getPoints()
					for (int p = 0; p < point.length; p++) {
						xList.add(point[p][0])
						yList.add(point[p][1])
						zList.add(z)
					}
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}
				}
			} else {
				for (ShapeFileRecord record : input.records) {
					if (shapeType.getBaseType() == ShapeType.POINT) {
						PointZ ptz = (PointZ)(record.getGeometry())
	            		xList.add(ptz.getX())
						yList.add(ptz.getY())
						zList.add(ptz.getZ())
					} else if (shapeType.getBaseType() == ShapeType.MULTIPOINT) {
						MultiPointZ plz = (MultiPointZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = plz.getzArray()
						for (int p = 0; p < point.length; p++) {
							xList.add(point[p][0])
							yList.add(point[p][1])
							zList.add(zArray[p])
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYLINE) {
						PolyLineZ plz = (PolyLineZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = plz.getzArray()
						for (int p = 0; p < point.length; p++) {
							xList.add(point[p][0])
							yList.add(point[p][1])
							zList.add(zArray[p])
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYGON) {
						PolygonZ pz = (PolygonZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = pz.getzArray()
						for (int p = 0; p < point.length; p++) {
							xList.add(point[p][0])
							yList.add(point[p][1])
							zList.add(zArray[p])
						}
					}
					
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            		// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
	            	}	            	
				}
			}
			
			int numSamples = zList.size()
			//final KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numSamples));
			//pointsTree = new KdTree.SqrEuclid<Integer>(2, new Integer(numSamples));
			KdTree<Double> pointsTree = new KdTree.SqrEuclid<Double>(2, new Integer(numSamples));
			for (i = 0; i < numSamples; i++) {
				double[] entry = new double[2]
				entry[0] = xList.get(i)
				entry[1] = yList.get(i);
				pointsTree.addPoint(entry, zList.get(i));
			}
			
			xList.clear()
			yList.clear()
			zList.clear()
			
			north = input.getyMax() + cellSize / 2.0;
	        south = input.getyMin() - cellSize / 2.0;
	        east = input.getxMax() + cellSize / 2.0;
	        west = input.getxMin() - cellSize / 2.0;
	            
			// initialize the output raster
	        WhiteboxRaster image;
	        if ((cellSize > 0) || ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified")))) {
	            if ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified"))) {
	                cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0,
	                        (input.getxMax() - input.getxMin()) / 500.0);
	            }
	            rows = (int) (Math.ceil((north - south) / cellSize));
	            cols = (int) (Math.ceil((east - west) / cellSize));
	
	            // update west and south
	            east = west + cols * cellSize;
	            south = north - rows * cellSize;
	
	            image = new WhiteboxRaster(outputFile, north, south, east, west,
	                rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
	                WhiteboxRasterBase.DataType.FLOAT, nodata, nodata);
	        } else {
	            image = new WhiteboxRaster(outputFile, "rw",
	                    baseFileHeader, WhiteboxRasterBase.DataType.FLOAT, nodata);
	            rows = image.getNumberRows()
	            cols = image.getNumberColumns()
	        }
	        
	        image.setPreferredPalette("spectrum.pal")

			double[] entry;
			List<KdTree.Entry<Double>> results;
			double value;
            double minVal;
	        for (row = 0; row < rows; row++) {
	        	y = image.getYCoordinateFromRow(row);
                for (col = 0; col < cols; col++) {
                    x = image.getXCoordinateFromColumn(col);
                    entry = [x, y];
                    results = pointsTree.neighborsWithinRange(entry, maxDist);
                    if (results.size() > 1) {
						// now find the minimum value
                        minVal = Double.POSITIVE_INFINITY
                        for (i = 0; i < results.size(); i++) {
                        	z = (double)results.get(i).value
                        	if (z < minVal) { 
                            	minVal = z
                        	}
                        }
	                     
                        image.setValue(row, col, minVal)
                    } else if (results.size() == 1) {
                    	image.setValue(row, col, (double)results.get(0).value)
                    } else {
                    	image.setValue(row, col, nodata)
                    }
                }
                progress = (int)(100f * row / (rows - 1))
                if (progress != oldProgress) {
                	oldProgress = progress
                	pluginHost.updateProgress(progress)
                	
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
                }
	        }
			
		    image.addMetadataEntry("Created by the " + descriptiveName + " tool.");
	        image.addMetadataEntry("Created on " + new Date());
	        image.addMetadataEntry("Max. Search Distance:\t${maxDist}");
	        image.addMetadataEntry("Resolution:\t${cellSize}");
	        image.close();
		    
			pluginHost.returnData(outputFile);
	
			long end = System.currentTimeMillis()  
			double duration = (end - start) / 1000.0
			pluginHost.showFeedback("Interpolation completed in " + duration + " seconds.")
		
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("Progress", 0)
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
	def myClass = new MinimumInterpolation(pluginHost, args, name, descriptiveName)
}
