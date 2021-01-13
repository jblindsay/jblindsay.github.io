/*
 * Copyright (C) 2015 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
 
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.concurrent.Future;
import java.util.concurrent.*;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.geospatialfiles.LASReader;
import whitebox.geospatialfiles.LASReader.PointRecord;
import whitebox.ui.plugin_dialog.ScriptDialog;
import whitebox.utilities.StringUtilities;
import whitebox.geospatialfiles.LASReader.VariableLengthRecord;
import whitebox.geospatialfiles.ShapeFile;
import whitebox.geospatialfiles.shapefile.*;
import whitebox.geospatialfiles.shapefile.attributes.*;
import whitebox.geospatialfiles.VectorLayerInfo;
import whitebox.utilities.Topology;
import groovy.transform.CompileStatic;
import whitebox.structures.BoundingBox;
import com.vividsolutions.jts.geom.PrecisionModel;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;
import com.vividsolutions.jts.triangulate.VoronoiDiagramBuilder;
import whitebox.structures.BooleanBitArray1D;

/*
 * This tool can be used to identify points within a point cloud
 * contained within a LAS file that correspond with the ground 
 * surface. The points are then output into a MultiPoint shapefile.
 */
def name = "LAS2MultiPointZShape"
def descriptiveName = "Convert LAS to MultiPointZ Shapefile"
def description = "Isolates points associated with the ground surface in a LiDAR point cloud."
def toolboxes = ["LidarTools"]

public class LAS2MultiPointZShape {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public LAS2MultiPointZShape(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.descriptiveName = descriptiveName
			
		if (args.length > 0) {
			execute(args)
		} else {
			// create an ActionListener to handle the return from the dialog
            def ac = new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    if (event.getActionCommand().equals("ok")) {
                        args = sd.collectParameters()
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
        	};
			// Create a dialog for this tool to collect user-specified
			// tool parameters.
			sd = new ScriptDialog(pluginHost, descriptiveName, ac)	
		
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
			sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
            sd.addDialogFile("Output file", "Output Vector File:", "close", "Vector Files (*.shp), SHP", true, false)
            
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
	  		// make sure there are the appropriate number of arguments
		  	if (args.length != 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}

			// declare some variables
			int numPoints;
			long numPointsLong;
			int progress, oldProgress = -1;
			PointRecord point;
	    	
			// Read the input parameters
			String inputFile = args[0];
            String outputFile = args[1];
			
	    	// Create the LAS object
            LASReader las = new LASReader(inputFile);
			numPointsLong = las.getNumPointRecords();
			
			// Create the output file
			DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("ID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);
            
			ShapeFile output = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields);

			if (numPointsLong < (long)2147483647) {
				numPoints = (int)numPointsLong;
				
				double[][] xyData = new double[numPoints][2]
	            double[] zData = new double[numPoints]
				double[] mData = new double[numPoints]
				for (int a = 0; a < numPoints; a++) {
					point = las.getPointRecord(a);
					xyData[a][0] = point.getX();
					xyData[a][1] = point.getY();
					zData[a] = point.getZ();
					mData[a] = point.getIntensity();
					progress = (int) (100f * (a + 1) / numPoints);
	                if (progress != oldProgress) {
	                    oldProgress = progress;
	                    pluginHost.updateProgress("Saving output:", progress);
	                    if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
	                }
	            }
	
	            MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
	            Object[] rowData = new Object[1];
	            rowData[0] = new Double(1.0);
	            output.addRecord(wbPoint, rowData);

			} else {
				pluginHost.showFeedback("This tool doesn't currently support LAS files with greater than \n2,147,483,647 points. Please contact the plugin author to request \nsupport.")
				return
//				/* A multipoint record can fit a maximum of 
//				 *  2,147,483,647 points so we will have to use  
//				 *  multiple records.
//				 */
//				long start = 0;
//				long currentVal = 0;
//				numPoints = (int)2147483647;
//				int recNum = 1;
//				boolean flag = true;
//				while (flag) {
//					currentVal += 2147483647L;
//					if (currentVal > numPointsLong) {
//						flag = false;
//						currentVal -= 2147483647L;
//						numPoints -= (int)(numPointsLong - currentVal);
//						currentVal = numPointsLong;
//					}
//					
//					double[][] xyData = new double[numPoints][2]
//		            double[] zData = new double[numPoints]
//					double[] mData = new double[numPoints]
//					int b = 0;
//					for (long a = start; a < currentVal; a++) {
//						point = las.getPointRecord(a);
//						xyData[b][0] = point.getX();
//						xyData[b][1] = point.getY();
//						zData[b] = point.getZ();
//						mData[b] = point.getIntensity();
//						progress = (int) (100f * (b + 1) / numPoints);
//						b++;
//		                if (progress != oldProgress) {
//		                    oldProgress = progress;
//		                    pluginHost.updateProgress("Saving output:", progress);
//		                    if (pluginHost.isRequestForOperationCancelSet()) {
//		                        pluginHost.showFeedback("Operation cancelled")
//								return
//		                    }
//		                }
//		            }
//		
//		            MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
//		            Object[] rowData = new Object[1];
//		            rowData[0] = new Double(recNum);
//		            output.addRecord(wbPoint, rowData);
//					recNum++;
//					start = currentVal;
//				}
			}
			
            output.write()


            // display the output image
			String paletteDirectory = pluginHost.getResourcesDirectory() + "palettes" + File.separator;
			VectorLayerInfo vli = new VectorLayerInfo(outputFile, paletteDirectory, 255i, -1);
			vli.setPaletteFile(paletteDirectory + "spectrum.pal");
			vli.setFilledWithOneColour(false);
			vli.setFillAttribute("Feature Z Value");
			vli.setPaletteScaled(true);
			vli.setDisplayMinValue(las.getMinZ());
			vli.setDisplayMaxValue(las.getMaxZ());
			vli.setMarkerSize(2.5f);
			vli.setRecordsColourData();
			pluginHost.returnData(vli);
      
	  	} catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
	    } catch (Exception e) {
	        pluginHost.showFeedback("$e"); //"An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress("", 0)
        }
	}
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myClass = new LAS2MultiPointZShape(pluginHost, args, name, descriptiveName)
}
