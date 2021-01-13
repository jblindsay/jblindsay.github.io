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
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import java.util.Collections
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.Color;
import javax.swing.JPanel;
import java.awt.Dimension
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.Month;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleAnchor;

import whitebox.plugins.PluginInfo
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "Profile"
def descriptiveName = "Profile"
def description = "Creates a topographic profile"
def toolboxes = ["TerrainAnalysis"]

// The following lines are necessary for the script 
// to be recognized as a menu extension.
parentMenu = "Tools"
menuLabel = "Profile"

public class Profile implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public Profile(WhiteboxPluginHost pluginHost, 
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
        	DialogFile dfIn1 = sd.addDialogFile("Input surface file", "Input Raster Surface File:", "open", "Raster Files (*.dep), DEP", true, false)
            DialogFile dfIn2 = sd.addDialogFile("Input vector lines file; must be of a PolyLine ShapeType.", "Input Vector Lines File:", "open", "Vector Files (*.shp), SHP", true, false)

			def btn = sd.addDialogButton("Or create a new line vector now...", "centre")
			btn.addActionListener(new ActionListener() {
 	            public void actionPerformed(ActionEvent e) {
 	            	// find an appropriate file name for it
					def outputFile = pluginHost.getWorkingDirectory() + "ProfileLine.shp";
					def file = new File(outputFile);
					if (file.exists()) {
						for (int i = 1; i < 101; i++) {
							outputFile = pluginHost.getWorkingDirectory() + "ProfileLine${i}.shp";
							file = new File(outputFile);
							if (!file.exists()) {
								break;
							}
						}
					}
					DBFField[] fields = new DBFField[1];
		            
		            fields[0] = new DBFField();
		            fields[0].setName("FID");
		            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
		            fields[0].setFieldLength(10);
		            fields[0].setDecimalCount(0);
		            
		            ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYLINE, fields);
		            output.write();
		            
		            pluginHost.returnData(outputFile);
		            
		            pluginHost.editVector();

		            pluginHost.showFeedback("Press the Digitize New Feature icon on the toolbar \n" +
		            "to add a line. Then toggle the Edit Vector icon when you \n" + 
		            "are done digitizing");
 	            }
 	        });  

 	        DialogFile dfOut = sd.addDialogFile("Output text file (optional)", "Output Text File (Blank for none):", "save", "Text Files (*.csv), CSV", true, true)
            
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
        	int progress, oldProgress
            int row, col, i
	        double x, y, z, nodata
	        int cols, rows 
	        double minZ = Double.POSITIVE_INFINITY
	        double maxZ = Double.NEGATIVE_INFINITY
	        
	        ArrayList<com.vividsolutions.jts.geom.Geometry> pointList = new ArrayList<>()
			double[][] points
            double[] zArray
			Object[] recData
			int startingPointInPart, endingPointInPart
			ArrayList<Integer> edgeList = new ArrayList<>()
			double x1, y1, x2, y2, xPrime
        	boolean foundIntersection
        	BoundingBox box 
        	
            if (args.length < 2) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputSurfaceFile = args[0]
            String inputLines = args[1]
            boolean textOutput = false
            String outputFileName = ""
            if (args.length == 3) {
            	if (!args[2].isEmpty() && 
            	     !args[2].toLowerCase().equals("not specified")) {
	            	textOutput = true
	            	outputFileName = args[2]
            	}
            }
			
			WhiteboxRaster surface = new WhiteboxRaster(inputSurfaceFile, "r")
			rows = surface.getNumberRows()
			cols = surface.getNumberColumns()
			nodata = surface.getNoDataValue()
			String domainLabel = "Distance"
			String xyUnits = surface.getXYUnits()
			if (xyUnits == null || !xyUnits.toLowerCase().equals("not specified")) {
				domainLabel = domainLabel + " (${xyUnits})"
			}
			String rangeLabel = "Value"
			String zUnits = surface.getZUnits()
			if (zUnits == null || !zUnits.toLowerCase().equals("not specified")) {
				rangeLabel = rangeLabel + " (${zUnits})"
			}
			

			ShapeFile input = new ShapeFile(inputLines)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POLYLINE) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYLINE ShapeType.")
            	return
            }

            XYSeriesCollection xyCollection = new XYSeriesCollection();

			StringBuilder sb = new StringBuilder()
        	
			int numFeatures = input.getNumberOfRecords()
            int count = 0
            int numSeries = 0
            oldProgress = -1
			for (ShapeFileRecord record : input.records) {
				int recNum = record.getRecordNumber()
                PolyLine pl = (PolyLine)(record.getGeometry())
				points = record.getGeometry().getPoints()
				int numPoints = points.length;
				int[] partData = record.getGeometry().getParts()
				int numParts = partData.length
				
				for (int part = 0; part < numParts; part++) {
					sb.append("X,Y,Distance,Z\n")
					XYSeries series
					if (numParts > 1) {
						series = new XYSeries("Profile ${recNum}-${part}")
					} else {
						series = new XYSeries("Profile ${recNum}")
					}
					numSeries++
					
//                	ArrayList<Double> zList = new ArrayList<>()
//					ArrayList<Double> distList = new ArrayList<>()
					
                    startingPointInPart = partData[part];
                    if (part < numParts - 1) {
                        endingPointInPart = partData[part + 1] - 1;
                    } else {
                        endingPointInPart = numPoints - 1;
                    }

					row = surface.getRowFromYCoordinate(points[0][1])
                    col = surface.getColumnFromXCoordinate(points[0][0])
                    z = surface.getValue(row, col)
//                    zList.add(z)
                    double dist = 0
//					distList.add(dist)
                    
					double dX, dY
					int endingRow, endingCol, startingRow, startingCol
					double pathDist, distStep
					int numSteps
					for (i = startingPointInPart; i < endingPointInPart; i++) {
                    	startingRow = surface.getRowFromYCoordinate(points[i][1])
                    	startingCol = surface.getColumnFromXCoordinate(points[i][0])
                    	
                    	endingRow = surface.getRowFromYCoordinate(points[i + 1][1])
                    	endingCol = surface.getColumnFromXCoordinate(points[i + 1][0])
                    	
                    	dX = endingCol - startingCol
                    	dY = endingRow - startingRow

                    	pathDist = Math.sqrt(dX * dX + dY * dY)
						numSteps = (int)(Math.ceil(pathDist))
						
						dX = dX / pathDist
						dY = dY / pathDist

						distStep = Math.sqrt((points[i][0] - points[i + 1][0]) * (points[i][0] - points[i + 1][0]) + (points[i][1] - points[i + 1][1]) * (points[i][1] - points[i + 1][1])) / pathDist

						if (numSteps > 0) {
							for (int j in 1..numSteps) {
								col = (int)(startingCol + j * dX)
								row = (int)(startingRow + j * dY)
								z = surface.getValue(row, col)
								dist += distStep
								if (z != nodata) {
									if (z < minZ) { minZ = z }
									if (z > maxZ) { maxZ = z }
//									zList.add(z)
//									distList.add(dist)
									series.add(dist, z)
									if (textOutput) {
										x = surface.getXCoordinateFromColumn(col)
										y = surface.getYCoordinateFromRow(row)
										sb.append("${x},${y},${dist},${z}\n")
									}
								}
							}
						}

                    }

					xyCollection.addSeries(series);
				}
				
				count++
                progress = (int)(100f * count / numFeatures)
            	if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
            }

			JFreeChart chart = ChartFactory.createXYLineChart(
            "Profile",
            domainLabel,
            rangeLabel,
            xyCollection,
            PlotOrientation.VERTICAL,  // Plot Orientation
            true,                      // Show Legend
            true,                      // Use tooltips
            false                      // Configure chart to generate URLs?
            );

			XYPlot plot = (XYPlot)chart.getPlot();
			plot.setBackgroundPaint(Color.WHITE); 
			plot.setDomainGridlinePaint(Color.BLACK); 
			plot.setRangeGridlinePaint(Color.BLACK); 

			if (numSeries == 1 || numSeries > 8) {
				chart.getLegend().setVisible(false);
			}

			double zRange = maxZ - minZ
			minZ = minZ - zRange / 10.0
			maxZ = maxZ + zRange / 10.0
	        NumberAxis range = (NumberAxis) plot.getRangeAxis();
	        range.setRange(minZ, maxZ);

            ChartPanel chartPanel = new ChartPanel(chart)
			chartPanel.setPreferredSize(new Dimension(700, 300))
			
        	pluginHost.returnData(chartPanel)

        	if (textOutput) {
        		def outFile = new File(outputFileName)
        		if (outFile.exists()) { outFile.delete() }
        		outFile.text = sb.toString()
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

//    // Return true if val is between theshold1 and theshold2.
//    private static boolean isBetween(double val, double threshold1, double threshold2) {
//        if (val == threshold1 || val == threshold2) {
//            return true;
//        }
//        return threshold2 > threshold1 ? val > threshold1 && val < threshold2 : val > threshold2 && val < threshold1;
//    }
//
//    private class RecordInfo implements Comparable<RecordInfo> {
//
//        public double maxY;
//        public int recNumber;
//        
//        public RecordInfo(double maxY, int recNumber) {
//            this.maxY = maxY;
//            this.recNumber = recNumber;
//        }
//
//        @Override
//        public int compareTo(RecordInfo other) {
//            final int BEFORE = -1;
//            final int EQUAL = 0;
//            final int AFTER = 1;
//
//            if (this.maxY < other.maxY) {
//                return BEFORE;
//            } else if (this.maxY > other.maxY) {
//                return AFTER;
//            }
//
//            if (this.recNumber < other.recNumber) {
//                return BEFORE;
//            } else if (this.recNumber > other.recNumber) {
//                return AFTER;
//            }
//
//            return EQUAL;
//        }
//    }
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new Profile(pluginHost, args, name, descriptiveName)
}
