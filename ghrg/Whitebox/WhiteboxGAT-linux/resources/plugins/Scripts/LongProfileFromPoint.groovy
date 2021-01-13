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
import java.awt.Dimension
import java.awt.Color
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.shapefile.attributes.DBFField;
import groovy.transform.CompileStatic

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
import org.jfree.chart.renderer.xy.XYItemRenderer
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

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "LongProfileFromPoint"
def descriptiveName = "Long Profile From Point"
def description = "Plots the stream long-profile from one or more vector points."
def toolboxes = ["StreamAnalysis", "HydroTools", "FlowpathTAs"]

public class LongProfileFromPoint {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public LongProfileFromPoint(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
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
			sd.addDialogFile("Input digital elevation model (DEM) file", "Input Digital Elevation Model (DEM):", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input seed points file", "Input Seed Points File:", "open", "Shapefiles (*.shp), SHP", true, false)
            def btn = sd.addDialogButton("Create a new seed point vector now...", "centre")
			btn.addActionListener(new ActionListener() {
 	            public void actionPerformed(ActionEvent e) {
 	            	// find an appropriate file name for it
					def outputFile = pluginHost.getWorkingDirectory() + "Points.shp";
					def file = new File(outputFile);
					if (file.exists()) {
						for (int i = 1; i < 101; i++) {
							outputFile = pluginHost.getWorkingDirectory() + "Points${i}.shp";
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
		            
		            ShapeFile output = new ShapeFile(outputFile, ShapeType.POINT, fields);
					output.write();
		            
		            pluginHost.returnData(outputFile);
		            
		            pluginHost.editVector();

		            pluginHost.showFeedback("Press the 'Digitize New Feature' icon on the toolbar \n" +
		            "to add a point. Then toggle the 'Edit Vector' icon when you \n" + 
		            "are done digitizing");
 	            }
 	        });  
			sd.addDialogFile("Output text file (optional)", "Output Text File (blank for none):", "save", "Text Files (*.csv), CSV", true, true)
            
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
	  		int progress, oldProgress, dir
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			double maxSlope, slope, dist
        	double x, y, z, zN
        	boolean flag
        	double minZ = Double.POSITIVE_INFINITY
	        double maxZ = Double.NEGATIVE_INFINITY
	        double maxDist = Double.NEGATIVE_INFINITY
	        int numPoints = 1;
	        StringBuilder sb = new StringBuilder()
        	
			if (args.length < 2) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String demFile = args[0]
			String inputSeedFile = args[1]
			boolean textOutput = false
            String outputFileName = ""
            if (args.length == 3) {
            	if (args[2] != null && !args[2].isEmpty() && 
            	     !args[2].toLowerCase().equals("not specified")) {
	            	textOutput = true
	            	outputFileName = args[2]
            	}
            }
            
			// read the input image and PP vector files
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double demNodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int cols = dem.getNumberColumns()
			double gridResX = dem.getCellSizeX();
            double gridResY = dem.getCellSizeY();
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = [diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY]

            XYSeriesCollection xyCollection = new XYSeriesCollection();
        	
			def input = new ShapeFile(inputSeedFile)
            ShapeType shapeType = input.getShapeType()
			if (shapeType.getBaseType() != ShapeType.POINT) {
				pluginHost.showFeedback("Error: The seed file must be of a POINT shape type.")
				return;
			}

			int featureNum = 0;
			int row, col;
        	numPoints = input.getNumberOfRecords();
        	double[][] point

        	for (ShapeFileRecord record : input.records) {
        	 	point = record.getGeometry().getPoints()
				col = dem.getColumnFromXCoordinate(point[0][0])
				row = dem.getRowFromYCoordinate(point[0][1])
				z = dem.getValue(row, col);
				
				if (z != demNodata) {
					// it's a valid seed point
					XYSeries series = new XYSeries("Profile ${record.getRecordNumber()}")
					def zList = new ArrayList<Double>()
					def distList = new ArrayList<Double>()
					z = dem.getValue(row, col)
					if (z < minZ) { minZ = z }
					if (z > maxZ) { maxZ = z }
					dist = 0;
					zList.add(z)
					distList.add(dist)
					//series.add(dist, z)
                	// Trace the flowpath from this point
					flag = false;
					while (!flag) {
						z = dem.getValue(row, col)
						
						// What's the flow direction
						dir = -1
						maxSlope = -99999999
						for (int i = 0; i < 8; i++) {
							zN = dem.getValue(row + dY[i], col + dX[i])
							if (zN != demNodata) {
								slope = (z - zN) / gridLengths[i]
								if (slope > maxSlope && slope > 0) {
									maxSlope = slope
									dir = i
								}
							}
						}
						if (dir > -1) {
							row += dY[dir];
							col += dX[dir];
							x = dem.getXCoordinateFromColumn(col);
							y = dem.getYCoordinateFromRow(row);
							z = dem.getValue(row, col);
							if (z < minZ) { minZ = z }
							if (z > maxZ) { maxZ = z }
							dist += gridLengths[dir];
							println("$dist $z")
							if (dist > maxDist) { maxDist = dist }
							//series.add(dist, z)
							zList.add(z)
							distList.add(dist)
						} else {
							flag = true;
							double profileMaxDist = distList.get(distList.size() - 1);
							for (int i = 0; i < zList.size(); i++) {
								series.add(profileMaxDist - distList.get(i), zList.get(i));
							}
							if (textOutput) {
								sb.append("Profile ${featureNum + 1}\n\n")
								sb.append("X,Y,Distance,Z\n")
								for (int i = 0; i < zList.size(); i++) {
									sb.append("${profileMaxDist - distList.get(i)},${zList.get(i)}\n");
								}
							}
							xyCollection.addSeries(series);
						}
					}
				}
				
				progress = (int)(100f * featureNum / (numPoints - 1))
        		if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}
        		
				featureNum++
        	}
			dem.close()
			
			String domainLabel = "Distance To Mouth"
			String xyUnits = dem.getXYUnits()
			if (xyUnits == null || !xyUnits.toLowerCase().equals("not specified")) {
				domainLabel = domainLabel + " (${xyUnits})"
			}
			String rangeLabel = "Elevation"
			String zUnits = dem.getZUnits()
			if (zUnits == null || !zUnits.toLowerCase().equals("not specified")) {
				rangeLabel = rangeLabel + " (${zUnits})"
			}

			JFreeChart chart = ChartFactory.createXYLineChart(
            	"Long Profile",
            	domainLabel,
            	rangeLabel,
	            xyCollection,
	            PlotOrientation.VERTICAL,  // Plot Orientation
	            false,                      // Show Legend
	            true,                      // Use tooltips
	            false                      // Configure chart to generate URLs?
            );

			XYPlot plot = (XYPlot)chart.getPlot();
			plot.setBackgroundPaint(Color.WHITE); 
			plot.setDomainGridlinePaint(Color.BLACK); 
			plot.setRangeGridlinePaint(Color.BLACK); 

			XYItemRenderer renderer = plot.getRenderer();
//			for (int i = 0; i < numPoints; i++) {
//				renderer.setSeriesPaint(i, Color.blue);
//			}

			double zRange = maxZ - minZ
			minZ = minZ - zRange / 10.0
			maxZ = maxZ + zRange / 10.0
	        NumberAxis range = (NumberAxis) plot.getRangeAxis();
	        range.setRange(minZ, maxZ);

	        NumberAxis domain = (NumberAxis) plot.getDomainAxis();
	        domain.setRange(0, maxDist);

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
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def tdf = new LongProfileFromPoint(pluginHost, args, name, descriptiveName)
}
