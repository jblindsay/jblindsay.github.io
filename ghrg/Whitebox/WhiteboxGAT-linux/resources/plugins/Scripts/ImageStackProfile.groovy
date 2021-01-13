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
def name = "ImageStackProfile"
def descriptiveName = "Image Stack Profile"
def description = "Creates an image stack profile (signature)"
def toolboxes = ["ImageProc"]

public class ImageStackProfile implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ImageStackProfile(WhiteboxPluginHost pluginHost, 
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
        	sd.addDialogMultiFile("Select some input raster files", "Input Raster Files:", "Raster Files (*.dep), DEP")
			DialogFile dfIn2 = sd.addDialogFile("Input vector points file; must be of a POINT ShapeType.", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)

			def btn = sd.addDialogButton("Or create a new point vector now...", "centre")
			btn.addActionListener(new ActionListener() {
 	            public void actionPerformed(ActionEvent e) {
 	            	// find an appropriate file name for it
					def outputFile = pluginHost.getWorkingDirectory() + "SignaturePoints.shp";
					def file = new File(outputFile);
					if (file.exists()) {
						for (int i = 1; i < 101; i++) {
							outputFile = pluginHost.getWorkingDirectory() + "SignaturePoints${i}.shp";
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

		            pluginHost.showFeedback("Press the Digitize New Feature icon on the toolbar \n" +
		            "to add a point. Then toggle the Edit Vector icon when you \n" + 
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
            String inputFileString = args[0]
            String inputPoints = args[1]
            boolean textOutput = false
            String outputFileName = ""
            if (args.length == 3) {
            	if (!args[2].isEmpty() && 
            	     !args[2].toLowerCase().equals("not specified")) {
	            	textOutput = true
	            	outputFileName = args[2]
            	}
            }

            String[] inputFiles = inputFileString.split(";")

			int numRasters = inputFiles.length
            WhiteboxRaster[] rasters = new WhiteboxRaster[numRasters]
			i = 0
			for (String inputFile : inputFiles) {
				rasters[i] = new WhiteboxRaster(inputFile, "r")
				i++
			}
			
			String domainLabel = "Image Number"
			
			String rangeLabel = "Value"
			String zUnits = rasters[0].getZUnits()
			if (zUnits == null || !zUnits.toLowerCase().equals("not specified")) {
				rangeLabel = rangeLabel + " (${zUnits})"
			}
			
			ShapeFile input = new ShapeFile(inputPoints)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POINT) {
            	pluginHost.showFeedback("The input shapefile should be of a POINT ShapeType.")
            	return
            }

            XYSeriesCollection xyCollection = new XYSeriesCollection();

			StringBuilder sb = new StringBuilder()
        	
			int numFeatures = input.getNumberOfRecords()
            int count = 0
            int numSeries = 0
            oldProgress = -1
            sb.append("X,Y,Image,Z\n")
			for (ShapeFileRecord record : input.records) {
				int recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				for (int p = 0; p < points.length; p++) {
					XYSeries series = new XYSeries("Profile ${numSeries + 1}")
					numSeries++

					x = points[p][0]
                    y = points[p][1]
                    
                    for (i = 0; i < numRasters; i++) {
                    	col = rasters[i].getColumnFromXCoordinate(x)
                    	row = rasters[i].getRowFromYCoordinate(y)
                    	z = rasters[i].getValue(row, col)
                    	if (z < minZ) { minZ = z }
						if (z > maxZ) { maxZ = z }
                    	series.add((i + 1), z)
                    	if (textOutput) {
							sb.append("${x},${y},${i+1},${z}\n")
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
            "Image Stack Profile",
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
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new ImageStackProfile(pluginHost, args, name, descriptiveName)
}
