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
import java.awt.Dimension
import java.awt.Color;
import javax.swing.JPanel;
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
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
def name = "LocalTopographicPositionScaleSignature"
def descriptiveName = "Local Topographic Position Scale Signature"
def description = "Creates a topographic position scale signature."
def toolboxes = ["ElevResiduals"]

public class LocalTopographicPositionScaleSignature implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public LocalTopographicPositionScaleSignature(WhiteboxPluginHost pluginHost, 
		String[] args, def descriptiveName) {
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
			def helpFile = "LocalTopographicPositionScaleSignature"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "LocalTopographicPositionScaleSignature.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			DialogFile dfInput = sd.addDialogFile("Input DEM file", "Input DEM:", "open", "Raster Files (*.dep), DEP", true, false)
			//DialogFile dfIn2 = sd.addDialogFile("Input vector points file; must be of a POINT ShapeType.", "Input Vector Points File:", "open", "Vector Files (*.shp), SHP", true, false)
			DialogFieldSelector dfIn2 = sd.addDialogFieldSelector("Input signature label field (blank for none).", "Input Label Field:", false)
            
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

            def minN = sd.addDialogDataInput("Minimum neighbourhood radius", "Minimum Neighbourhood Radius (cells):", "3", true, false)
            def maxN = sd.addDialogDataInput("Maximum neighbourhood radius", "Maximum Neighbourhood Radius (cells):", "", true, false)
            def stepSize = sd.addDialogDataInput("Step size", "Step Size (cells):", "1", true, false)
			
            //Listener           
            def lsn = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = dfInput.getValue()
            		if (value != null && !value.isEmpty()) { 
            			String fileName = value.trim()
            			File file = new File(fileName)
            			WhiteboxRasterInfo wri = new WhiteboxRasterInfo(fileName)
						int rows = wri.getNumberRows()
						int cols = wri.getNumberColumns()
						int maxSize = (int)(Math.min(rows, cols) * 0.5)
						maxN.setValue(maxSize.toInteger().toString())
            		}
            	} 
            } as PropertyChangeListener
            dfInput.addPropertyChangeListener(lsn)

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
			double minZ = Double.POSITIVE_INFINITY
	        double maxZ = Double.NEGATIVE_INFINITY
	        int row, col, x1, x2, y1, y2
			double outValue, v, s, m, N
			double[][] points
			
			if (args.length != 6) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			double resolution = (image.getCellSizeX() + image.getCellSizeY()) / 2.0
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()

			String[] inputPointsData = args[1].split(";")
            String inputPoints = inputPointsData[0] //args[1]
            String labelFieldName = ""
			if (inputPointsData.length == 2 && !inputPointsData[1].trim().isEmpty()) {
            	labelFieldName = inputPointsData[1].trim()
            }
            boolean textOutput = false
            String outputFileName = ""
            if (!args[2].isEmpty() && 
            	!args[2].toLowerCase().equals("not specified")) {
            	textOutput = true
            	outputFileName = args[2]
            }
			
			int minNeighbourhood = Integer.parseInt(args[3])
			if (minNeighbourhood < 1) { minNeighbourhood = 1 }
			int maxNeighbourhood = Integer.parseInt(args[4])
			if (maxNeighbourhood > Math.max(rows, cols)) { maxNeighbourhood = Math.max(rows, cols) }
			int neighbourhoodStep = Integer.parseInt(args[5])
			if (neighbourhoodStep < 1) { neighbourhoodStep = 1 }
			
			double minValue = image.getMinimumValue()
			double maxValue = image.getMaximumValue()
			double range = maxValue - minValue
			double K = minValue + range / 2.0
			
			double[][] I = new double[rows][cols]
			double[][] I2 = new double[rows][cols]
			int[][] IN = new int[rows][cols]
			double[][] maxVal = new double[rows][cols]
			
			// calculate the integral image
			int progress = 0
			int oldProgress = -1
			double z, sum, sumN, sumSqr
			for (row = 0; row < rows; row++) {
				sum = 0
				sumSqr = 0
  				sumN = 0
				for (col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z == nodata) {
  						z = 0
  					} else {
  						z = z - K
  						sumN++
  					}
  					sum += z
  					sumSqr += z * z
  					if (row > 0) {
  						I[row][col] = sum + I[row - 1][col]
  						I2[row][col] = sumSqr + I2[row - 1][col]
						IN[row][col] = (int)(sumN + IN[row - 1][col])
  					} else {
						I[row][col] = sum
						I2[row][col] = sumSqr
						IN[row][col] = (int)sumN
  					}
  					maxVal[row][col] = Double.NEGATIVE_INFINITY
  				}
  				progress = (int)(100f * row / rows)
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


			String domainLabel = "Neighbourhood Radius (cells)"
			
			String rangeLabel = "Deviation From Mean (SD)"
//			String zUnits = rasters[0].getZUnits()
//			if (zUnits == null || !zUnits.toLowerCase().equals("not specified")) {
//				rangeLabel = rangeLabel + " (${zUnits})"
//			}
			
			ShapeFile input = new ShapeFile(inputPoints)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POINT) {
            	pluginHost.showFeedback("The input shapefile should be of a POINT ShapeType.")
            	return
            }

            XYSeriesCollection xyCollection = new XYSeriesCollection();

			StringBuilder sb = new StringBuilder()
        	
			int numFeatures = input.getNumberOfRecords()
			int numPoints = 0
			for (ShapeFileRecord record : input.records) {
				int recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				for (int p = 0; p < points.length; p++) {
					numPoints++
				}
			}

			String[] labels = new String[numPoints]
			if (!labelFieldName.isEmpty()) {
				int k = 0
				for (ShapeFileRecord record : input.records) {
					int recNum = record.getRecordNumber()
                	points = record.getGeometry().getPoints()
                	String lbl = (String)input.getAttributeTable().getValue(recNum - 1, labelFieldName);
					for (int p = 0; p < points.length; p++) {
						labels[k] = lbl
						k++
					}
				}
			} else {
				for (int k = 0; k < numPoints; k++) {
					labels[k] = "Point ${k + 1}"
				}
			}
			

			int numLines = 0
			for (int neighbourhood = minNeighbourhood; neighbourhood <= maxNeighbourhood; neighbourhood += neighbourhoodStep) {
				numLines++
			}

			double[][] outputData = new double[numLines][numPoints + 2]
			
            int count = 0
            int numSeries = 0
            oldProgress = -1
            //sb.append("NeighbourhoodSize,DeviationFromMean\n")
			for (ShapeFileRecord record : input.records) {
				int recNum = record.getRecordNumber()
                points = record.getGeometry().getPoints()
				for (int p = 0; p < points.length; p++) {
					XYSeries series = new XYSeries(labels[numSeries]) //("Point ${numSeries + 1}")
					numSeries++

					col = image.getColumnFromXCoordinate(points[p][0])
                    row = image.getRowFromYCoordinate(points[p][1])
                    z = image.getValue(row, col)
                    if (z != nodata) {
                    	int lineNum = 0
                    	for (int neighbourhood = minNeighbourhood; neighbourhood <= maxNeighbourhood; neighbourhood += neighbourhoodStep) {
                    		y1 = row - neighbourhood
							if (y1 < 0) { y1 = 0 }
							if (y1 >= rows) { y1 = rows - 1 }
			
							y2 = row + neighbourhood
							if (y2 < 0) { y2 = 0 }
							if (y2 >= rows) { y2 = rows - 1 }

							x1 = col - neighbourhood
							if (x1 < 0) { x1 = 0 }
							if (x1 >= cols) { x1 = cols - 1 }
		
							x2 = col + neighbourhood
							if (x2 < 0) { x2 = 0 }
							if (x2 >= cols) { x2 = cols - 1 }
								
							N = IN[y2][x2] + IN[y1][x1] - IN[y1][x2] - IN[y2][x1]
							if (N > 0) {
								sum = I[y2][x2] + I[y1][x1] - I[y1][x2] - I[y2][x1]
								sumSqr = I2[y2][x2] + I2[y1][x1] - I2[y1][x2] - I2[y2][x1]
								v = (sumSqr - (sum * sum) / N) / N
								if (v > 0) {
									s = Math.sqrt(v)
									m = sum / N
									outValue = ((z - K) - m) / s
									series.add(neighbourhood, outValue)

									if (textOutput) {
										//sb.append("${neighbourhood * 2 + 1},${outValue}\n")
										if (numSeries == 1) {
											outputData[lineNum][0] = neighbourhood
											outputData[lineNum][1] = outputData[lineNum][0] * resolution
										}
										outputData[lineNum][numSeries + 1] = outValue
									}
								}
							}
							lineNum++
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
            "Topographic Position Scale Signature",
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

//	        NumberAxis range2 = (NumberAxis) plot.getRangeAxis();
//	        range2.setRange(-3.0, 3.0);

            ChartPanel chartPanel = new ChartPanel(chart)
			chartPanel.setPreferredSize(new Dimension(700, 300))
			
        	pluginHost.returnData(chartPanel)

        	if (textOutput) {
        		sb.append("Window Radius,Distance")
        		for (int j = 0; j < numPoints; j++) {
        			sb.append(",${labels[j]}")
        		}
        		sb.append("\n")
        		
        		for (int i = 0; i < numLines; i++) {
        			for (int j = 0; j < numPoints + 2; j++) {
        				sb.append(outputData[i][j]).append(",")
        			}
        			sb.append("\n")
        		}
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
	def f = new LocalTopographicPositionScaleSignature(pluginHost, args, descriptiveName)
}
