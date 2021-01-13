/* Early draft of the Roughness From Point Cloud Profiles tool.
 *  Kept as back-up. Dr. John Lindsay.
 */
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.Date
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.awt.Dimension
import java.awt.Color;
import javax.swing.JPanel;
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.*
import groovy.transform.CompileStatic
import whitebox.structures.BoundingBox
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.ShapeType
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.geospatialfiles.shapefile.attributes.DBFField
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.structures.KdTree
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
def name = "ProfileRoughnessLiDAR"
def descriptiveName = "Profile Roughness (LiDAR)"
//def description = "Extracts profile roughness from LiDAR data."
//def toolboxes = ["LidarTools"]

public class ProfileRoughnessLiDAR implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private String descriptiveName
	private ScriptDialog sd;
//	private boolean something = false;
	
	public ProfileRoughnessLiDAR(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
		this.pluginHost = pluginHost;
		this.descriptiveName = descriptiveName;
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
			sd.addDialogFile("Input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
			sd.addDialogFile("Input Vector Polygon for clipping point cloud (optional)", "Input Vector Clipping Polygon File (Optional):", "open", "Vector Files (*.shp), SHP", true, true)
			sd.addDialogDataInput("Enter Profile Azimuth, in degrees 0-360", "Profile Azimuth (degrees):", "90.0", true, false)
            sd.addDialogDataInput("Enter Profile Spacing", "Profile Spacing (m):", "0.25", true, false)
            sd.addDialogDataInput("Enter Profile Sampling Density", "Profile Sampling Density (m):", "0.1", true, false)
            DialogComboBox comboBox = sd.addDialogComboBox("Enter Roughness Statistic ", "Roughness Statistic:", ["RMS Height", "Correlation Length", "Shadow Roughness Index"], 0)
            DialogDataInput incidenceBox = sd.addDialogDataInput("Enter Satellite Incidence Angle", "Incidence Angle (degrees):", "30", true, true)
            incidenceBox.visible = false
			DialogDataInput wavelengthBox = sd.addDialogDataInput("Enter Satellite Wavelength", "Wavelength (cm):", "5", true, true)
            wavelengthBox.visible = false

            def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = comboBox.getValue()
            		if (!value.isEmpty() && value != null) { 
            			if (value.toLowerCase().contains("shadow")) {
            				incidenceBox.visible = true
            				wavelengthBox.visible = true
		            	} else {
		            		incidenceBox.visible = false
		            		wavelengthBox.visible = false
		            	}
            		} 
            	} 
            } as PropertyChangeListener
            comboBox.addPropertyChangeListener(lstr)

            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}



	// The execute function is the main part of the tool, where the actual
	// work is completed.
	@CompileStatic
	private void execute(String[] args) {
		try {
			int progress, oldProgress
			int a 
	  		double x, y, z, rotatedX, rotatedY
	  		
			if (args.length != 8) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			boolean isClipFileSpecified = false
			String clipFile 
			if (!args[1].toLowerCase().contains("not specified") && args[1].trim().length() > 0) {
				isClipFileSpecified = true
				clipFile = args[1]
			}
			double profileAzimuth = Math.toRadians(Double.parseDouble(args[2]))
			double profileSpacing = Double.parseDouble(args[3])
			double samplingDensity = Double.parseDouble(args[4])
			String statType = "rms"
			if (args[5].toLowerCase().contains("cor")) {
				statType = "cor"
			} else if (args[5].toLowerCase().contains("shadow")) {
				statType = "shadow"
			}
			double incidenceAngle = 0
			if (!args[6].toLowerCase().contains("not specified")) {
				incidenceAngle = Math.toRadians(Double.parseDouble(args[6]))
			}
			double wavelength = 0
			if (!args[7].toLowerCase().contains("not specified")) {
				wavelength = (Double.parseDouble(args[7]) / 100)
			}
			
			if (samplingDensity >= profileSpacing) { 
				pluginHost.showFeedback("Warning: The Profile Spacing should be larger than the Sampling Density")
			}

			pluginHost.updateProgress("Reading LAS Dataset:", 0)
			LASReader las = new LASReader(inputFile)
			BoundingBox extent = new BoundingBox(las.getMinX(), las.getMinY(), las.getMaxX(), las.getMaxY());
			int numPoints = (int)las.getNumPointRecords()
			double midX = las.getMinX() + (las.getMaxX() - las.getMinX()) / 2.0
			double midY = las.getMinY() + (las.getMaxY() - las.getMinY()) / 2.0
			def outputFile = pluginHost.getWorkingDirectory() + "rotatedPoints2.shp"
			DBFField[] fields = new DBFField[1];

            fields[0] = new DBFField();
            fields[0].setName("FID");
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
            fields[0].setFieldLength(10);
            fields[0].setDecimalCount(0);

			def output = new ShapeFile(outputFile, ShapeType.MULTIPOINT, fields);
			PointsList points = new PointsList()
							
			double[] entry;
			double minXRotated = Double.POSITIVE_INFINITY
			double maxXRotated = Double.NEGATIVE_INFINITY
			double minYRotated = Double.POSITIVE_INFINITY
			double maxYRotated = Double.NEGATIVE_INFINITY
			ArrayList<PointRecord> points2 = las.getPointRecordsInBoundingBox(extent)
	 		if (!isClipFileSpecified) {
	 			oldProgress = 0
	 			a = 0
		 		for (PointRecord point : points2) {
		 			a++
		 			x = point.getX()
		 			y = point.getY()
		 			z = point.getZ()
			 		x = x - midX;
	                y = y - midY;
	                rotatedX = (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
	                rotatedY = (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
//	                points.addPoint(rotatedX + midX, rotatedY + midY) 
					if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
		 					
					
		 			progress = (int)(100f * a / numPoints)
					if (progress != oldProgress) {
						pluginHost.updateProgress("Reading LAS Dataset:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
		 		}
	 		} else {
	 			
				def clipPoly = new ShapeFile(clipFile)
				ShapeType shapeType = clipPoly.getShapeType()
	
				if (shapeType.getBaseType() != ShapeType.POLYGON) {
					pluginHost.showFeedback("The input shapefile must be of a POLYGON base ShapeType.")
	                return
				}
		
				double[][] nodes = clipPoly.getRecord(0).getGeometry().getPoints()
				int numNodes = nodes.length;
				int m;
			    int n;
			    boolean result = false;
				oldProgress = 0
				a = 0
	 			for (PointRecord point : points2) {
	 				a++
		 			x = point.getX()
		 			y = point.getY()
		 			z = point.getZ()
		 			n = numNodes - 1
		 			result = false;
		 			for (m = 0; m < numNodes; n = m++) {
				    	if ((nodes[m][1] > y) != (nodes[n][1] > y) &&
				            (x < (nodes[n][0] - nodes[m][0]) * (y - nodes[m][1]) / (nodes[n][1]-nodes[m][1]) + nodes[m][0])) {
				          result = !result;
				        }
				    }
				    if (result) { 
				 		x = x - midX;
		                y = y - midY;
		                rotatedX = (x * Math.cos(profileAzimuth)) - (y * Math.sin(profileAzimuth));
		                rotatedY = (x * Math.sin(profileAzimuth)) + (y * Math.cos(profileAzimuth));
//						points.addPoint(rotatedX + midX, rotatedY + midY) 
if (rotatedX < minXRotated) { minXRotated = rotatedX }
					if (rotatedX > maxXRotated) { maxXRotated = rotatedX }
					if (rotatedY < minYRotated) { minYRotated = rotatedY }
					if (rotatedY > maxYRotated) { maxYRotated = rotatedY }
		 			
				    }
				    progress = (int)(100f * a / numPoints)
					if (progress != oldProgress) {
						pluginHost.updateProgress("Reading LAS Dataset:", progress)
						oldProgress = progress
						if (pluginHost.isRequestForOperationCancelSet()) {
	                        pluginHost.showFeedback("Operation cancelled")
							return
	                    }
					}
		 		}
	 		}

	 		Object[] recData = new Object[1]
			recData[0] = new Double(1)
		
//			whitebox.geospatialfiles.shapefile.Geometry g = new MultiPoint(points.getPointsArray())
//        	output.addRecord(g, recData)
//				
//			output.write()
//			pluginHost.returnData(outputFile)
		

			double deltaXRotated = maxXRotated - minXRotated
			double deltaYRotated = maxYRotated - minYRotated
			int numProfiles = (int) (deltaXRotated / profileSpacing)
			int numSamples = (int) (deltaYRotated / samplingDensity)
			double dist
			double noData = -32768.0
			double[] statValues = new double[numProfiles]
			oldProgress = 0
			
			for (int column = 0; column < numProfiles; column++) {
				double[] zValues = new double[numSamples]
				for (int row = 0; row < numSamples; row++) {
					x = minXRotated + column * profileSpacing
					y = maxYRotated - row * samplingDensity
					 rotatedX = (x * Math.cos(-profileAzimuth)) - (y * Math.sin(-profileAzimuth));
	                rotatedY = (x * Math.sin(-profileAzimuth)) + (y * Math.cos(-profileAzimuth));

//					entry = [y, x];
//                    results = pointsTree.nearestNeighbor(entry, 1, false)
//                    dist = results.get(0).distance
//                    if (samplingDensity > dist) {
//	                    zValues[row] = (double) results.get(0).value
//                    } else {
//                    	zValues[row] = noData
//                    }
                   points.addPoint(rotatedX + midX, rotatedY + midY) 
	 
				}
//
//				switch (statType) {
//					case "rms":
//						//calculate RMS for profile
//						statValues[column] = getRmsh(zValues, samplingDensity)
//						break;
//
//					case "cor":
//						//calculate CL for profile
//
//						break;
//
//					case "shadow":
//						// calculate the proportional shadow based rougness for profile
//						statValues[column] = getSri(zValues, samplingDensity, incidenceAngle, wavelength)
//						if (column == 10) {
//							flag = true;
//						}
//						break;
//				}
				progress = (int)(100f * column / (numProfiles - 1))
				if (progress != oldProgress) {
					pluginHost.updateProgress("Calculating Statistic:", progress)
					oldProgress = progress
					if (pluginHost.isRequestForOperationCancelSet()) {
                        pluginHost.showFeedback("Operation cancelled")
						return
                    }
				}
			}
//			
//			
//
//		XYSeriesCollection xyCollection = new XYSeriesCollection();
//		XYSeries series = new XYSeries("")
//		//String outStr = ""
//		for (int row = 0; row < numProfiles; row++) {
//            x = row;
//            y = statValues[row];
//            //outStr += statValues[a]+ "\t" + x + "\t" + y + "\n"
//            if (y != noData) {
//            	series.add(x, y)	
//            } else {
//            	if (series.getItemCount() > 0) {
//	            	xyCollection.addSeries(series);	
//	            	series = new XYSeries("")
//            	}	
//            }
//        }
//        xyCollection.addSeries(series);
//        //pluginHost.returnData(outStr)
//
//       String chartTitle
//       String rangeLabel
//       String domainLabel = "Profile Number"
//		switch (statType) {
//			case "rms":
//				chartTitle = "Root Mean Squared Height"
//				rangeLabel = "RMS Height (Z units)"
//				break
//
//			case "cor":
//				chartTitle = "Roughness Correlation Length"
//				rangeLabel = "Correlaion Length (Map units)"
//				break
//				//TO DO: Add case for proportional shadow length stat
//		}
//
//       
//		JFreeChart chart = ChartFactory.createXYLineChart(
//        chartTitle,
//        domainLabel,
//        rangeLabel,
//        xyCollection,
//        PlotOrientation.VERTICAL,  // Plot Orientation
//        false,                      // Show Legend
//        true,                      // Use tooltips
//        false                      // Configure chart to generate URLs?
//        );
//
//		XYPlot plot = (XYPlot)chart.getPlot();
//		plot.setBackgroundPaint(Color.WHITE); 
//		plot.setDomainGridlinePaint(Color.BLACK); 
//		plot.setRangeGridlinePaint(Color.BLACK); 
//
		// TODO: change the line colour of each series to red.

//		ChartPanel chartPanel = new ChartPanel(chart)
//		chartPanel.setPreferredSize(new Dimension(700, 300))
//		
//    	pluginHost.returnData(chartPanel)
whitebox.geospatialfiles.shapefile.Geometry g = new MultiPoint(points.getPointsArray())
        	output.addRecord(g, recData)
				
			output.write()
			pluginHost.returnData(outputFile)
		

					
		} catch (Exception e) {
	    	pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
        } finally {
        	// reset the progress bar
        	pluginHost.updateProgress(0)
        }
	}

	boolean flag = false;
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
    
	@CompileStatic
    public double getRmsh(double[] zValues, double samplingDensity) {
    	double result
        double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;
        long N = 0;
        double noData = -32768.0
        double x, y
        int row
        double slope, intercept, yHat
           
		for (row = 0; row < zValues.length; row++) {
            x = row * samplingDensity;
            y = zValues[row];
            if (y != noData) {
                sumX += x;
                sumY += y;
                sumXY += x * y;
                sumXX += x * x;
                sumYY += y * y;
                N++;
            }
        }
        if (N > 0) {
        	slope = (N * sumXY - (sumX * sumY)) / (N * sumXX - (sumX * sumX));
        	intercept = (sumY - slope * sumX) / N;
    		sumXX = 0
    		for (row = 0; row < zValues.length; row++) {
    			if (zValues[row] != noData) {
    				yHat = (slope * row * samplingDensity + intercept)
    				sumXX += (zValues[row] - yHat) * (zValues[row] - yHat)
    			}
    		}
    		result = Math.sqrt(sumXX / (N - 1))
        } else {
        	result = noData
        }
       return result
    }

    
	@CompileStatic
    public double getSri(double[] zValues, double samplingDensity, double incidenceAngle, double wavelength) {
    	double noData = -32768.0
        double x, y, z, zN, xN
        int row, rowN
        double dist, deltaZ, angleN
		double[] horizonAngle = new double[zValues.length]
		double[] hALocation = new double[zValues.length]
		double[] elevation = new double[zValues.length]
		
        // Calculate horizon angle for each location
        for (row = 0; row < zValues.length; row++) {
        	horizonAngle[row] = -99;
            x = row * samplingDensity;
            z = zValues[row];
            if (z != noData) {
               for (rowN = row - 1; rowN >= 0; rowN--) {
           			if (rowN > 0) {
               	 		zN = zValues[rowN];
               	 		if (zN != noData) {
		               	 	xN = rowN * samplingDensity;
		               	 	dist = (row - rowN) * samplingDensity; //x - xN;
		               	 	deltaZ = zN - z;
		               	 	angleN = Math.tan(deltaZ/dist)
		               	 	if (angleN > horizonAngle[row]) {
		               	 		horizonAngle[row] = angleN;
		               	 		hALocation[row] = xN
		               	 		elevation[row] = zN
               	 			}
               	 		}

               		}
               }
            }
        }

        if (flag) {
        	flag = false
        	String str = ""
        	for (int j = 0; j < horizonAngle.length; j++) {
        		str += "${horizonAngle[j]}\t${j*samplingDensity - hALocation[j]}\t ${elevation[j]}\t ${elevation[j]-zValues[j]}\n"
        	}
        	pluginHost.returnData(str)
        }

		// Calculate roughness element height
		double numValidElement, numShadowElement
		
		double[] roughnessElementHeight = new double[zValues.length]
		int searchWindowSize = (int)((10 * (wavelength / samplingDensity)) / 2) 
			if (searchWindowSize < 1) {
				searchWindowSize = 1
			}
		for (row = 0; row < zValues.length; row++) {
			z = zValues[row];
			if (z != noData) {
				numValidElement++
				for (rowN = row - searchWindowSize; rowN <= row + searchWindowSize; rowN++) {
					if (rowN >= 0 && rowN < zValues.length) {
						if (zValues[rowN] != noData) {
							if (zValues[rowN] < z) {
								z = zValues[rowN]
							}
						}
					}
				}
				roughnessElementHeight[row] = zValues[row] - z
			}

		}
		// Find Profile locations within shadows
		 boolean[] inShadow = new boolean[zValues.length]
		 for (row = 0; row < zValues.length; row++) {
		 	if (horizonAngle[row] > incidenceAngle && zValues[row] != noData) {
//		 		if (roughnessElementHeight[row] >= (0.5 * wavelength)) {
		 			inShadow[row] = true
		 			numShadowElement++
//		 		}	
		 	}
		 }

		return numShadowElement / numValidElement

    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def myTool = new ProfileRoughnessLiDAR(pluginHost, args, name, descriptiveName)
}
