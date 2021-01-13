/*
 * Copyright (C) 2014 Dr. Ehsan Roshani and Dr. John Lindsay <jlindsay@uoguelph.ca>
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
import whitebox.geospatialfiles.ShapeFileInfo
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import com.vividsolutions.jts.geom.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.structures.KdTree
import groovy.transform.CompileStatic
import whitebox.stats.Kriging
import whitebox.stats.Kriging.*
import whitebox.stats.Kriging.Variogram
import whitebox.stats.KrigingPoint

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "KrigingInterpolation"
def descriptiveName = "Kriging Interpolation"
def description = "Performs an Kriging interpolation"
def toolboxes = ["Interpolation"]

public class KrigingInterpolation implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
    private Kriging k = null
    private Variogram v = null

	public KrigingInterpolation(WhiteboxPluginHost pluginHost, 
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
        	DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and Value field.", "Input Value Field:", false)															
			DialogCheckBox dcb = sd.addDialogCheckBox("Use z-values", "Use z-values", false)
            dcb.setVisible(false)
            DialogFile dfOutput = sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, false)								
			
			//DialogCheckBox chxError = sd.addDialogCheckBox("Saves the Kriging Variance Error Map", "Save Variance Error Raster", false)												
			DialogFile dfError = sd.addDialogFile("Error Output file", "Output Error Raster File (optional):", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, true)						
			//dfError.visible = false
			
			DialogDataInput txtCellSize = sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, true)														
			sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)

			DialogDataInput txtNNeighbor = sd.addDialogDataInput("Enter number of neighbors.", "Number of Neighbors:", "10", true, false)											
			
			sd.addDialogLabel(" ")
            sd.addDialogLabel("<html><b>Variogram Parameters</b></html>")
            
			DialogComboBox cboxSVType = sd.addDialogComboBox("Enter variogram model type", "Model Type:", ["Gaussian", "Exponential", "Spherical"], 0)	
			DialogDataInput txtNlag = sd.addDialogDataInput("Number of Lags.", "Number of Lags:", "12", true, false)																
			DialogDataInput txtLagSize = sd.addDialogDataInput("Lag Size.", "Lag Size:", "", true, false)																			
			
			DialogCheckBox chxNugget = sd.addDialogCheckBox("Apply nugget in the theoretical variogram calculation ", "Apply Nugget:", true)									

			DialogCheckBox chxAnIsotropic = sd.addDialogCheckBox("Is the data Anisotropic", "Use Anisotropic Model:", false)														//15
			DialogDataInput txtAngle = sd.addDialogDataInput("Angle (rad).", "Angle (Rad):", "", true, true)																		//16
			DialogDataInput txtTolerance = sd.addDialogDataInput("Tolerance.", "Tolerance (Rad):", "", true, true)															//17
			DialogDataInput txtBandWidth = sd.addDialogDataInput("Band Width.", "Band Width:", "", true, true)																	//18
			txtAngle.visible = false
			txtTolerance.visible = false
			txtBandWidth.visible = false

			def btn = sd.addDialogButton("View Variogram", "Center")
			
			DialogCheckBox chxCurve = sd.addDialogCheckBox("Show variogram curve", "Show Variogram:", true)															//13
			DialogCheckBox chxMap = sd.addDialogCheckBox("Show variogram map", "Show Variogram Map:", true)																//14

			btn.addActionListener(new ActionListener() {
 	            public void actionPerformed(ActionEvent e) {
				    k = new Kriging()

					addPropertyListenerToKriging(k)
					
					if (Boolean.parseBoolean(chxNugget.getValue())) {
						k.ConsiderNugget = true;
					} else {
						k.ConsiderNugget = false;
					}

					String[] inputData = dfs.getValue().split(";")
					if (inputData[0] == null || inputData[0].isEmpty()) {
						pluginHost.showFeedback("Please specify a shapefile and attribute.")
						return
					}
					if (!Boolean.parseBoolean(dcb.getValue())) {
						if (inputData.length < 2 || inputData[1] == null || inputData[1].isEmpty()) {
							pluginHost.showFeedback("Please specify a shapefile and attribute.")
							return
						}
					}
					if (!Boolean.parseBoolean(dcb.getValue()) &&dfs.getValue().length()>2) {
						k.readPointFile(inputData[0], inputData[1]);
					} else if (Boolean.parseBoolean(dcb.getValue())) {
						k.readPointFile(inputData[0]);
					}
			        
			        k.LagSize =  (txtLagSize.getValue()).toDouble()
			        if (Boolean.parseBoolean(chxAnIsotropic.getValue())) {
						k.Anisotropic = true;
						k.Angle = Double.parseDouble(txtAngle.getValue())
						k.BandWidth = Double.parseDouble(txtBandWidth.getValue())
						k.Tolerance = Double.parseDouble(txtTolerance.getValue())
					} else {
						k.Anisotropic = false;
					}
			
			        int numLags = Integer.parseInt(txtNlag.getValue())
			        boolean anisotropic = Boolean.parseBoolean(chxAnIsotropic.getValue())
			        
			        if ((cboxSVType.getValue()).toLowerCase().contains("gaussian")) {
						v = k.getSemivariogram(Kriging.SemivariogramType.GAUSSIAN, 1d, numLags , anisotropic, true);
					} else if ((cboxSVType.getValue()).toLowerCase().contains("exponential")) {
						v = k.getSemivariogram(Kriging.SemivariogramType.EXPONENTIAL, 1d, numLags , anisotropic, true);
					} else {
						v = k.getSemivariogram(Kriging.SemivariogramType.SPHERICAL, 1d, numLags , anisotropic, true);
					}
			        
					if (Boolean.parseBoolean(chxCurve.getValue())) {
			        	k.DrawSemivariogram(k.bins, v);
					}
					
			        if (Boolean.parseBoolean(chxMap.getValue())) {        		
						k.calcBinSurface(k.SVType, 1, k.NumberOfLags, Boolean.parseBoolean(chxAnIsotropic.getValue()));
						k.DrawSemivariogramSurface(k.LagSize * (k.NumberOfLags), Boolean.parseBoolean(chxAnIsotropic.getValue()));
					}
	            }
        	});      
            
			//Listener for chxError            
            def lstrAnIso = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = chxAnIsotropic.getValue()
            		if (!value.isEmpty()&& value != null) { 
            			if ( chxAnIsotropic.getValue() == "true") {
            				txtAngle.setValue("0.0")
            				txtTolerance.setValue("1.5707")
            				txtBandWidth.setValue("5")
							txtAngle.visible = true
							txtTolerance.visible = true
							txtBandWidth.visible = true
		            	} else {
		            		txtAngle.setValue("")
            				txtTolerance.setValue("")
            				txtBandWidth.setValue("")
							txtAngle.visible = false
							txtTolerance.visible = false
							txtBandWidth.visible = false
		            	}
            		}
            	} 
            } as PropertyChangeListener
            chxAnIsotropic.addPropertyChangeListener(lstrAnIso)

            //Listener for dfs, It updates the lag size and boundaries 
            def lsnField = { evt -> if (evt.getPropertyName().equals("value")) { 
            		
            		String value = dfs.getValue()
            		if (value != null && !value.isEmpty()) { 
            			value = value.trim()
            			String[] strArray = dfs.getValue().split(";")
            			String fileName = strArray[0]
            			File file = new File(fileName)
            			ShapeFile shapefile = new ShapeFile(fileName)

						double difY = shapefile.getyMax() - shapefile.getyMin()
						double difX = shapefile.getxMax() - shapefile.getxMin()
						double maxL = Math.max(difX, difY)
						txtLagSize.setValue((maxL / (txtNlag.getValue()).toDouble()).toString())
						if (file.exists()) {
	            			ShapeFileInfo shapefileInfo = new ShapeFileInfo(fileName)
		            		if (shapefileInfo.getShapeType().getDimension() == ShapeTypeDimension.Z) {
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
            dfs.addPropertyChangeListener(lsnField)
            
            sd.visible = true
        }
    }

    @CompileStatic
    private void execute(String[] args) {
		try {
			
			double cellSize = -1.0
			String str
			
			// read the input variables
			String[] inputData = args[0].split(";")
			if (inputData[0] == null || inputData[0].isEmpty()) {
				pluginHost.showFeedback("Input shapefile and attribute not specified.")
				return
			}
//			if (inputData.length < 2 || inputData[1] == null || inputData[1].isEmpty()) {
//				pluginHost.showFeedback("Input shapefile and attribute not specified.")
//				return
//			}
			String inputFile = inputData[0]
			boolean useZValues = Boolean.parseBoolean(args[1])
			String outputFile = args[2]
			if (outputFile == null || outputFile.isEmpty()) {
				pluginHost.showFeedback("Output file not specified.")
				return
			}
			String outputErrorFile = args[3]
			str = args[4]
			if (str == null || str.isEmpty()) {
	        	str = "not specified"
	        }
	        if (!str.toLowerCase().contains("not")) {
	            cellSize = Double.parseDouble(args[5]);
	        }
	        String baseFileHeader = args[6]
	        if (baseFileHeader == null || baseFileHeader.isEmpty()) {
	        	baseFileHeader = "not specified"
	        }
	        int numNeighbours = Integer.parseInt(args[7])
			if (numNeighbours < 0) {
				numNeighbours = 0
			}
			String modelType = "spherical"
			if (args[8].toLowerCase().contains("gauss")) {
				modelType = "gaussian"
			} else if (args[8].toLowerCase().contains("expon")) {
				modelType = "exponential"
			}
			int numLags = Integer.parseInt(args[7])
			double lagSize = Double.parseDouble(args[8])
			boolean applyNugget = Boolean.parseBoolean(args[9])
			boolean anisotropic = Boolean.parseBoolean(args[10])
			double angle = 0.0
			double tolerance = 0.0
			double bandWidth = 0.0
			if (anisotropic) {
				str = args[11]
				if (str != null && str.isEmpty() && !str.toLowerCase().contains("not")) {
					angle = Double.parseDouble(args[11])
				} else {
					pluginHost.showFeedback("Error: Anisotropic parameter not set properly")
					return
				}
				str = args[12]
				if (str != null && str.isEmpty() && !str.toLowerCase().contains("not")) {
					tolerance = Double.parseDouble(args[12])
				} else {
					pluginHost.showFeedback("Error: Anisotropic parameter not set properly")
					return
				}
				str = args[13]
				if (str != null && str.isEmpty() && !str.toLowerCase().contains("not")) {
					bandWidth = Double.parseDouble(args[13])
				} else {
					pluginHost.showFeedback("Error: Anisotropic parameter not set properly")
					return
				}
			}
			boolean showSemivariogram = Boolean.parseBoolean(args[14])
			boolean showSemivariogramMap = Boolean.parseBoolean(args[15])

			if (k == null || v == null) {
				pluginHost.updateProgress("Calculating Semivariogram...", 0)
				k = new Kriging()

				addPropertyListenerToKriging(k)
				k.ConsiderNugget = applyNugget
				k.readPointFile(inputData[0], inputData[1])
		        k.LagSize =  lagSize
		        k.Anisotropic = anisotropic
		        if (anisotropic) {
					k.Angle = angle
					k.BandWidth = bandWidth
					k.Tolerance = tolerance
				}

				switch (modelType) {
					case "gaussian":
						v = k.getSemivariogram(Kriging.SemivariogramType.GAUSSIAN, 1d, numLags, anisotropic, true);
						break
					case "exponential":
						v = k.getSemivariogram(Kriging.SemivariogramType.EXPONENTIAL, 1d, numLags, anisotropic, true);
						break
					default: // spherical
						v = k.getSemivariogram(Kriging.SemivariogramType.SPHERICAL, 1d, numLags, anisotropic, true);
				}
			} else if (k != null && v != null) {
				if (k.LagSize != lagSize || k.NumberOfLags != numLags) {
					pluginHost.updateProgress("Calculating Semivariogram...", 0)
					k = new Kriging()
	
					addPropertyListenerToKriging(k)
					k.ConsiderNugget = applyNugget
					k.readPointFile(inputData[0], inputData[1])
			        k.LagSize =  lagSize
			        k.Anisotropic = anisotropic
			        if (anisotropic) {
						k.Angle = angle
						k.BandWidth = bandWidth
						k.Tolerance = tolerance
					}
	
					switch (modelType) {
						case "gaussian":
							v = k.getSemivariogram(Kriging.SemivariogramType.GAUSSIAN, 1d, numLags, anisotropic, true);
							break
						case "exponential":
							v = k.getSemivariogram(Kriging.SemivariogramType.EXPONENTIAL, 1d, numLags, anisotropic, true);
							break
						default: // spherical
							v = k.getSemivariogram(Kriging.SemivariogramType.SPHERICAL, 1d, numLags, anisotropic, true);
					}
				}
			}

			ShapeFile input = new ShapeFile(inputFile)

			double north, south, east, west

			north = input.getyMax() + cellSize / 2.0;
	        south = input.getyMin() - cellSize / 2.0;
	        east = input.getxMax() + cellSize / 2.0;
	        west = input.getxMin() - cellSize / 2.0;

            // initialize the output raster
            WhiteboxRaster output;
            double nodata = -32768.0d
	        if (baseFileHeader.toLowerCase().contains("not specified")) {
	        	if (cellSize < 0) {
	        		cellSize = Math.max((input.getxMax() - input.getxMin()) / 998, (input.getyMax() - input.getyMin()) / 998)
	        	}
	        	
	        	int rows = (int) (Math.ceil((north - south) / cellSize));
                int cols = (int) (Math.ceil((east - west) / cellSize));

                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;

				output = new WhiteboxRaster(outputFile, north, south, east, west,
                    rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                    WhiteboxRasterBase.DataType.FLOAT, nodata, nodata);
                output.setPreferredPalette("spectrum.plt")
	        } else {
	        	output = new WhiteboxRaster(outputFile, "rw",
	                baseFileHeader, WhiteboxRasterBase.DataType.FLOAT, 
	                nodata);
	        }

			if (!outputErrorFile.toLowerCase().equals("not specified")) {
        		WhiteboxRaster errorOutput = new WhiteboxRaster(outputErrorFile, "rw",
	                outputFile, WhiteboxRasterBase.DataType.FLOAT, 
	                nodata);
	            errorOutput.setPreferredPalette("spectrum.plt")
				
				k.interpolateRaster(v, numNeighbours, output, errorOutput)
				
				errorOutput.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
		        errorOutput.addMetadataEntry("Created on " + new Date())
				
				errorOutput.addMetadataEntry("Semivariogram Model = " + v.Type)
				errorOutput.addMetadataEntry("Range = " + v.Range)
				errorOutput.addMetadataEntry("Sill = " + v.Sill)
				errorOutput.addMetadataEntry("Nugget = " + v.Nugget)
				errorOutput.addMetadataEntry("Lag Size = " + lagSize)
				errorOutput.addMetadataEntry("Num. Bins = " + numLags)
				errorOutput.addMetadataEntry("RMSE = " + Math.sqrt(v.mse))
				errorOutput.close()
			} else {
				k.interpolateRaster(v, numNeighbours, output, false)
			}
			
			if (showSemivariogram) {
	        	k.DrawSemivariogram(k.bins, v)
			}
			
	        if(showSemivariogramMap) {        		
				k.calcBinSurface(k.SVType, 1, numLags, anisotropic)
				k.DrawSemivariogramSurface(lagSize * numLags, anisotropic)
			}			
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			
			output.addMetadataEntry("Semivariogram Model = " + v.Type)
			output.addMetadataEntry("Range = " + v.Range)
			output.addMetadataEntry("Sill = " + v.Sill)
			output.addMetadataEntry("Nugget = " + v.Nugget)
			output.addMetadataEntry("Lag Size = " + lagSize)
			output.addMetadataEntry("Num. Bins = " + numLags)
			output.addMetadataEntry("RMSE = " + Math.sqrt(v.mse))
			output.close()

			pluginHost.returnData(outputFile)
			if (!outputErrorFile.toLowerCase().equals("not specified")) {
				pluginHost.returnData(outputErrorFile)
			}

			DecimalFormat df = new DecimalFormat("###,###,###,##0.000")
			
			StringBuilder ret = new StringBuilder()
			ret.append("<!DOCTYPE html>")
			ret.append('<html lang="en">')

			ret.append("<head>")
//            ret.append("<meta content=\"text/html; charset=iso-8859-1\" http-equiv=\"content-type\">")
            ret.append("<title>Kriging Report</title>").append("\n")
            
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
            ret.append("<body><h1>Kriging Report</h1>").append("\n")

        
			ret.append("<p><b>Input File Name:</b> &nbsp ").append(new File(inputFile).getName()).append("<br>\n")
			ret.append("<b>Output File Name:</b> &nbsp ").append(new File(outputFile).getName()).append("\n")
			if (!outputErrorFile.toLowerCase().equals("not specified")) {
				ret.append("<br><b>Output Error File Name:</b> &nbsp ").append(new File(outputErrorFile).getName()).append("</p>\n")
			} else {
				ret.append("</p>")
			}
			
			ret.append("<br><table border=\"1\" cellspacing=\"0\" cellpadding=\"3\">").append("\n")
			
			ret.append("<tr><th>Characteristic</th><th>Value</th></tr>")

			ret.append("<tr><td>Semivariogram Model</td><td class=\"numberCell\">")
			ret.append(v.Type).append("</td></tr>").append("\n")

			ret.append("<tr><td>Range</td><td class=\"numberCell\">")
			ret.append(df.format(v.Range)).append("</td></tr>").append("\n")

			ret.append("<tr><td>Sill</td><td class=\"numberCell\">")
			ret.append(df.format(v.Sill)).append("</td></tr>").append("\n")

			ret.append("<tr><td>Nugget</td><td class=\"numberCell\">")
			ret.append(df.format(v.Nugget)).append("</td></tr>").append("\n")

			ret.append("<tr><td>Lag Size</td><td class=\"numberCell\">")
			ret.append(df.format(lagSize)).append("</td></tr>").append("\n")

			ret.append("<tr><td>Number of Bins</td><td class=\"numberCell\">")
			ret.append(numLags).append("</td></tr>").append("\n")

			ret.append("<tr><td>RMSE</td><td class=\"numberCell\">")
			ret.append(df.format(Math.sqrt(v.mse))).append("</td></tr>").append("\n")

			ret.append("</table>")

			ret.append("</body></html>")
			pluginHost.returnData(ret.toString());
			
		} catch (Exception e) {
			pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
	        pluginHost.logException("Error in " + descriptiveName, e)
		} finally {
			pluginHost.updateProgress("Progress", 0)
		}
    }

	private void addPropertyListenerToKriging(Kriging krig) {
		def progressListener = { evt -> if (evt.getPropertyName().equals("progress")) { 
        		int progress = (int)evt.getNewValue()
        		pluginHost.updateProgress("Interpolating Data:", progress)
        	} 
        } as PropertyChangeListener
        krig.addPropertyChangeListener(progressListener)
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
    def f = new KrigingInterpolation(pluginHost, args, name, descriptiveName)
}
