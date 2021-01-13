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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.structures.BooleanBitArray2D
import groovy.time.*
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "MapGullyDepth"
def descriptiveName = "Map Gully Depth"
def description = "Maps gullies in a DEM and estimates their depth."
def toolboxes = ["RelativeLandscapePosition"]


//String[] args

public class MapGullyDepth implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public MapGullyDepth(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("DEM file", "Digital Elevation Model (DEM) File:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Maximum gully cross-sectional width", "Maximum Gully Width:", "25.0", true, false)
			sd.addDialogDataInput("Maximum gully cross-sectional depth", "Maximum Gully Depth:", "3.0", true, false)
			sd.addDialogDataInput("Minimum gully cross-sectional depth", "Minimum Gully Depth:", "0.25", true, false)
			sd.addDialogDataInput("Difference from mean elevation (DFME) threshold", "DFME Threshold:", "0.0", true, false)
			sd.addDialogDataInput("Plan curvature threshold", "Plan Curvature Threshold:", "500.0", true, false)
			sd.addDialogDataInput("Smoothing parameter", "Smoothing Parameter:", "5", true, false)
			sd.addDialogCheckBox("Background value is NoData", "Use NoData Background Value", true)
            
			// resize the dialog to the standard size and display it
			sd.setSize(800, 400)
			sd.visible = true
		}
	}

	@CompileStatic
	private void execute(String[] args) {
		try {

			if (args.length != 9) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// input parameters
			String demFile = args[0]
			String outputFile = args[1]
			double maxGullyWidth = Double.parseDouble(args[2])
			double maxGullyDepth = Double.parseDouble(args[3])
			double minGullyDepth = Double.parseDouble(args[4])
			double dfmeThreshold = Double.parseDouble(args[5])
			double highPlanCurvValue = Double.parseDouble(args[6])
			int smoothingParameter = Integer.parseInt(args[7])
			boolean nodataForBackground = Boolean.parseBoolean(args[8])
			
			int progress, oldProgress
			def wd = pluginHost.getWorkingDirectory()
			
			// calculate the difference from mean elevation
//			def start = new Date()
			def dfmeFile = wd + "dfme.dep"
			def neighbourhoodSize = String.valueOf((int)(maxGullyWidth / 2))
			args = [demFile, dfmeFile, neighbourhoodSize]
			pluginHost.runPlugin("DifferenceFromMeanElevation", args, false, true)
			
			// calculate the plan curvature on a smoothed DEM
			def planCurvFile
			if (smoothingParameter > 2) {
				def smoothDEMFile = wd + "smooth DEM.dep" 
				String dim = String.valueOf(smoothingParameter)
				def rounded = "true" 
				def reflectEdges = "true" 
				args = [demFile, smoothDEMFile, dim, dim, rounded, reflectEdges] 
				pluginHost.runPlugin("FilterMean", args, false, true)
			
				planCurvFile = wd + "plan curv.dep" 
				args = [smoothDEMFile, planCurvFile, "1.0"]
				pluginHost.runPlugin("PlanCurv", args, false, true)

				new File(smoothDEMFile).delete()
				new File(smoothDEMFile.replace(".dep", ".tas")).delete()
			} else {
				planCurvFile = wd + "plan curv.dep" 
				args = [demFile, planCurvFile, "1.0"]
				pluginHost.runPlugin("PlanCurv", args, false, true)
			
			}
			
			// read in the DFME, DEM, and PlanCurv rasters
			WhiteboxRasterInfo dfme = new WhiteboxRasterInfo(dfmeFile)
			WhiteboxRasterInfo planCurv = new WhiteboxRasterInfo(planCurvFile)
			int rows = dfme.getNumberRows()
			int cols = dfme.getNumberColumns()
			double dfmeNodata = dfme.getNoDataValue()
			double pcNodata = planCurv.getNoDataValue()
			
			// find the areas of low DFME and high plan curv
			BooleanBitArray2D lowDfme = new BooleanBitArray2D(rows, cols)
			BooleanBitArray2D highPc = new BooleanBitArray2D(rows, cols)
			double[] dfmeData
			double[] planCurvData
			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				dfmeData = dfme.getRowValues(row)
				planCurvData = planCurv.getRowValues(row)
				for (int col = 0; col < cols; col++) {
					if (dfmeData[col] <= dfmeThreshold &&
					     dfmeData[col] != dfmeNodata) {
						lowDfme.setValue(row, col, true)
					}
					if (planCurvData[col] >= highPlanCurvValue &&
					     planCurvData[col] != pcNodata) {
						highPc.setValue(row, col, true)
					}
				}
				progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Mapping Gully Depth (Loop 1):", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}

			// close the rasters and then delete the files
			String dfmeDataFile = dfme.getDataFile()
			dfme.close()
			dfme = null
			String planCurvDataFile = planCurv.getDataFile()
			planCurv.close()
			planCurv = null
			
			new File(dfmeFile).delete()
			new File(dfmeDataFile).delete()
			new File(planCurvFile).delete()
			new File(planCurvDataFile).delete()

			// create the DEM raster object
			WhiteboxRaster dem = new WhiteboxRaster(demFile, "r")
			double nodata = dem.getNoDataValue()
			double gridResX = dem.getCellSizeX()
			double gridResY = dem.getCellSizeY()
			double gridResXY = Math.sqrt(gridResX * gridResX + gridResY * gridResY)

			double backgroundVal = 0.0
			if (nodataForBackground) {
				backgroundVal = nodata
			}
			
			// create an output grid
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", demFile, DataType.FLOAT, backgroundVal)
			output.setPreferredPalette("spectrum.pal")
			output.setZUnits(dem.getZUnits())
			
			// scan each potential cross-section
			double z, stZ, endZ, cellZ, maxDepth, deltaZ, dist, depth
			double minCSLength
			boolean flag
			int j, k
			int startingCellX, startingCellY, endingCellX, endingCellY
			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
					if (lowDfme.getValue(row, col)) {
						cellZ = dem.getValue(row, col)
						boolean[] csFlags = new boolean[4]
						double[] csLengths = new double[4]
						double[] csDepths = new double[4]
						int csNum
						
						// N-S cross-section
						csNum = 0
						csLengths[csNum] = gridResY
						csFlags[csNum] = highPc.getValue(row, col)
						j = row
						k = col
						flag = true
						while (flag) {
							j--
							if (!lowDfme.getValue(j, k)) {
								startingCellY = j //+ 1
								csLengths[csNum] += gridResY
								flag = false
							} else {
								csLengths[csNum] += gridResY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}
			
						j = row
						k = col
						flag = true
						while (flag) {
							j++
							if (!lowDfme.getValue(j, k)) {
								endingCellY = j //- 1
								csLengths[csNum] += gridResY
								flag = false
							} else {
								csLengths[csNum] += gridResY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}

						if (csFlags[csNum]) {
							// The cross-section is narrow enough and contains
							// at least one cell of high plan curv. Now check
							// if it does not exceed the max gully depth.
							stZ = dem.getValue(startingCellY, col)
							endZ = dem.getValue(endingCellY, col)
							deltaZ = endZ - stZ
							maxDepth = 0
							for (int i = startingCellY; i <= endingCellY; i++) {
								dist = (i - startingCellY) * gridResY
								z = dem.getValue(i, col)
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - z
								if (depth > maxDepth) {
									maxDepth = depth 
								}
							}
							if (maxDepth <= maxGullyDepth && maxDepth > minGullyDepth) {
								dist = (row - startingCellY) * gridResY
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - cellZ
								if (depth <= 0) {
									depth = backgroundVal
								}
								csDepths[csNum] = depth
							} else {
								csFlags[csNum] = false
								csDepths[csNum] = backgroundVal
							}
						} else {
							csDepths[csNum] = backgroundVal
						}


						// E-W cross-section
						csNum = 1
						csLengths[csNum] = gridResX
						csFlags[csNum] = highPc.getValue(row, col)
						j = row
						k = col
						flag = true
						while (flag) {
							k--
							if (!lowDfme.getValue(j, k)) {
								startingCellX = k //+ 1
								csLengths[csNum] += gridResX
								flag = false
							} else {
								csLengths[csNum] += gridResX
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}
			
						j = row
						k = col
						flag = true
						while (flag) {
							k++
							if (!lowDfme.getValue(j, k)) {
								endingCellX = k //- 1
								csLengths[csNum] += gridResX
								flag = false
							} else {
								csLengths[csNum] += gridResX
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}

						if (csFlags[csNum]) {
							// The cross-section is narrow enough and contains
							// at least one cell of high plan curv. Now check
							// if it does not exceed the max gully depth.
							stZ = dem.getValue(row, startingCellX)
							endZ = dem.getValue(row, endingCellX)
							deltaZ = endZ - stZ
							maxDepth = 0
							for (int i = startingCellX; i <= endingCellX; i++) {
								dist = (i - startingCellX) * gridResX
								z = dem.getValue(row, i)
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - z
								if (depth > maxDepth) {
									maxDepth = depth 
								}
							}
							if (maxDepth <= maxGullyDepth && maxDepth > minGullyDepth) {
								dist = (col - startingCellX) * gridResX
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - cellZ
								if (depth <= 0) {
									depth = backgroundVal
								}
								csDepths[csNum] = depth
							} else {
								csFlags[csNum] = false
								csDepths[csNum] = backgroundVal
							}							
						} else {
							csDepths[csNum] = backgroundVal
						}


						// NW-SE cross-section
						csNum = 2
						csLengths[csNum] = gridResXY
						csFlags[csNum] = highPc.getValue(row, col)
						j = row
						k = col
						flag = true
						while (flag) {
							j--
							k--
							if (!lowDfme.getValue(j, k)) {
								startingCellY = j //+ 1
								startingCellX = k //+ 1
								csLengths[csNum] += gridResXY
								flag = false
							} else {
								csLengths[csNum] += gridResXY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}
			
						j = row
						k = col
						flag = true
						while (flag) {
							j++
							k++
							if (!lowDfme.getValue(j, k)) {
								endingCellY = j //- 1
								endingCellX = k //- 1
								csLengths[csNum] += gridResXY
								flag = false
							} else {
								csLengths[csNum] += gridResXY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}

						if (csFlags[csNum]) {
							// The cross-section is narrow enough and contains
							// at least one cell of high plan curv. Now check
							// if it does not exceed the max gully depth.
							stZ = dem.getValue(startingCellY, startingCellX)
							endZ = dem.getValue(endingCellY, endingCellX)
							deltaZ = endZ - stZ
							maxDepth = 0
							j = startingCellY
							k = startingCellX
							for (int i = 0; i < (int)(csLengths[csNum] / gridResXY); i++) {
								dist = Math.sqrt((j - startingCellY) * (j - startingCellY) + (k - startingCellX) * (k - startingCellX)) * gridResXY
								z = dem.getValue(j, k)
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - z
								if (depth > maxDepth) {
									maxDepth = depth 
								}
								j++
								k++
							}
							if (maxDepth <= maxGullyDepth && maxDepth > minGullyDepth) {
								dist = Math.sqrt((row - startingCellY) * (row - startingCellY) + (col - startingCellX) * (col - startingCellX)) * gridResXY
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - cellZ
								if (depth <= 0) {
									depth = backgroundVal
								}
								csDepths[csNum] = depth
							} else {
								csFlags[csNum] = false
								csDepths[csNum] = backgroundVal
							}
						} else {
							csDepths[csNum] = backgroundVal
						}

						
						// NE-SW cross-section
						csNum = 3
						csLengths[csNum] = gridResXY
						csFlags[csNum] = highPc.getValue(row, col)
						j = row
						k = col
						flag = true
						while (flag) {
							j--
							k++
							if (!lowDfme.getValue(j, k)) {
								startingCellY = j //+ 1
								startingCellX = k //- 1
								csLengths[csNum] += gridResXY
								flag = false
							} else {
								csLengths[csNum] += gridResXY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}
			
						j = row
						k = col
						flag = true
						while (flag) {
							j++
							k--
							if (!lowDfme.getValue(j, k)) {
								endingCellY = j //- 1
								endingCellX = k //+ 1
								csLengths[csNum] += gridResXY
								flag = false
							} else {
								csLengths[csNum] += gridResXY
								if (csLengths[csNum] > maxGullyWidth) {
									flag = false
									csFlags[csNum] = false
									break
								}
								csFlags[csNum] = csFlags[csNum] | highPc.getValue(j, k)
							}
						}

						if (csFlags[csNum]) {
							// The cross-section is narrow enough and contains
							// at least one cell of high plan curv. Now check
							// if it does not exceed the max gully depth.
							stZ = dem.getValue(startingCellY, startingCellX)
							endZ = dem.getValue(endingCellY, endingCellX)
							deltaZ = endZ - stZ
							maxDepth = 0
							j = startingCellY
							k = startingCellX
							for (int i = 0; i < (int)(csLengths[csNum] / gridResXY); i++) {
								dist = Math.sqrt((j - startingCellY) * (j - startingCellY) + (k - startingCellX) * (k - startingCellX)) * gridResXY
								z = dem.getValue(j, k)
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - z
								if (depth > maxDepth) {
									maxDepth = depth 
								}
								j++
								k--
							}
							if (maxDepth <= maxGullyDepth && maxDepth > minGullyDepth) {
								dist = Math.sqrt((row - startingCellY) * (row - startingCellY) + (col - startingCellX) * (col - startingCellX)) * gridResXY
								depth = (stZ + deltaZ * dist / csLengths[csNum]) - cellZ
								if (depth <= 0) {
									depth = backgroundVal
								}
								csDepths[csNum] = depth
							} else {
								csFlags[csNum] = false
								csDepths[csNum] = backgroundVal
							}
						} else {
							csDepths[csNum] = backgroundVal
						}

						// find the shortest valid cross section
						int i = -1
						minCSLength = Double.POSITIVE_INFINITY
						for (j = 0; j < 4; j++) {
							if (csFlags[j] && csLengths[j] < minCSLength) {
								minCSLength = csLengths[j]
								i = j
							}
						}
						if (i >= 0) {
							output.setValue(row, col, csDepths[i])
						}
					}
				}
				progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress("Mapping Gully Depth (Loop 2):", progress)
					oldProgress = progress
				}
				// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
			
			lowDfme = null
			planCurv = null
			dem.close()
			
			output.addMetadataEntry("Created by the " + MapGullyDepth + " tool.")
			output.addMetadataEntry("Created on " + new Date())
			String zunits = ""
			if (!output.getZUnits().toLowerCase().equals("not specified")) {
				zunits = output.getZUnits()
			}
			String xyunits = ""
			if (!output.getXYUnits().toLowerCase().equals("not specified")) {
				xyunits = output.getXYUnits()
			}
			output.addMetadataEntry("Max gully width = " + maxGullyWidth + " " + xyunits)
			output.addMetadataEntry("Max gully depth = " + maxGullyDepth + " " + zunits)
			output.addMetadataEntry("Min gully depth = " + minGullyDepth + " " + zunits)
			output.addMetadataEntry("Diff from mean elev. threshold = " + dfmeThreshold + " " + zunits)
			output.addMetadataEntry("Min cross-section plan curv = " + highPlanCurvValue)
			output.addMetadataEntry("Smoothing parameter = " + smoothingParameter)
			
			output.close()
	
			// display the output image
			pluginHost.returnData(outputFile)

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
	def f = new MapGullyDepth(pluginHost, args, name, descriptiveName)
}
