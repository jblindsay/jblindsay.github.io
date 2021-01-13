/*
 * Copyright (C) 2016 Dr. John Lindsay <jlindsay@uoguelph.ca>
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

import javax.swing.*
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.util.PriorityQueue
import java.util.Date
import java.util.ArrayList
import java.util.Queue
import java.util.LinkedList
import java.util.ArrayDeque
import java.text.DecimalFormat
import java.util.stream.IntStream
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.structures.BooleanBitArray2D
import whitebox.structures.DoubleArray2D
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ImpoundmentIndex"
def descriptiveName = "Impoundment Index"
def description = "Calculates the impoundment size resulting from damming a DEM."
def toolboxes = ["FlowpathTAs", "WetlandTools"]

public class ImpoundmentIndex implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName

	private AtomicInteger numSolved = new AtomicInteger(0)
	private AtomicInteger curProgress = new AtomicInteger(-1)
	
	public ImpoundmentIndex(WhiteboxPluginHost pluginHost, 
		String[] args, String name, String descriptiveName) {
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
			sd.addDialogFile("Input DEM raster", "Input DEM Raster:", "open", "Raster Files (*.dep), DEP", true, false)
            sd.addDialogFile("Output impoundment index raster file", "Output Impoundment Size Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogComboBox("What measure of impoundment size should be output?", "Output impoundment", ["area", "volume"], 0)
			sd.addDialogDataInput("Max. dam height (in z-units)", "Max. dam height (in z-units)", "", true, false)
            sd.addDialogDataInput("Max. dam length (in grid cells)", "Max. dam length (in grid cells)", "", true, false)
            
//			DialogCheckBox outDEM = sd.addDialogCheckBox("Output dammed DEM?", "Output dammed DEM?", false)
//			DialogDataInput lowerThreshold = sd.addDialogDataInput("Lower threshold", "Lower threshold", "", true, false)
//            DialogDataInput upperThreshold = sd.addDialogDataInput("Upper threshold", "Upper threshold", "", true, false)
//			lowerThreshold.visible = false
//			upperThreshold.visible = false
//							
//            //Listener for chxError            
//            def lstr = { evt -> if (evt.getPropertyName().equals("value")) 
//            	{ 
//            		String value = outDEM.getValue()
//            		if (!value.isEmpty() && value != null) {
//            			if (outDEM.getValue() == "true") {
//            				lowerThreshold.visible = true
//							upperThreshold.visible = true
//		            	} else {
//		            		lowerThreshold.visible = false
//							upperThreshold.visible = false
//		            	}
//            		}
//            	} 
//            } as PropertyChangeListener
//            outDEM.addPropertyChangeListener(lstr)
//
			
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

			//Date start = new Date()
			
			int progress, oldProgress, i, row, col, rN, cN, rN2, cN2, dir, perpDir1, perpDir2, n

			/*
			 * 6  7  0
			 * 5  X  1
			 * 4  3  2
			 */
	  		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]

			int[] perpendicular1 = [ 2, 3, 4, 5, 6, 7, 0, 1 ]
			int[] perpendicular2 = [ 6, 7, 0, 1, 2, 3, 4, 5 ]
			
			int[] inflowingVals = [ 4, 5, 6, 7, 0, 1, 2, 3 ]
        	double z, zN, height, outVal
        	DecimalFormat df = new DecimalFormat("###.#")
			
			if (args.length < 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			// read the input parameters
			String inputFile = args[0]
			String outputFile = args[1]
			int mode = 1 // 0 = output volume; 1 = output area
			if (args[2].toLowerCase().contains("vol")) {
				mode = 0
			}
			double maxDamHeight = Double.parseDouble(args[3])
			int damLength = Integer.parseInt(args[4])
			damLength = (int)(Math.floor(damLength / 2))
			boolean createDEM = false

			int damProfileLength = damLength * 2 + 1
			double[] damProfile = new double[damProfileLength]
			double[] damProfileFilled = new double[damProfileLength]
							
			// read the input image
			WhiteboxRaster dem = new WhiteboxRaster(inputFile, "r")
			dem.setForceAllDataInMemory(true)
			
			double nodata = dem.getNoDataValue()
			int rows = dem.getNumberRows()
			int rowsLessOne = rows - 1
			int cols = dem.getNumberColumns()
			int numPixels = rows * cols
			double gridResX = dem.getCellSizeX();
            double gridResY = dem.getCellSizeY();
            if (dem.getXYUnits().toLowerCase().contains("deg")) {
            	// estimate the cell size 
                double p1 = 111412.84;		// longitude calculation term 1
                double p2 = -93.5;			// longitude calculation term 2
                double p3 = 0.118;			// longitude calculation term 3
                double lat = Math.toRadians((dem.getNorth() - dem.getSouth()) / 2.0);
                double longlen = (p1 * Math.cos(lat)) + (p2 * Math.cos(3 * lat)) + (p3 * Math.cos(5 * lat));

                gridResX = gridResX * longlen;
                gridResY = gridResY * longlen;
            }
            double diagGridRes = Math.sqrt(gridResX * gridResX + gridResY * gridResY);
            double[] gridLengths = [diagGridRes, gridResX, diagGridRes, gridResY, diagGridRes, gridResX, diagGridRes, gridResY]
			double cellArea = gridResX * gridResY

			DoubleArray2D maxDamHeightArray = new DoubleArray2D(rows, cols, nodata)

			/*
			 * It doesn't make much sense to have users input a 
			 * hydrologically conditioned, depressionless DEM to a 
			 * tool that is effectively adding depressions, i.e. 
			 * dams, into the DEM. So instead of the usual way 
			 * of calcuating flow directions from the DEM, I have 
			 * opted for using a flood-order operation, which can
			 * calculate flow directions (back-links) for DEMs that
			 * contain depressions and other sinks. Flood order 
			 * calculation is based on a priority-flood algorithm and
			 * requires two steps; first an initialization of the 
			 * grids and second the actual priority-flood operation.
			 */
			
			byte[] backLink = [4, 5, 6, 7, 0, 1, 2, 3]
			byte[][] flowDir = new byte[rows][cols]
			boolean isPit, isEdgeCell
			GridCell gc
			PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>((2 * rows + 2 * cols) * 2);
			BooleanBitArray2D inQueue = new BooleanBitArray2D(rows, cols)
			int numSolvedCells = 0
			int numCellsTotal = rows * cols
			int colN, rowN

			int numValidCells = 0
			
			// initialize the grids
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col)
					maxDamHeightArray.setValue(row, col, z)
					flowDir[row][col] = -2
					if (z != nodata) {
						isPit = true
						isEdgeCell = false
						for (n = 0; n < 8; n++) {
							zN = dem.getValue(row + dY[n], col + dX[n])
							if (isPit && zN != nodata && zN < z) {
								isPit = false
							} else if (zN == nodata) {
								isEdgeCell = true
							}
						}
						if (isPit && isEdgeCell) {
							queue.add(new GridCell(row, col, z, -1))
							flowDir[row][col] = -1
							inQueue.setValue(row, col, true)
						}
						numValidCells++
					} else {
                        numSolvedCells++
                    }
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Calculating flow patterns (Loop 1 of 2)", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			oldProgress = (int) (100f * numSolvedCells / numCellsTotal);
            
			int[] floodOrderRow = new int[numValidCells]
			int[] floodOrderCol = new int[numValidCells]
			int cellOrderVal = 0
			
			int flatIndex
			int k
			while (queue.isEmpty() == false) {
                gc = queue.poll();
                row = gc.row;
                col = gc.col;
                flatIndex = gc.flatIndex;
                
				floodOrderRow[cellOrderVal] = row
                floodOrderCol[cellOrderVal] = col
                cellOrderVal++
                
                for (n = 0; n < 8; n++) {
                    rowN = row + dY[n];
                    colN = col + dX[n];
                    zN = dem.getValue(rowN, colN)
                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
                        numSolvedCells++;
                        flowDir[rowN][colN] = backLink[n]
                        k = 0
						if (zN == dem.getValue(row, col)) {
							k = flatIndex + 1
						}
                        queue.add(new GridCell(rowN, colN, zN, k))
						
                        inQueue.setValue(rowN, colN, true)
                    }
                }
                progress = (int) (100f * numSolvedCells / numCellsTotal);
                if (progress > oldProgress) {
                    pluginHost.updateProgress("Calculating flow patterns (Loop 2 of 2)", progress)
                    oldProgress = progress;
                    // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
                }
            }

			// count the number of inflowing neighbours
			byte[][] numInflow = new byte[rows][cols]
			numSolved.set(-1)
			curProgress.set(0)
			IntStream.range(0, rows).parallel().forEach{ countInflowingCells(it, flowDir, numInflow) };

			/*
			 * If we were using a hydrologically conditioned DEM, then
			 * we could guarantee that all downslope cells are of lower
			 * elevation (i.e. monotonic downward slopes). All cutoff  
			 * values of height used to parse the elevation distribution 
			 * of upstream sites would remove values that cannot be 
			 * useful for downstream sites. However, we cannot make 
			 * this assumption with non-conditioned DEMs and it is
			 * likely that there will be downstream cells that are of
			 * higher elevation. This would be a problem if we used
			 * z + maxDamHeight as the cutoff value for clipping
			 * values from the distribution since z upstream may be 
			 * lower than the current z, thus removing potentially 
			 * useful values from the distribution. To solve this  
			 * we need to calculate the height cutoff values in 
			 * the flood order such that it only ever increases 
			 * upstream. This is effectively the same as filling 
			 * the DEM. 
			 * 
			 * Rather than performing a second priority-flood
			 * operation (the first being used to calculate flow 
			 * directions above) we simply saved the order of cells 
			 * during the first operation. This requires more memory 
			 * but is much faster than the original priority-flood op.
			 */
			double[][] cutoffHeights = new double[rows][cols]
			oldProgress = -1
			for (i = 0; i < numValidCells; i++) {
				row = floodOrderRow[i]
				col = floodOrderCol[i]
				z = dem.getValue(row, col)
				dir = flowDir[row][col]
				if (dir >= 0) {
					rN = row + dY[dir]
					cN = col + dX[dir]
					if (z + maxDamHeight > cutoffHeights[rN][cN]) {
						cutoffHeights[row][col] = z + maxDamHeight
					} else {
						cutoffHeights[row][col] = cutoffHeights[rN][cN]
					}
				} else {
					cutoffHeights[row][col] = z + maxDamHeight
				}
				progress = (int) (100f * i / (numValidCells - 1));
                if (progress != oldProgress) {
                    pluginHost.updateProgress("Calculating Cutoff Heights:", progress)
                    oldProgress = progress;
                    // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
                }
			}

			// initialize the output image
			pluginHost.updateProgress("Creating output file:", 0)
			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFile, DataType.FLOAT, 0.0)
			output.setPreferredPalette("spectrum_white_background.pal");
            output.setDataScale(DataScale.CONTINUOUS);
            output.setZUnits("upslope grid cells");
            output.setNonlinearity(0.2)

			double[][] damHeight = new double[rows][cols]
			UpslopeValues[][] upslopeVals = new UpslopeValues[rows][cols]
			oldProgress = -1
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					z = dem.getValue(row, col);
					if (z != nodata) {
						upslopeVals[row][col] = new UpslopeValues(z)
						// what's the flow direction?
						dir = flowDir[row][col]

						// what's the perpendicular flow direction?
						perpDir1 = perpendicular1[dir]
						perpDir2 = perpendicular2[dir]

						damProfile = new double[damProfileLength]

						damProfile[damLength] = z

						// find the dam height
						rN = row
						cN = col
						rN2 = row
						cN2 = col
						for (i = 1; i <= damLength; i++) {
							rN += dY[perpDir1]
							cN += dX[perpDir1]
							zN = dem.getValue(rN, cN)
							if (zN != nodata) {
								damProfile[damLength + i] = zN
							} else {
								damProfile[damLength + i] = -100000.0
							}

							rN2 += dY[perpDir2]
							cN2 += dX[perpDir2]
							zN = dem.getValue(rN2, cN2)
							if (zN != nodata) {
								damProfile[damLength - i] = zN
							} else {
								damProfile[damLength - i] = -100000.0
							}
						}

						damProfileFilled = new double[damProfileLength]
						damProfileFilled[0] = damProfile[0]
						for (i = 1; i < damProfileLength; i++) {
							if (damProfile[i] >= damProfileFilled[i - 1]) {
								damProfileFilled[i] = damProfile[i]
							} else {
								damProfileFilled[i] = damProfileFilled[i - 1]
							}
							// check to see it isn't any more than the max dam height
							if (damProfileFilled[i] - damProfile[i] > maxDamHeight) {
								damProfileFilled[i] = damProfile[i] + maxDamHeight
							}
						}

						damProfileFilled[damProfileLength - 1] = damProfile[damProfileLength - 1]
						for (i = damProfileLength - 2; i >= 0; i--) {
							if (damProfile[i] < damProfileFilled[i]) {
								if (damProfile[i] > damProfileFilled[i + 1]) {
									damProfileFilled[i] = damProfile[i]
								} else if (damProfileFilled[i] > damProfileFilled[i + 1]) {
									damProfileFilled[i] = damProfileFilled[i + 1]
								}
							}
						}
						
						damHeight[row][col] = damProfileFilled[damLength] - damProfile[damLength]


						z += damHeight[row][col]

						zN = maxDamHeightArray.getValue(row, col)
						if (zN != nodata && z > zN) {
							maxDamHeightArray.setValue(row, col, z)
						}

						rN = row
						cN = col
						rN2 = row
						cN2 = col
						for (i = 1; i <= damLength; i++) {
							rN += dY[perpDir1]
							cN += dX[perpDir1]
							zN = maxDamHeightArray.getValue(rN, cN)
							if (zN != nodata && z > zN) {
								maxDamHeightArray.setValue(rN, cN, z)
							}

							rN2 += dY[perpDir2]
							cN2 += dX[perpDir2]
							zN = maxDamHeightArray.getValue(rN2, cN2)
							if (zN != nodata && z > zN) {
								maxDamHeightArray.setValue(rN2, cN2, z)
							}
						}
					} else {
						output.setValue(row, col, nodata)
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Calculating Index (1 of 2):", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

//			for (row = 0; row < rows; row++) {
//				for (col = 0; col < cols; col++) {
//					output.setValue(row, col, maxDamHeightArray.getValue(row, col));
//				}
//				progress = (int)(100f * row / rowsLessOne)
//				if (progress != oldProgress) {
//					pluginHost.updateProgress("Calculating Index (2 of 2):", progress)
//					oldProgress = progress
//
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}
//
//			output.addMetadataEntry("Created by the "
//	                    + descriptiveName + " tool.")
//	        output.addMetadataEntry("Created on " + new Date())
//			output.close()
//
//			return;

			
			// close the DEM, it won't be needed any further.
			dem.close()

			// Now perform the accumulation operation.
			int r, c
			boolean flag
			ArrayList<Double> vals
			for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					if (numInflow[row][col] == 0) {
						r = row
						c = col
						numInflow[r][c] = -1
						flag = true
						while (flag) {
							height = damHeight[r][c]
							vals = upslopeVals[r][c].getValues()
							if (mode == 0) { // volume
								outVal = upslopeVals[r][c].volumeBelow(height, cellArea)
							} else { // area
								outVal = upslopeVals[r][c].numLessThan(height) * cellArea
							}
							upslopeVals[r][c].values = null
							upslopeVals[r][c] = null
							output.setValue(r, c, outVal)
							
							dir = flowDir[r][c]
							if (dir >= 0) {
								rN = r + dY[dir]
								cN = c + dX[dir]
								if (flowDir[rN][cN] >= 0) {
									upslopeVals[rN][cN].addValues(vals, cutoffHeights[rN][cN])
									numInflow[rN][cN] -= 1
									if (numInflow[rN][cN] > 0) {
										flag = false
									} else {
										numInflow[rN][cN] = -1
										r = rN
										c = cN
									}
								} else if (flowDir[rN][cN] == -1) {
									flag = false
								}
							} else {
								flag = false
							}
						}
					}
					if (flowDir[row][col] == -2) {
						output.setValue(row, col, nodata)
					}
				}
				progress = (int)(100f * row / rowsLessOne)
				if (progress != oldProgress) {
					pluginHost.updateProgress("Calculating Index (2 of 2):", progress)
					oldProgress = progress

					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

//			for (row = 0; row < rows; row++) {
//					for (col = 0; col < cols; col++) {
//						output.setValue(row, col, damHeight[row][col])
//					}
//			}
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
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

//	private void calcFlowDirection(int row, double[][] data, byte[][] flowDir, double[] gridLengths, double nodata) {
//		// check to see if the user has requested a cancellation
//		if (pluginHost.isRequestForOperationCancelSet()) {
//			return
//		}
//		
//		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
//		int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
//		int cols = flowDir[0].length
//        int rows = flowDir.length
//        int dir, cN
//		double z, zN, maxSlope, slope
//		for (int col = 0; col < cols; col++) {
//			z = data[1][col]
//			if (z != nodata) {
//				dir = -1
//				maxSlope = -99999999
//				for (int i = 0; i < 8; i++) {
//					cN = col + dX[i]
//					if (cN >= 0 && cN < cols) {
//						zN = data[1 + dY[i]][cN]
//						if (zN != nodata) {
//							slope = (z - zN) / gridLengths[i]
//							if (slope > maxSlope && slope > 0) {
//								maxSlope = slope
//								dir = i
//							}
//						}
//					}
//				}
//				flowDir[row][col] = (byte)dir
//			} else {
//				flowDir[row][col] = (byte)-2
//			}
//			
//		}
//
//		int solved = numSolved.incrementAndGet()
//		int progress = (int) (100f * solved / (rows - 1))
//		if (progress > curProgress.get()) {
//			curProgress.incrementAndGet()
//			pluginHost.updateProgress("Calculating Flow Directions:", progress)
//		}
//	}

	private void countInflowingCells(int row, byte[][] flowDir, byte[][] numInflow) {
		// check to see if the user has requested a cancellation
		if (pluginHost.isRequestForOperationCancelSet()) {
			return
		}
		
		int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
		int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
		int[] inflowingVals = [ 4, 5, 6, 7, 0, 1, 2, 3 ]
        int cols = flowDir[0].length
        int rows = flowDir.length
        int cN, rN, n
		for (int col = 0; col < cols; col++) {
			if (flowDir[row][col] >= 0) {
				n = 0
				for (int i = 0; i < 8; i++) {
					rN = row + dY[i]
					cN = col + dX[i]
					if (rN >= 0 && cN >= 0 && rN < rows && cN < cols) {
						if (flowDir[rN][cN] == inflowingVals[i]) {
							n++
						}
					}
				}
				numInflow[row][col] = (byte)n
			} else {
				//output.setValue(row, col, nodata)
				numInflow[row][col] = (byte)-1
			}
		}

		int solved = numSolved.incrementAndGet()
		int progress = (int) (100f * solved / (rows - 1))
		if (progress > curProgress.get()) {
			curProgress.incrementAndGet()
			pluginHost.updateProgress("Counting Inflowing Neighbours:", progress)
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

	@CompileStatic
    class UpslopeValues {
		ArrayList<Double> values
		double z
    	
    	UpslopeValues(double z) {
    		this.z = z
    	}

    	ArrayList<Double> getValues() {
    		if (this.values == null) {
    			values = new ArrayList<Double>()
    			values.add(z)
    		}
    		return values
    	}

    	void addValues(ArrayList<Double> otherValues, double cutoff) {
    		if (this.values == null) {
    			values = new ArrayList<Double>()
    			values.add(z)
    		}
    		/* Any value greater than cutoff is no longer
    		   a candidate since all downslope cells will 
    		   be lower than z and no cells higher than 
    		   the maxDamHeight value will be included. Note 
    		   the caveat that we are using DEMs containing
    		   depressiosn and therefore it is possible for 
    		   downstream locations to be higher. This fact 
    		   is accounted for in the calculation of the 
    		   cutoff values.
    		*/
    		for (Double zN : otherValues) {
    			if (zN <= cutoff) {
    				values.add(zN)
    			}
    		}
    	}

    	int numLessThan(double height) {
    		if (this.values == null) {
    			values = new ArrayList<Double>()
    			values.add(z)
    		}
    		int numLess = 0
    		if (height > 0) {
    			double cutoff = z + height
	    		for (double z2 : values) {
	    			if (z2 <= cutoff) { 
	    				numLess++
	    			}
	    		}
    		}
    		return numLess
    	}

    	double volumeBelow(double height, double areaMultiplier) {
    		if (this.values == null) {
    			values = new ArrayList<Double>()
    			values.add(z)
    		}
    		double volume = 0
    		if (height > 0) {
    			double cutoff = z + height
	    		for (double z2 : values) {
	    			if (z2 <= cutoff) { 
	    				volume += (cutoff - z2)
	    			}
	    		}
    		}
    		return volume * areaMultiplier
    	}
    }

    @CompileStatic
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;
        public int flatIndex;

        public GridCell(int row, int col, double z, int flatIndex) {
            this.row = row;
            this.col = col;
            this.z = z;
            this.flatIndex = flatIndex;
        }

        @Override
        public int compareTo(GridCell other) {
        	if (this.z > other.z) {
        		return 1
        	} else if (this.z < other.z) {
        		return -1
        	} else {
        		if (this.flatIndex > other.flatIndex) {
        			return 1
        		} else if (this.flatIndex < other.flatIndex) {
        			return -1
        		}
				return 0
        	}
    	}
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def di = new ImpoundmentIndex(pluginHost, args, name, descriptiveName)
}
