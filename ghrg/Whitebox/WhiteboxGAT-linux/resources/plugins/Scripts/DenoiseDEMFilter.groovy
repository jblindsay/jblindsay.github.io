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
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import whitebox.structures.BooleanBitArray2D
import whitebox.structures.DoubleArray2D
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "DenoiseDEM"
def descriptiveName = "De-noise DEM"
def description = "Removes short-scale variability from DEM."
def toolboxes = ["TerrainAnalysis", "LidarTools"]

public class DenoiseDEM implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public DenoiseDEM(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Input raster image file", "Input Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Output file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("Neighbourhood size", "Neighbourhood Size (cells):", "5", true, false)
            sd.addDialogCheckBox("Remove peaks only?", "Remove peaks only?", false)
            sd.addDialogCheckBox("Output a hillshade (shaded relief) raster?", "Output a hillshade raster?", false)
            sd.addDialogDataInput("Number of significant decimal places", "Significant Decimal Places:", "2", true, true)
            
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
			if (args.length < 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			int progress = 0
			int oldProgress = -1
			double z
			long sum
			int sumN, numCells
			
			Date start = new Date();
			
			// read the input parameters
			String inputFile = args[0];
			String outputFile = args[1];
			int filterSize = Integer.parseInt(args[2]);
			if (filterSize < 1) { filterSize = 1; }
			boolean peaksOnly = false;
			int numLoops = 5;
			if (args[3].toLowerCase().contains("t")) {
				peaksOnly = true;
				numLoops = 4;
			}
			boolean outputHillshade = false;
			if (args.length >= 5) {
				if (args[4].toLowerCase().contains("t")) {
					outputHillshade = true;
				}
			}
			int numSigDecimalPlaces = 3;
			if (args.length >= 6) {
				numSigDecimalPlaces = Integer.parseInt(args[5]);
				if (numSigDecimalPlaces < 0) { numSigDecimalPlaces = 0; }
				if (numSigDecimalPlaces > 8) { numSigDecimalPlaces = 8; }
			}

			int neighbourhoodSize = filterSize; //(int)(filterSize * 1.5)

			// read the input image
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()
			double minValue = image.getMinimumValue()
			double range = image.getMaximumValue() - minValue
			String demShortName = image.getShortHeaderFile()
			
//			double[][] integralImage = new double[rows][cols];
//			int[][] integralImageN = new int[rows][cols];
//
//			// calculate the integral image
//			
//			for (int row = 0; row < rows; row++) {
//				sum = 0
//				sumN = 0
//  				for (int col = 0; col < cols; col++) {
//  					z = image.getValue(row, col)
//  					if (z == nodata) {
//  						z = 0
//  					} else {
//  						z = (z - minValue) / range
//  						sumN++
//  					}
//  					sum += z
//  					if (row > 0) {
//  						integralImage[row][col] = sum + integralImage[row - 1][col]
//  						integralImageN[row][col] = sumN + integralImageN[row - 1][col]
//  					} else {
//  						integralImage[row][col] = sum
//  						integralImageN[row][col] = sumN
//  					}
//  				}
//  				progress = (int)(100f * row / rows)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Loop 1 of $numLoops:", progress)
//					oldProgress = progress
//					// check to see if the user has requested a cancellation
//					if (pluginHost.isRequestForOperationCancelSet()) {
//						pluginHost.showFeedback("Operation cancelled")
//						return
//					}
//				}
//			}

			long[][] integralImage = new long[rows][cols];
			int[][] integralImageN = new int[rows][cols];

			// calculate the integral image
//			int progress = 0
//			int oldProgress = -1
//			double z
//			long sum
//			int sumN
			double multiplier = Math.pow(10, numSigDecimalPlaces);
			
			for (int row = 0; row < rows; row++) {
				sum = 0
				sumN = 0
  				for (int col = 0; col < cols; col++) {
  					z = image.getValue(row, col)
  					if (z == nodata) {
  						z = 0
  					} else {
  						z = (z - minValue) * multiplier
  						sumN++
  					}
  					sum += (long)(Math.round(z))
  					if (row > 0) {
  						integralImage[row][col] = sum + integralImage[row - 1][col]
  						integralImageN[row][col] = sumN + integralImageN[row - 1][col]
  					} else {
  						integralImage[row][col] = sum
  						integralImageN[row][col] = sumN
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 1 of 2:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	  inputFile, DataType.FLOAT, nodata)
			output.setPreferredPalette(image.getPreferredPalette()) //"blue_white_red.plt")
			output.setForceAllDataInMemory(true)

			// calculate the detrended DEM
			//double[][] detrendedDEM = new double[rows][cols];
			DoubleArray2D detrendedDEM = new DoubleArray2D(rows, cols, nodata)
			
			oldProgress = -1
			int x1, x2, y1, y2, x, y
			double a, b, c, d, outValue
			int numNoDataVals
			for (int row = 0; row < rows; row++) {
				y1 = row - neighbourhoodSize
				if (y1 < 0) { y1 = 0 }
				if (y1 >= rows) { y1 = rows - 1 }

				y2 = row + neighbourhoodSize
				if (y2 < 0) { y2 = 0 }
				if (y2 >= rows) { y2 = rows - 1 }
				
				for (int col = 0; col < cols; col++) {
					z = image.getValue(row, col)

  					if (z != nodata) {
  						output.setValue(row, col, z)
  						x1 = col - neighbourhoodSize
						if (x1 < 0) { x1 = 0 }
						if (x1 >= cols) { x1 = cols - 1 }
	
						x2 = col + neighbourhoodSize
						if (x2 < 0) { x2 = 0 }
						if (x2 >= cols) { x2 = cols - 1 }

						a = integralImage[y1][x1]
						b = integralImage[y1][x2]
						c = integralImage[y2][x2]
						d = integralImage[y2][x1]

						numCells = (int)(integralImageN[y2][x2] + integralImageN[y1][x1] - integralImageN[y1][x2] - integralImageN[y2][x1])

						outValue = z - (((c + a - b - d) / numCells) / multiplier + minValue) //z - ((c + a - b - d) / numCells * range + minValue)
						detrendedDEM.setValue(row, col, outValue)
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 2 of $numLoops:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			integralImage = new double[0][0];
			integralImageN = new int[0][0];

			// find peaks and pits
			BooleanBitArray2D isPeak = new BooleanBitArray2D(rows, cols)
			BooleanBitArray2D isPit = new BooleanBitArray2D(rows, cols)
			int[] dX = [ 1, 1, 1, 0, -1, -1, -1, 0 ]
			int[] dY = [ -1, 0, 1, 1, 1, 0, -1, -1 ]
			int n, rowN, colN
			double zN, diff
			boolean isHighest, isLowest
			int windowDiameter = filterSize * 2 + 1
			int windowDiameterLessOne = windowDiameter - 1
			int numWindowCells = windowDiameter * windowDiameter
			GridCell gc

			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
  					z = detrendedDEM.getValue(row, col)
  					if (z != nodata) {
	  					isHighest = true
	  					isLowest = true
	  					for (n = 0; n < 8; n++) {
	  						rowN = row + dY[n];
	                    	colN = col + dX[n];
	                    	if (rowN >= 0 && rowN < rows && colN >= 0 && colN < cols) {
								zN = detrendedDEM.getValue(rowN, colN)
								if (zN != nodata) {
									if (zN > z) {
										isHighest = false;
									}
									if (zN < z) {
										isLowest = false;
									}
								}
	                    	}
	  					}
	  					if (isHighest) {
							isPeak.setValue(row, col, true)
						}
						if (isLowest) {
							isPit.setValue(row, col, true)
						}
  					}
				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 3 of $numLoops:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			// find peaks
			oldProgress = -1
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
  					z = detrendedDEM.getValue(row, col)
  					if (z != nodata) {
						if (isPeak.getValue(row, col)) {
	  						// create a subgrid of the neighbouring area
	  						double lowValue = Double.NEGATIVE_INFINITY
							DoubleArray2D data = new DoubleArray2D(windowDiameter, windowDiameter, nodata)
	  						x1 = col - filterSize
							x2 = col + filterSize
							y1 = row - filterSize
							y2 = row + filterSize
							
							for (x = x1; x <= x2; x++) {
								for (y = y1; y <= y2; y++) {
									zN = detrendedDEM.getValue(y, x)
									if (zN != nodata) {
										data.setValue(y - y1, x - x1, -zN)
									} else {
										data.setValue(y - y1, x - x1, nodata)
									}
								}
							}

							// remove depressions (equivalent to cleaving the peaks, since the surface has been inverted)
							PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>(numWindowCells);
							BooleanBitArray2D inQueue = new BooleanBitArray2D(windowDiameter, windowDiameter)

							for (y = 0; y < windowDiameter; y++) {
								for (x = 0; x < windowDiameter; x++) {
									z = data.getValue(y, x)
									if (z != nodata) {
										for (n = 0; n < 8; n++) {
					  						rowN = y + dY[n];
					                    	colN = x + dX[n];
					                    	zN = data.getValue(rowN, colN)
					                    	if (zN == nodata) {
					                    		queue.add(new GridCell(y, x, z))
					                    		inQueue.setValue(y, x, true)
					                    	}
										}
									} else {
										inQueue.setValue(y, x, true)
									}
								}
							}

							while (!queue.isEmpty()) {
				                gc = queue.poll();
				                y = gc.row;
				                x = gc.col;
				                z = gc.z;
				                
				                for (int i = 0; i < 8; i++) {
				                    rowN = y + dY[i];
				                    colN = x + dX[i];
				                    zN = data.getValue(rowN, colN);
				                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
				                        if (zN < z) {
				                        	diff = z - zN
				                        	zN = z
				                        	output.decrementValue(y1 + rowN, x1 + colN, diff)
				                        	detrendedDEM.setValue(y1 + rowN, x1 + colN, -zN)
				                        	data.setValue(rowN, colN, zN)
				                        }
				                        gc = new GridCell(rowN, colN, zN);
				                        queue.add(gc);
				                        inQueue.setValue(rowN, colN, true)
				                    }
				                }
	            			}
	  					}
  					}
  				}
  				progress = (int)(100f * row / rows)
				if (progress > oldProgress) {
					pluginHost.updateProgress("Loop 4 of $numLoops:", progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			if (!peaksOnly) {
				// find sinks
				oldProgress = -1
				for (int row = 0; row < rows; row++) {
					for (int col = 0; col < cols; col++) {
	  					z = detrendedDEM.getValue(row, col)
	  					if (z != nodata) {
							if (isPit.getValue(row, col)) {
		  						// create a subgrid of the neighbouring area
		  						double lowValue = Double.NEGATIVE_INFINITY
								DoubleArray2D data = new DoubleArray2D(windowDiameter, windowDiameter, nodata)
		  						x1 = col - filterSize
								x2 = col + filterSize
								y1 = row - filterSize
								y2 = row + filterSize
								
								for (x = x1; x <= x2; x++) {
									for (y = y1; y <= y2; y++) {
										zN = detrendedDEM.getValue(y, x)
										data.setValue(y - y1, x - x1, zN)
									}
								}
	
								// remove depressions (equivalent to cleaving the peaks, since the surface has been inverted)
								PriorityQueue<GridCell> queue = new PriorityQueue<GridCell>(numWindowCells);
								BooleanBitArray2D inQueue = new BooleanBitArray2D(windowDiameter, windowDiameter)
	
								for (y = 0; y < windowDiameter; y++) {
									for (x = 0; x < windowDiameter; x++) {
										z = data.getValue(y, x)
										if (z != nodata) {
											for (n = 0; n < 8; n++) {
						  						rowN = y + dY[n];
						                    	colN = x + dX[n];
						                    	zN = data.getValue(rowN, colN)
						                    	if (zN == nodata) {
						                    		queue.add(new GridCell(y, x, z))
						                    		inQueue.setValue(y, x, true)
						                    	}
											}
										} else {
											inQueue.setValue(y, x, true)
										}
									}
								}
	
								while (!queue.isEmpty()) {
					                gc = queue.poll();
					                y = gc.row;
					                x = gc.col;
					                z = gc.z;
					                
					                for (int i = 0; i < 8; i++) {
					                    rowN = y + dY[i];
					                    colN = x + dX[i];
					                    zN = data.getValue(rowN, colN);
					                    if ((zN != nodata) && (!inQueue.getValue(rowN, colN))) {
					                        if (zN < z) {
					                        	diff = z - zN
					                        	zN = z
					                        	output.incrementValue(y1 + rowN, x1 + colN, diff)
					                        	detrendedDEM.setValue(y1 + rowN, x1 + colN, zN)
					                        	data.setValue(rowN, colN, zN)
					                        }
					                        gc = new GridCell(rowN, colN, zN);
					                        queue.add(gc);
					                        inQueue.setValue(rowN, colN, true)
					                    }
					                }
		            			}
		  					}
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress > oldProgress) {
						pluginHost.updateProgress("Loop 5 of 5:", progress)
						oldProgress = progress
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}
			}

			output.flush()
			output.findMinAndMaxVals()
			output.setDisplayMaximum(image.getDisplayMaximum())
			output.setDisplayMinimum(image.getDisplayMinimum())
			output.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			Date stop = new Date()
			TimeDuration td = TimeCategory.minus(stop, start)
			output.addMetadataEntry("Elapsed time: $td")
			output.addMetadataEntry("Input DEM: $demShortName")
			output.addMetadataEntry("Filter Half-size: $filterSize")

			if (outputHillshade) {
				double outNoData = -32768;
				WhiteboxRaster hillshade = new WhiteboxRaster(outputFile.replace(".dep", "_hillshade.dep"), "rw", 
  		  	  		inputFile, DataType.INTEGER, outNoData)
				hillshade.setPreferredPalette("grey.pal")
				hillshade.setNoDataValue(outNoData);
            
				final double radToDeg = 180 / Math.PI;
		        final double degToRad = Math.PI / 180;
		        double azimuth = (315 - 90) * degToRad;
		        double altitude = 30 * degToRad;
		        double sinTheta = Math.sin(altitude);
		        double cosTheta = Math.cos(altitude);    
                double tanSlope;
		        double fx, fy, aspect;
		        double gridRes, eightGridRes;
		        double[] N = new double[8];
		        double term1, term2, term3;
		        
		        gridRes = image.getCellSizeX();
            	eightGridRes = 8 * gridRes;
            
				long[] histo = new long[256];
            	numCells = 0;

            	oldProgress = -1
				for (int row = 0; row < rows; row++) {
					for (int col = 0; col < cols; col++) {
	  					z = output.getValue(row, col)
	  					if (z != nodata) {
	  						// get the neighbouring cell Z values
	                        for (int i = 0; i < 8; i++) {
	                            N[i] = output.getValue(row + dY[i], col + dX[i]);
	                            if (N[i] == nodata) {
	                               N[i] = z;
	                            }
	                        }
	                        // calculate slope and aspect
	                        fy = (N[6] - N[4] + 2 * (N[7] - N[3]) + N[0] - N[2]) / eightGridRes;
	                        fx = (N[2] - N[4] + 2 * (N[1] - N[5]) + N[0] - N[6]) / eightGridRes;
	                        if (fx != 0) {
	                            tanSlope = Math.sqrt(fx * fx + fy * fy);
	                            aspect = (180 - Math.atan(fy / fx) * radToDeg + 90 * (fx / Math.abs(fx))) * degToRad;
	                            term1 = tanSlope / Math.sqrt(1 + tanSlope * tanSlope);
	                            term2 = sinTheta / tanSlope;
	                            term3 = cosTheta * Math.sin(azimuth - aspect);
	                            z = term1 * (term2 - term3);
	                        } else {
	                            z = 0.5;
	                        }
	                        z = (int)(z * 255);
	                        if (z < 0) {
	                            z = 0;
	                        }
	                        histo[(int) z]++;
	                        numCells++;
	                        hillshade.setValue(row, col, z);
	  					}
	  				}
	  				progress = (int)(100f * row / rows)
					if (progress > oldProgress) {
						pluginHost.updateProgress("Calculating hillshade:", progress)
						oldProgress = progress
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}

				// trim the display min and max values by clipPercent
	            double clipPercent = 0.01;
	            if (args.length >= 6) {
	                clipPercent = Double.parseDouble(args[5]) / 100.0;
	            }
	            int newMin = 0;
	            int newMax = 0;
	            double targetCellNum = numCells * clipPercent;
	            sum = 0;
	            for (int i = 0; i < 256; i++) {
	                sum += histo[i];
	                if (sum >= targetCellNum) {
	                    newMin = i;
	                    break;
	                }
	            }
	
	            sum = 0;
	            for (int i = 255; i >= 0; i--) {
	                sum += histo[i];
	                if (sum >= targetCellNum) {
	                    newMax = i;
	                    break;
	                }
	            }

	            hillshade.flush()
				hillshade.findMinAndMaxVals()
			
	            if (newMax > newMin) {
	                hillshade.setDisplayMinimum((double) newMin);
	                hillshade.setDisplayMaximum((double) newMax);
	            }
				
	            hillshade.addMetadataEntry("Created by the " + descriptiveName + " tool.")
	        	hillshade.addMetadataEntry("Created on " + new Date());
	
	            hillshade.close();

	            pluginHost.returnData(outputFile.replace(".dep", "_hillshade.dep"))
			}
			
			output.close()
			image.close()
	
	
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

    @CompileStatic
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        public double z;
        
        public GridCell(int row, int col, double z) {
            this.row = row;
            this.col = col;
            this.z = z;
        }

        @Override
        public int compareTo(GridCell other) {
    	  	if (this.z > other.z) {
        		return 1
        	} else if (this.z < other.z) {
        		return -1
        	} else {
        		return 0
        	}
        }
    }
}

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def f = new DenoiseDEM(pluginHost, args, name, descriptiveName)
}
