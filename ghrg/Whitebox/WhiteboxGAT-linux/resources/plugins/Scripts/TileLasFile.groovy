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
import java.util.concurrent.Future
import java.util.concurrent.*
import java.text.DecimalFormat
import java.util.Date
import java.util.GregorianCalendar
import java.util.ArrayList
import java.util.Arrays
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.LASReader
import whitebox.geospatialfiles.LASReader.PointRecord
import whitebox.geospatialfiles.LASReader.PointRecColours
//import whitebox.geospatialfiles.shapefile.PointZ
import whitebox.structures.KdTree
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.geospatialfiles.LASReader.VariableLengthRecord
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.geospatialfiles.shapefile.attributes.*
import groovy.transform.CompileStatic
import whitebox.structures.BoundingBox

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "TileLasFile"
def descriptiveName = "Tile LAS File"
def description = "Breaks a LAS file into a series of tiles."
def toolboxes = ["LidarTools"]

public class TileLasFile implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public TileLasFile(WhiteboxPluginHost pluginHost, 
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
			sd.addDialogFile("Select the input LAS file", "Input LAS File:", "open", "LAS Files (*.las), LAS", true, false)
			sd.addDialogDataInput("Enter a tile size in the x dimension", "x-Dimension Tile Size:", "", true, false)
            sd.addDialogDataInput("Enter a tile size in the y dimension", "y-Dimension Tile Size:", "", true, false)
            sd.addDialogDataInput("Enter grid origin x-coordinate", "Origin X-Coordinate:", "0.0", true, false)
            sd.addDialogDataInput("Enter grid origin y-coordinate", "Origin Y-Coordinate", "0.0", true, false)
            sd.addDialogDataInput("Enter the minimum number of points in a tile (Optional).", "Min. Number of Points:", "0", true, true)
            
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

			String lasFile = args[0]
            double widthX = Double.parseDouble(args[1])
            double widthY = Double.parseDouble(args[2])
            double originX = Double.parseDouble(args[3])
            double originY = Double.parseDouble(args[4])
            int minPoints = 0
            if (args.length >= 6) {
            	if (args[5] != null && !args[5].isEmpty() &&
            	     !args[5].toLowerCase().equals("not specified")) {
            		minPoints = Integer.parseInt(args[5])
            	}
            }
            
			LASReader las = new LASReader(lasFile)
			int totalPoints = (int)las.getNumPointRecords()

			BoundingBox extent = new BoundingBox(las.getMinX(), las.getMinY(), las.getMaxX(), las.getMaxY());

			int startXGrid = (int)(Math.floor((extent.getMinX() - originX) / widthX))
			int endXGrid = (int)(Math.ceil((extent.getMaxX() - originX) / widthX))
			int startYGrid = (int)(Math.floor((extent.getMinY() - originY) / widthY))
			int endYGrid = (int)(Math.ceil((extent.getMaxY() - originY) / widthY))
			int cols = (int)(Math.abs(endXGrid - startXGrid))
			int rows = (int)(Math.abs(endYGrid - startYGrid))
			int numTiles = rows * cols

			if (numTiles > 32767) {
				pluginHost.showFeedback("There are too many output tiles.\nChoose a larger grid width.");
				return
			}
			
			// set up the output files of the shapefile and the dbf
            DBFField[] fields = new DBFField[1]

			fields[0] = new DBFField()
            fields[0].setName("FID")
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC)
            fields[0].setFieldLength(10)
            fields[0].setDecimalCount(0)
            
//            fields[0] = new DBFField()
//            fields[0].setName("Z")
//            fields[0].setDataType(DBFField.DBFDataType.NUMERIC)
//            fields[0].setFieldLength(10)
//            fields[0].setDecimalCount(5)
//
//            fields[1] = new DBFField()
//            fields[1].setName("INTENSITY")
//            fields[1].setDataType(DBFField.DBFDataType.NUMERIC)
//            fields[1].setFieldLength(10)
//            fields[1].setDecimalCount(5)

            //ShapeFile[] outputFiles = new ShapeFile[numTiles];
            
            double x, y, z, intensity
            int row, col
			int tileNum = 0
			int FID = 0
			int progress
			int oldProgress = -1

			// first figure out the tile number of each point
			short[] tileData = new short[totalPoints]
			
			for (int p = 0; p < totalPoints; p++) {
				PointRecord point = las.getPointRecord(p)
				x = point.getX()
				y = point.getY()
				col = (int)Math.floor((x - originX) / widthX) - startXGrid // relative to the grid edge
				row = (int)Math.floor((y - originY) / widthY) - startYGrid // relative to the grid edge
				tileData[p] = (short)(row * cols + col)
                progress = (int)(100f * p / (totalPoints - 1))
    			if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 1 of 3:", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}        		
			}

			// figure out the last point for each tile, so that the shapefile can be closed afterwards
			int[] firstPointNum = new int[numTiles]
			Arrays.fill(firstPointNum, (int)-1)
			int[] lastPointNum = new int[numTiles]
			int[] numPointsInTile = new int[numTiles]
			for (int p = 0; p < totalPoints; p++) {
				lastPointNum[tileData[p]] = p
				numPointsInTile[tileData[p]]++
				if (firstPointNum[tileData[p]] < 0) {
					firstPointNum[tileData[p]] = p
				}
				progress = (int)(100f * p / (totalPoints - 1))
    			if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 2 of 3:", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}        		
			}

			boolean[] outputTile = new boolean[numTiles]
			for (tileNum = 0; tileNum < numTiles; tileNum++) {
				if (numPointsInTile[tileNum] > minPoints) {
					outputTile[tileNum] = true
				}
			}

			int minRow = 999999
			int minCol = 999999
			for (tileNum = 0; tileNum < numTiles; tileNum++) {
				if (outputTile[tileNum]) {
					row = (int)Math.floor(tileNum / cols)
					col = tileNum % cols
					if (row < minRow) { minRow = row }
					if (col < minCol) { minCol = col }
				}
			}

			int numPointsOutput = 0
			int numTilesCreated = 0
			FID = 1
			for (tileNum = 0; tileNum < numTiles; tileNum++) {
				if (outputTile[tileNum]) {
					row = (int)Math.floor(tileNum / cols)
					col = tileNum % cols
					String outputFile = lasFile.replace(".las", "_Row${row - minRow + 1}_Col${col - minCol + 1}.shp")
					ShapeFile sf = new ShapeFile(outputFile, ShapeType.MULTIPOINTZ, fields)
					double[][] xyData = new double[numPointsInTile[tileNum]][2]
					double[] zData = new double[numPointsInTile[tileNum]]
					double[] mData = new double[numPointsInTile[tileNum]]
					int q = 0
					for (int p = firstPointNum[tileNum]; p <= lastPointNum[tileNum]; p++) {
						if (tileData[p] == tileNum) {
							PointRecord point = las.getPointRecord(p)
							//x = point.getX()
							//y = point.getY()
							xyData[q][0] = point.getX()
							xyData[q][1] = point.getY()
							zData[q] = point.getZ()
							mData[q] = point.getIntensity()
							q++
							//z = point.getZ()
							//whitebox.geospatialfiles.shapefile.Point wbPoint = new whitebox.geospatialfiles.shapefile.Point(x, y);                  
	                		//PointZ wbPoint = new PointZ(x, y, z, point.getIntensity());                  
	                		//Object[] rowData = new Object[1]
			                //rowData[0] = new Double()
			                //rowData[0] = new Double(FID)
			                //FID++
			                //sf.addRecord(wbPoint, rowData)
						}
					}

					MultiPointZ wbPoint = new MultiPointZ(xyData, zData, mData);
            		Object[] rowData = new Object[1]
	                rowData[0] = new Double(FID)
	                FID++
	                sf.addRecord(wbPoint, rowData)
			                
					sf.write()
					numTilesCreated++
				}
				numPointsOutput += numPointsInTile[tileNum]
				progress = (int)(100f * numPointsOutput / (totalPoints - 1))
    			if (progress != oldProgress) {
					pluginHost.updateProgress("Loop 3 of 3:", progress)
        			oldProgress = progress
        			// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
        		}
			}
			
			// create the index shapefile
			fields = new DBFField[4]

            fields[0] = new DBFField()
            fields[0].setName("TILE_NUM")
            fields[0].setDataType(DBFField.DBFDataType.NUMERIC)
            fields[0].setFieldLength(10)
            fields[0].setDecimalCount(0)

			fields[1] = new DBFField()
            fields[1].setName("ROW")
            fields[1].setDataType(DBFField.DBFDataType.NUMERIC)
            fields[1].setFieldLength(10)
            fields[1].setDecimalCount(0)

			fields[2] = new DBFField()
            fields[2].setName("COLUMN")
            fields[2].setDataType(DBFField.DBFDataType.NUMERIC)
            fields[2].setFieldLength(10)
            fields[2].setDecimalCount(0)

            fields[3] = new DBFField()
            fields[3].setName("NUM_POINTS")
            fields[3].setDataType(DBFField.DBFDataType.NUMERIC)
            fields[3].setFieldLength(10)
            fields[3].setDecimalCount(0)

			String outputFile = lasFile.replace(".las", "_index.shp")
			ShapeFile indexFile = new ShapeFile(outputFile, ShapeType.POLYGON, fields)
			int[] parts = [0];
			double north, south, east, west
            for (row = 0; row < rows; row++) {
				for (col = 0; col < cols; col++) {
					tileNum = row * cols + col
					if (numPointsInTile[tileNum] > minPoints) {
						west = (col + startXGrid) * widthX + originX
						east = (col + 1 + startXGrid) * widthX + originX
						south = (row + startYGrid) * widthY + originY
						north = (row + 1 + startYGrid) * widthY + originY
						
						PointsList points = new PointsList();
						points.addPoint(west, north);
			            points.addPoint(east, north);
			            points.addPoint(east, south);
			            points.addPoint(west, south);
			            points.addPoint(west, north);
			            
			            Polygon poly = new Polygon(parts, points.getPointsArray());
			            Object[] rowData = new Object[4];
			            rowData[0] = new Double(tileNum + 1);
			            rowData[1] = new Double(row - minRow + 1);
			            rowData[2] = new Double(col - minCol + 1);
			            rowData[3] = new Double(numPointsInTile[tileNum]);
			            indexFile.addRecord(poly, rowData);
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
			indexFile.write()
			pluginHost.returnData(outputFile)
			//pluginHost.showFeedback("Operation complete.\n${numTilesCreated} tiles were created.")
		
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
	def myClass = new TileLasFile(pluginHost, args, name, descriptiveName)
}
