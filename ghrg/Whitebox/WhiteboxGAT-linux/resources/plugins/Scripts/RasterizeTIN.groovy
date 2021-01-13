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
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import com.vividsolutions.jts.geom.*
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import whitebox.structures.BoundingBox
import whitebox.structures.RowPriorityGridCell
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D
import org.apache.commons.math3.geometry.euclidean.threed.Plane
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "RasterizeTIN"
def descriptiveName = "Rasterize TIN"
def description = "Rasterizes a triangular irregular network (TIN)"
def toolboxes = ["Interpolation"]

public class RasterizeTIN implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public RasterizeTIN(WhiteboxPluginHost pluginHost, 
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
            sd.setHelpFile("RasterizeTIN")
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "RasterizeTIN.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
        	DialogFile dfIn = sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            sd.addDialogFile("Output file", "Output Raster File:", "saveAs", "Whitebox Raster Files (*.dep), DEP", true, false)
            sd.addDialogDataInput("Output raster cell size.", "Cell Size (optional):", "", true, false)
            sd.addDialogFile("Input base file", "Base Raster File (optional):", "open", "Whitebox Raster Files (*.dep), DEP", true, true)
            
            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
        }
    }

    // The CompileStatic annotation can be used to significantly
    // improve the performance of a Groovy script to nearly 
    // that of native Java code.
    //@CompileStatic
    private void execute(String[] args) {
        try {
        	int progress, oldProgress
            int row, col, i
	        double rowYCoord, value, x, y, z
	        double cellSize = -1.0
	        int cols, rows, topRow, bottomRow, numEdges, stCol, endCol
	        double nodata = -32768.0
	        double east, west, north, south
	        ArrayList<com.vividsolutions.jts.geom.Geometry> pointList = new ArrayList<>()
			double[][] points
            double[] zArray
			Object[] recData
			ArrayList<Integer> edgeList = new ArrayList<>()
			double x1, y1, x2, y2, xPrime
        	boolean foundIntersection
        	DecimalFormat df = new DecimalFormat("###,###,###,###")
	        
            if (args.length != 4) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            
            // read the input parameters
            String inputFile = args[0]
            String outputFile = args[1]
			if (!args[2].toLowerCase().contains("not specified")) {
	            cellSize = Double.parseDouble(args[2]);
	        }
	        String baseFileHeader = args[3]
	        if (baseFileHeader == null || baseFileHeader.isEmpty()) {
	        	baseFileHeader = "not specified"
	        }

			// initialize the input TIN
	        ShapeFile input = new ShapeFile(inputFile)
			ShapeType shapeType = input.getShapeType()
            if (shapeType != ShapeType.POLYGONZ) {
            	pluginHost.showFeedback("The input shapefile should be of a POLYGONZ shapetype.")
            	return
            }

			// initialize the output raster
            WhiteboxRaster output;
            if ((cellSize > 0) || ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified")))) {
                if ((cellSize < 0) & (baseFileHeader.toLowerCase().contains("not specified"))) {
                    cellSize = Math.min((input.getyMax() - input.getyMin()) / 500.0,
                            (input.getxMax() - input.getxMin()) / 500.0);
                }
                north = input.getyMax() + cellSize / 2.0;
                south = input.getyMin() - cellSize / 2.0;
                east = input.getxMax() + cellSize / 2.0;
                west = input.getxMin() - cellSize / 2.0;
                rows = (int) (Math.ceil((north - south) / cellSize));
                cols = (int) (Math.ceil((east - west) / cellSize));

                // update west and south
                east = west + cols * cellSize;
                south = north - rows * cellSize;

                output = new WhiteboxRaster(outputFile, north, south, east, west,
                        rows, cols, WhiteboxRasterBase.DataScale.CONTINUOUS,
                        WhiteboxRasterBase.DataType.FLOAT, nodata, nodata);
            } else {
                output = new WhiteboxRaster(outputFile, "rw",
                        baseFileHeader, WhiteboxRasterBase.DataType.FLOAT, nodata);
                rows = output.getNumberRows()
                cols = output.getNumberColumns()
            }

            output.setPreferredPalette("high_relief.pal")

            // first sort the records based on their maxY coordinate. This will
            // help reduce the amount of disc IO for larger rasters.
            ArrayList<RecordInfo> myList = new ArrayList<>();
            
            for (ShapeFileRecord record : input.records) {
                i = record.getRecordNumber();
                y = ((PolygonZ)(record.getGeometry())).getBox().getMaxY()
                myList.add(new RecordInfo(y, i));
            }
            
            Collections.sort(myList);

            long heapSize = Runtime.getRuntime().totalMemory()
            int flushSize = (int)(heapSize / 32)
            int j, numCellsToWrite
            Queue<RowPriorityGridCell> pq = new PriorityQueue<>(flushSize)
            
            RowPriorityGridCell cell
            int numFeatures = input.getNumberOfRecords()
            int count = 0
            ShapeFileRecord record
            oldProgress = -1
			for (RecordInfo ri : myList) {
                record = input.getRecord(ri.recNumber - 1);
				PolygonZ pz = (PolygonZ)(record.getGeometry())
				points = record.getGeometry().getPoints()
				if (points.length != 4) {
					pluginHost.showFeedback("The input shapefile does not appear to be a TIN created using the Create TIN tool.")
					return
				}
				
				zArray = pz.getzArray()

				Vector3D pt1 = new Vector3D(points[0][0], points[0][1], zArray[0])
				Vector3D pt2 = new Vector3D(points[1][0], points[1][1], zArray[1])
				Vector3D pt3 = new Vector3D(points[2][0], points[2][1], zArray[2])
				Plane plane = new Plane(pt1, pt2, pt3)
				
				Vector3D normal = plane.getNormal()
		
				def A = normal.getX()
				def B = normal.getY()
				def C = normal.getZ()
				def D = -(A * pt1.getX() + B * pt1.getY() + C * pt1.getZ())


				topRow = output.getRowFromYCoordinate(pz.getYMax());
                bottomRow = output.getRowFromYCoordinate(pz.getYMin());

				for (row = topRow; row <= bottomRow; row++) {
                    edgeList.clear();
                    foundIntersection = false;
                    rowYCoord = output.getYCoordinateFromRow(row);
                    // find the x-coordinates of each of the edges that 
                    // intersect this row's y coordinate

                    for (i = 0; i < 3; i++) {
                        if (isBetween(rowYCoord, points[i][1], points[i + 1][1])) {
                            y1 = points[i][1];
                            y2 = points[i + 1][1];
                            if (y2 != y1) {
                                x1 = points[i][0];
                                x2 = points[i + 1][0];

                                // calculate the intersection point
                                xPrime = (x1 + (rowYCoord - y1) / (y2 - y1) * (x2 - x1));
                                edgeList.add(new Integer(output.getColumnFromXCoordinate(xPrime)));
                                foundIntersection = true;
                            }
                        }
                    }

                    if (foundIntersection) {
                        numEdges = edgeList.size();
                        if (numEdges == 2) {
                            stCol = Math.min(edgeList.get(0), edgeList.get(1));
                            endCol = Math.max(edgeList.get(0), edgeList.get(1));
                            for (col = stCol; col <= endCol; col++) {
                            	x = output.getXCoordinateFromColumn(col)
                            	y = output.getYCoordinateFromRow(row)
                            	z = -(A * x + B * y + D) / C
                            	cell = new RowPriorityGridCell(row, col, z)
                                pq.add(cell)
                            }
                        } else {
                            //sort the edges.
                            int[] edgeArray = new int[numEdges];
                            edgeList.toArray(edgeArray);
                            Arrays.sort(edgeArray);

                            boolean fillFlag = true;
                            for (i = 0; i < numEdges - 1; i++) {
                                stCol = edgeArray[i];
                                endCol = edgeArray[i + 1];
                                if (fillFlag) {
                                    for (col = stCol; col <= endCol; col++) {
                                        x = output.getXCoordinateFromColumn(col)
                            			y = output.getYCoordinateFromRow(row)
                            			z = -(A * x + B * y + D) / C
                                        pq.add(new RowPriorityGridCell(row, col, z))
                                    }
                                }
                                fillFlag = !fillFlag;
                            }
                        }
                    }
                }

				if (pq.size() >= flushSize) {
                    j = 0;
                    numCellsToWrite = pq.size();
                    while (pq.size() > 0) {
                        cell = pq.poll();
                        output.setValue(cell.row, cell.col, cell.z);
                        j++;
                        if (j % 1000 == 0) {
                            pluginHost.updateProgress("Writing to Output (" + df.format(j) + " of " + df.format(numCellsToWrite) + "):", progress)
                        }
                    }
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

            j = 0;
            numCellsToWrite = pq.size();
            while (pq.size() > 0) {
                cell = pq.poll()
                output.setValue(cell.row, cell.col, cell.z)
                j++
                if (j % 1000 == 0) {
                    pluginHost.updateProgress("Writing to Output (" + df.format(j) + " of " + df.format(numCellsToWrite) + "):", progress)
                }
            }

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

    // Return true if val is between theshold1 and theshold2.
    private static boolean isBetween(double val, double threshold1, double threshold2) {
        if (val == threshold1 || val == threshold2) {
            return true;
        }
        return threshold2 > threshold1 ? val > threshold1 && val < threshold2 : val > threshold2 && val < threshold1;
    }

    private class RecordInfo implements Comparable<RecordInfo> {

        public double maxY;
        public int recNumber;
        
        public RecordInfo(double maxY, int recNumber) {
            this.maxY = maxY;
            this.recNumber = recNumber;
        }

        @Override
        public int compareTo(RecordInfo other) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.maxY < other.maxY) {
                return BEFORE;
            } else if (this.maxY > other.maxY) {
                return AFTER;
            }

            if (this.recNumber < other.recNumber) {
                return BEFORE;
            } else if (this.recNumber > other.recNumber) {
                return AFTER;
            }

            return EQUAL;
        }
    }
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new RasterizeTIN(pluginHost, args, descriptiveName)
}
