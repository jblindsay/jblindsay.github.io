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
import java.io.File
import java.util.Date
import java.util.ArrayList
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities;
import com.vividsolutions.jts.geom.*
import com.vividsolutions.jts.triangulate.DelaunayTriangulationBuilder
import whitebox.geospatialfiles.VectorLayerInfo
import whitebox.geospatialfiles.shapefile.attributes.*
import whitebox.geospatialfiles.shapefile.ShapeFileRecord
import whitebox.utilities.Topology
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D
import org.apache.commons.math3.geometry.euclidean.threed.Plane
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ConstructTIN"
def descriptiveName = "Construct TIN"
def description = "Constructs a triangular irregular network (TIN) model from vector points"
def toolboxes = ["Interpolation"]

public class ConstructTIN implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName
	
    public ConstructTIN(WhiteboxPluginHost pluginHost, 
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
            sd.setHelpFile("ConstructTIN")
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "ConstructTIN.groovy"
            sd.setSourceFile(scriptFile)
			
            // add some components to the dialog
        	//DialogFile dfs = sd.addDialogFile("Input file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", true, false)
            DialogFieldSelector dfs = sd.addDialogFieldSelector("Input file and height field.", "Input Height Field:", false)
            DialogCheckBox dcb = sd.addDialogCheckBox("Use z-values", "Use z-values", false)
            dcb.setVisible(false)
            sd.addDialogFile("Output file", "Output Vector File:", "saveAs", "Vector Files (*.shp), SHP", true, false)
            
            def listener = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = dfs.getValue()
            		if (value != null && !value.isEmpty()) {
            			value = value.trim()
            			String[] strArray = dfs.getValue().split(";")
            			String fileName = strArray[0]
            			File file = new File(fileName)
            			if (file.exists()) {
	            			ShapeFile shapefile = new ShapeFile(fileName)
		            		if (shapefile.getShapeType().getDimension() == ShapeTypeDimension.Z) {
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
            dfs.addPropertyChangeListener(listener)
            
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
            double fx, fy, tanSlope, aspect, hillshade, z, x, y, d
			double term1, term2, term3
			double azimuth = Math.toRadians(315.0 - 90)
			double altitude = Math.toRadians(30.0)
			double sinTheta = Math.sin(altitude)
			double cosTheta = Math.cos(altitude)
            if (args.length != 3) {
                pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
            }
            // read the input parameters
            String[] inputData = args[0].split(";")
            boolean useZValues = Boolean.parseBoolean(args[1])
            String outputFile = args[2]
			String inputFile = inputData[0]
            
            ShapeFile input = new ShapeFile(inputFile)
			AttributeTable table = new AttributeTable(inputFile.replace(".shp", ".dbf"))
			ShapeType shapeType = input.getShapeType()
            if (shapeType.getDimension() != ShapeTypeDimension.Z && useZValues) {
            	useZValues = false
            }
			int heightField = -1
            if (inputData.length == 2 && !inputData[1].trim().isEmpty()) {
            	String heightFieldName = inputData[1].trim()
            	String[] fieldNames = input.getAttributeTableFields()
            	for (int i = 0; i < fieldNames.length; i++) {
            		if (fieldNames[i].trim().equals(heightFieldName)) {
            			heightField = i
            			break
            		}
            	}
            } else if (!useZValues) {
            	pluginHost.showFeedback("A field within the input file's attribute table must be selected to assign point heights.")
            	return
            }
			
			ArrayList<com.vividsolutions.jts.geom.Geometry> pointList = new ArrayList<>()
			
			double[][] point
			Object[] recData
			Coordinate c
			GeometryFactory geomFactory = new GeometryFactory()
			int i = 0
			int numFeatures = input.getNumberOfRecords()
			oldProgress = -1
			if (!useZValues) {
				for (ShapeFileRecord record : input.records) {
					recData = table.getRecord(i)
					z = (Double)(recData[heightField])
					point = record.getGeometry().getPoints()
					for (int p = 0; p < point.length; p++) {
						x = point[p][0]
						y = point[p][1]
						c = new Coordinate(x, y, z)
						pointList.add(geomFactory.createPoint(c))
					}
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            	}
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			} else {
				for (ShapeFileRecord record : input.records) {
					if (shapeType.getBaseType() == ShapeType.POINT) {
						PointZ ptz = (PointZ)(record.getGeometry())
                		z = ptz.getZ()
                		x = ptz.getX()
						y = ptz.getY()
						c = new Coordinate(x, y, z)
						pointList.add(geomFactory.createPoint(c))
					} else if (shapeType.getBaseType() == ShapeType.MULTIPOINT) {
						MultiPointZ plz = (MultiPointZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = plz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							c = new Coordinate(x, y, z)
							pointList.add(geomFactory.createPoint(c))
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYLINE) {
						PolyLineZ plz = (PolyLineZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = plz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							c = new Coordinate(x, y, z)
							pointList.add(geomFactory.createPoint(c))
						}
					} else if (shapeType.getBaseType() == ShapeType.POLYGON) {
						PolygonZ pz = (PolygonZ)(record.getGeometry())
						point = record.getGeometry().getPoints()
						double[] zArray = pz.getzArray()
						for (int p = 0; p < point.length; p++) {
							x = point[p][0]
							y = point[p][1]
							z = zArray[p]
							c = new Coordinate(x, y, z)
							pointList.add(geomFactory.createPoint(c))
						}
					}
					
					i++
	                progress = (int)(100f * i / numFeatures)
	            	if (progress != oldProgress) {
						pluginHost.updateProgress("Reading Points:", progress)
	            		oldProgress = progress
	            	}
	            	// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}
			
			com.vividsolutions.jts.geom.Geometry geom = geomFactory.buildGeometry(pointList)
			
			DelaunayTriangulationBuilder dtb = new DelaunayTriangulationBuilder()
			dtb.setSites(geom)
			com.vividsolutions.jts.geom.Geometry polys = dtb.getTriangles(geomFactory)

			// set up the output files of the shapefile and the dbf
			DBFField[] fields = new DBFField[9];
			fields[0] = new DBFField();
			fields[0].setName("FID");
			fields[0].setDataType(DBFField.DBFDataType.NUMERIC);
			fields[0].setFieldLength(10);
			fields[0].setDecimalCount(0);
			
			fields[1] = new DBFField();
			fields[1].setName("EQ_AX");
			fields[1].setDataType(DBFField.DBFDataType.FLOAT);
			fields[1].setFieldLength(10);
			fields[1].setDecimalCount(5);
			
			fields[2] = new DBFField();
			fields[2].setName("EQ_BY");
			fields[2].setDataType(DBFField.DBFDataType.FLOAT);
			fields[2].setFieldLength(10);
			fields[2].setDecimalCount(5);
			
			fields[3] = new DBFField();
			fields[3].setName("EQ_CZ");
			fields[3].setDataType(DBFField.DBFDataType.FLOAT);
			fields[3].setFieldLength(10);
			fields[3].setDecimalCount(5);
			
			fields[4] = new DBFField();
			fields[4].setName("EQ_D");
			fields[4].setDataType(DBFField.DBFDataType.FLOAT);
			fields[4].setFieldLength(10);
			fields[4].setDecimalCount(5);
			
			fields[5] = new DBFField();
			fields[5].setName("SLOPE");
			fields[5].setDataType(DBFField.DBFDataType.FLOAT);
			fields[5].setFieldLength(10);
			fields[5].setDecimalCount(4);
			
			fields[6] = new DBFField();
			fields[6].setName("ASPECT");
			fields[6].setDataType(DBFField.DBFDataType.FLOAT);
			fields[6].setFieldLength(10);
			fields[6].setDecimalCount(3);
			
			fields[7] = new DBFField();
			fields[7].setName("HILLSHADE");
			fields[7].setDataType(DBFField.DBFDataType.FLOAT);
			fields[7].setFieldLength(10);
			fields[7].setDecimalCount(0);
			
			fields[8] = new DBFField();
			fields[8].setName("CNTR_PT_Z");
			fields[8].setDataType(DBFField.DBFDataType.FLOAT);
			fields[8].setFieldLength(10);
			fields[8].setDecimalCount(3);
			
			ShapeFile output = new ShapeFile(outputFile, ShapeType.POLYGONZ, fields);
			
			long[] histo = new long[256];
			int FID = 1
			numFeatures = polys.getNumGeometries()
			oldProgress = -1
			for (int a = 0; a < polys.getNumGeometries(); a++) {
			    com.vividsolutions.jts.geom.Geometry g = polys.getGeometryN(a);
			    if (g instanceof com.vividsolutions.jts.geom.Polygon) {
			        com.vividsolutions.jts.geom.Polygon p = (com.vividsolutions.jts.geom.Polygon) g;
			        ArrayList<ShapefilePoint> pnts = new ArrayList<>();
			        int[] parts = new int[p.getNumInteriorRing() + 1];
			        
			        Coordinate[] buffCoords = p.getExteriorRing().getCoordinates();
			        if (!Topology.isLineClosed(buffCoords)) {
			            pluginHost.showFeedback("Exterior ring not closed.");
			        }
			        if (Topology.isClockwisePolygon(buffCoords)) {
			            for (i = 0; i < buffCoords.length; i++) {
			                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y, buffCoords[i].z, 0.0));
			            }
			        } else {
			            for (i = buffCoords.length - 1; i >= 0; i--) {
			                pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y, buffCoords[i].z, 0.0));
			            }
			        }
			
			        for (int b = 0; b < p.getNumInteriorRing(); b++) {
			            parts[b + 1] = pnts.size();
			            buffCoords = p.getInteriorRingN(b).getCoordinates();
			            if (!Topology.isLineClosed(buffCoords)) {
			                pluginHost.showFeedback("Interior ring not closed.");
			            }
			            if (Topology.isClockwisePolygon(buffCoords)) {
			                for (i = buffCoords.length - 1; i >= 0; i--) {
			                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y, buffCoords[i].z, 0.0));
			                }
			            } else {
			                for (i = 0; i < buffCoords.length; i++) {
			                    pnts.add(new ShapefilePoint(buffCoords[i].x, buffCoords[i].y, buffCoords[i].z, 0.0));
			                }
			            }
			        }
			
			        PointsList pl = new PointsList(pnts);
					PolygonZ wbPoly = new PolygonZ(parts, pl.getPointsArray(), pl.getZArray());
			
					double centroidX = 0
					double centroidY = 0
					ShapefilePoint sfp = pl.getPoint(0)
					Vector3D pt1 = new Vector3D(sfp.x, sfp.y, sfp.z)
					centroidX += sfp.x
					centroidY += sfp.y
					sfp = pl.getPoint(1)
					Vector3D pt2 = new Vector3D(sfp.x, sfp.y, sfp.z)
					centroidX += sfp.x
					centroidY += sfp.y
					sfp = pl.getPoint(2)
					Vector3D pt3 = new Vector3D(sfp.x, sfp.y, sfp.z)
					centroidX += sfp.x
					centroidY += sfp.y
					Plane plane = new Plane(pt1, pt2, pt3)
					
					centroidX = centroidX / 3.0
					centroidY = centroidY / 3.0
			
					Vector3D normal = plane.getNormal()
			
					def A = normal.getX()
					def B = normal.getY()
					def C = normal.getZ()
					def D = -(A * pt1.getX() + B * pt1.getY() + C * pt1.getZ())
					
					if (C != 0) {
						fx = -A / C
						fy = -B / C
						if (fx != 0) {
				            tanSlope = Math.sqrt(fx * fx + fy * fy);
				            aspect = Math.toRadians(180 - Math.toDegrees(Math.atan(fy / fx)) + 90 * (fx / Math.abs(fx)))
				            term1 = tanSlope / Math.sqrt(1 + tanSlope * tanSlope);
				            term2 = sinTheta / tanSlope;
				            term3 = cosTheta * Math.sin(azimuth - aspect);
				            hillshade = term1 * (term2 - term3);
				        } else {
				            hillshade = 0.5;
				        }
				        hillshade = (int)(hillshade * 255);
				        if (hillshade < 0) {
				            hillshade = 0;
				        }
					} else {
						hillshade = 0.0
					}
					histo[(int)hillshade]++;
			
					z = -(A * centroidX + B * centroidY + D) / C
					
			        Object[] rowData = new Object[9]
			        rowData[0] = new Double(FID)
			        rowData[1] = new Double(A)
			        rowData[2] = new Double(B)
			        rowData[3] = new Double(C)
			        rowData[4] = new Double(D)
			        rowData[5] = new Double(Math.toDegrees(Math.atan(tanSlope)))
			        rowData[6] = new Double(Math.toDegrees(aspect))
			        rowData[7] = new Double(hillshade)
			        rowData[8] = new Double(z)
			        
			        output.addRecord(wbPoly, rowData);
			        FID++
			    }
			    progress = (int)(100f * a / (numFeatures - 1))
            	if (progress != oldProgress) {
					pluginHost.updateProgress("Saving Data:", progress)
            		oldProgress = progress
            	}
            	// check to see if the user has requested a cancellation
				if (pluginHost.isRequestForOperationCancelSet()) {
					pluginHost.showFeedback("Operation cancelled")
					return
				}
			}
			
			// trim the display min and max values by 2%
			int newMin = 0;
			int newMax = 0;
			double targetCellNum = FID * 0.02;
			long sum = 0;
			for (int k = 0; k < 256; k++) {
			    sum += histo[k];
			    if (sum >= targetCellNum) {
			        newMin = k;
			        break;
			    }
			}
			
			sum = 0;
			for (int k = 255; k >= 0; k--) {
			    sum += histo[k];
			    if (sum >= targetCellNum) {
			        newMax = k;
			        break;
			    }
			}
			
			output.write()
			table = output.getAttributeTable()
			double range = newMax - newMin
			for (int r = 0; r < output.getNumberOfRecords(); r++) {
				recData = table.getRecord(r)
                hillshade = (Double)(recData[7])
                hillshade = 255.0 * (hillshade - newMin) / range
                if (hillshade < 0) { hillshade = 0 }
                if (hillshade > 255) { hillshade = 255 }
				recData[7] = new Double(hillshade)
				table.updateRecord(r, recData)
			}

			output.write()
			
			def paletteDirectory = pluginHost.getResourcesDirectory() + "palettes" + File.separator
			VectorLayerInfo vli = new VectorLayerInfo(outputFile, paletteDirectory, 255i, -1)
			vli.setFilledWithOneColour(false)
			vli.setFillAttribute("HILLSHADE")
			vli.setPaletteFile(paletteDirectory + "grey.pal")
			vli.setPaletteScaled(true)
//			vli.setMinimumValue(newMin)
//			vli.setMaximumValue(newMax)
			vli.setRecordsColourData()
			if (numFeatures > 5000) {
				vli.setOutlined(false)
			}
			
			pluginHost.returnData(vli)
			
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
    def f = new ConstructTIN(pluginHost, args, descriptiveName)
}
