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
import java.beans.PropertyChangeEvent
import java.beans.PropertyChangeListener
import java.io.File
import java.net.URL
import java.util.zip.*
import java.util.Date
import java.util.ArrayList
import java.util.Arrays
import javax.swing.*
import org.apache.commons.io.FilenameUtils
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.ui.plugin_dialog.*
import whitebox.utilities.FileUtilities
import whiteboxgis.WhiteboxGui
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "RetrieveSRTMData"
def descriptiveName = "Retrieve SRTM DEM Data"
def description = "Downloads Shuttle Radar Topography Mission DEM data."
def toolboxes = ["IOTools"]

// The following lines are necessary for the script 
// to be recognized as a menu extension.
parentMenu = "Data Layers"
menuLabel = "Retrieve SRTM DEM Data"

public class RetrieveSRTMData implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd;
    private String descriptiveName

	public RetrieveSRTMData(WhiteboxPluginHost pluginHost, 
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
            sd.addDialogComboBox("SRTM Dataset", "SRTM Dataset:", ["3-arcsecond (Global)", "1-arcsecond (US only)"], 0)
            
            String[] lats = new String[119]
            int lat = 61
            for (int a in 0..118) {
				lats[a] = "${lat}"
				lat--
			}
            sd.addDialogComboBox("Starting latitude", "Starting Latitude:", lats, 18)
            sd.addDialogComboBox("Ending latitude", "Ending Latitude:", lats, 17)
			
            String[] longs = new String[361]
            for (int a in 0..360) {
				int k = -1 * (a - 180)
				longs[a] = "${k}"
			}
            sd.addDialogComboBox("Starting longitude", "Starting Longitude:", longs, 261)
            sd.addDialogComboBox("Ending longitude", "Ending Longitude:", longs, 260)
            
            sd.addDialogCheckBox("Fill missing data holes?", "Fill missing data holes?", true)
			sd.addDialogCheckBox("Replace ocean values with nodata?", "Replace ocean values with nodata", true)
			DialogCheckBox cb = sd.addDialogCheckBox("Mosaic DEM tiles?", "Mosaic DEM tiles?", true)
			DialogFile df = sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)

            //Listener for cb            
            def lstr = { evt -> if (evt.getPropertyName().equals("value")) { 
            		String value = cb.getValue()
            		if (!value.isEmpty() && value != null) { 
            			if (cb.getValue().equals("true")) {
            				df.setVisible(true)
		            	} else {
		            		df.setVisible(false)
		            	}
            		}
            	} 
            } as PropertyChangeListener
            cb.addPropertyChangeListener(lstr)

            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
        }
    }

    @CompileStatic
    private void execute(String[] args) {
        try {
            WhiteboxGui wg = (WhiteboxGui)(pluginHost)
			int north, south, east, west, progress, oldProgress
			String wd = pluginHost.getWorkingDirectory()
			
			if (args.length < 5) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
                return
			}

			boolean is3Arcsecond = true
			if (args[0].contains("1") || args[0].toLowerCase().contains("one")) {
				is3Arcsecond = false
			}
			int startingLat = Integer.parseInt(args[1])
			int endingLat = Integer.parseInt(args[2])
			int startingLon = Integer.parseInt(args[3])
			int endingLon = Integer.parseInt(args[4])
			
			boolean fillDataHoles = true
			boolean mosaicImages = true
			boolean fixElevationsBool = true

			if (args.length >= 6) {
				fillDataHoles = Boolean.parseBoolean(args[5])
			}
			if (args.length >= 7) {
				fixElevationsBool = Boolean.parseBoolean(args[6])
			}
			if (args.length >= 8) {
				mosaicImages = Boolean.parseBoolean(args[7])
			}
			String outputFile = wd + "SRTM mosaic.dep"
			if (args.length >= 9) {
				outputFile = args[8]
			}

			if (startingLat < endingLat) {
				north = endingLat - 1
				south = startingLat
			} else {
				south = endingLat
				north = startingLat - 1
			}
			
			if (startingLon < endingLon) {
				east = endingLon - 1
				west = startingLon
			} else {
				west = endingLon
				east = startingLon - 1
			}

			int numTiles = (north - south + 1) * (east - west + 1)

			if (numTiles > 500) {
				if (pluginHost.showFeedback("You will be importing more than 500 tiles." +
				"\nAre you sure this is what you want to do?", 
				JOptionPane.YES_NO_OPTION, 
				JOptionPane.QUESTION_MESSAGE) == JOptionPane.NO_OPTION) {
                	return
				}
			}
			
			int regionNum = -1
			String region = ""
			boolean foundRegion = false
			oldProgress = -1
			int numDownloadedTiles = 0

			String startDir = "http://dds.cr.usgs.gov/srtm/version2_1/SRTM3"
			String[] directories = ["North_America", "South_America", "Eurasia", "Africa", "Australia", "Islands"]
			if (!is3Arcsecond) {
				startDir = "http://dds.cr.usgs.gov/srtm/version2_1/SRTM1"
				directories = ["Region_01", "Region_02", "Region_03", "Region_04", "Region_05", "Region_06", "Region_07"]
			}
					
			def importedFiles = new ArrayList<String>()
			for (int lat in south..north) {
				for (int lon in west..east) {
					String northStr, westStr
					if (lat >= 0) {
						if (lat > 9) {
							northStr = "N${lat}"
						} else {
							northStr = "N0${lat}"
						}
					} else {
						if (lat < -9) {
							northStr = "S${-lat}"
						} else {
							northStr = "S0${-lat}"
						}
					}
			
					if (lon >= 0) {
						if (lon > 99) {
							westStr = "E${lon}"
						} else if (lon > 9) {
							westStr = "E0${lon}"
						} else {
							westStr = "E00${lon}"
						}
					} else {
						if (lon < -99) {
							westStr = "W${-lon}"
						} else if (lon < -9) {
							westStr = "W0${-lon}"
						} else {
							westStr = "W00${-lon}"
						}
					}
					
					def shortName = "${northStr}${westStr}"
					
					if (!foundRegion) {
						// which region is the data in?
						for (int i = 0; i < directories.size(); i++) {
							region = directories[i]
							String downloadArtifact = "${startDir}/${region}/${shortName}.hgt.zip"
							String zipFileName = wd + "${shortName}.zip" 
							String outDir = (new File(zipFileName)).getParentFile()
							if (!outDir.endsWith(File.separator)) { outDir = outDir + File.separator }
							
							Downloader d = new Downloader(wg, zipFileName, downloadArtifact)
							int ret = d.download(shortName)
							if (ret == 0) {
								regionNum = i
								region = directories[regionNum]
								foundRegion = true
								//pluginHost.showFeedback(directories[i])
								String extractedFile = extractFolder(zipFileName)
								// import the file
								args = [extractedFile, outDir] 
								pluginHost.runPlugin("ImportSRTM", args, false, true)
								deleteFolder(new File("${outDir}${shortName}${File.separator}"))
								(new File(zipFileName)).delete()
								importedFiles.add("${outDir}${shortName}.dep")
								break
							}
						}
	
//						if (regionNum < 0) {
//							pluginHost.showFeedback("The data could not be located on the server.")
//	                		return
//						}
						
//						region = directories[regionNum]
//						foundRegion = true
					} else {
						String downloadArtifact = "${startDir}/${region}/${shortName}.hgt.zip"
						String zipFileName = wd + "${shortName}.zip" 
						String outDir = (new File(zipFileName)).getParentFile()
						if (!outDir.endsWith(File.separator)) { outDir = outDir + File.separator }
						
						Downloader d = new Downloader(wg, zipFileName, downloadArtifact)
						int ret = d.download(shortName)
						if (ret == 0) {
							String extractedFile = extractFolder(zipFileName)
							// import the file
							args = [extractedFile, outDir] 
							pluginHost.runPlugin("ImportSRTM", args, false, true)
							deleteFolder(new File("${outDir}${shortName}${File.separator}"))
							(new File(zipFileName)).delete()
							importedFiles.add("${outDir}${shortName}.dep")
							
						} else if (ret == 1) {
							// couldn't find the file
						} else if (ret == 2) {
							// try it a second time
							if (d.download(shortName) == 0) {
								String extractedFile = extractFolder(zipFileName)
								// import the file
								args = [extractedFile, outDir] 
								pluginHost.runPlugin("ImportSRTM", args, false, true)
								deleteFolder(new File("${outDir}${shortName}${File.separator}"))
								(new File(zipFileName)).delete()
								importedFiles.add("${outDir}${shortName}.dep")
								
							} else {
								pluginHost.showFeedback("There was a timeout error when trying to download ${shortName}")
							}
							//return
						} else if (ret == 3) {
							// try it a second time
//							if (d.download(shortName) == 0) {
//								String extractedFile = extractFolder(zipFileName)
//								// import the file
//								args = [extractedFile, outDir] 
//								pluginHost.runPlugin("ImportSRTM", args, false, true)
//								deleteFolder(new File("${outDir}${shortName}${File.separator}"))
//								(new File(zipFileName)).delete()
//								importedFiles.add("${outDir}${shortName}.dep")
//								
//							} else {
//								//pluginHost.showFeedback("There was a server access error when trying to download ${shortName}.\nThe SRTM FTP server may be down. Try again later.")
//							}
							//return
						} else {
							pluginHost.showFeedback("There was a problem downloading the data for ${shortName}")
							return
						}
					}

					numDownloadedTiles++
					progress = (int)(100f * numDownloadedTiles / numTiles)
					if (progress > oldProgress) {
						pluginHost.updateProgress("Downloading tiles from FTP site:", progress)
						oldProgress = progress
						
						// check to see if the user has requested a cancellation
						if (pluginHost.isRequestForOperationCancelSet()) {
							pluginHost.showFeedback("Operation cancelled")
							return
						}
					}
				}
			}

			if (importedFiles.size() == 0) {
				if (is3Arcsecond) {
					pluginHost.showFeedback("No SRTM tiles were downloaded. Either the coordinates of the area of interest \nis incorrect or the server is down.")
					return
				} else {
					pluginHost.showFeedback("No SRTM tiles were downloaded. Perhaps try the 3-arcsecond SRTM dataset with the area of interest.")
					return
				}
			}

			if (!fillDataHoles && !mosaicImages) {
				for (String ifile : importedFiles) {
					if (fixElevationsBool) {
					  fixElevations(ifile);
					}
					pluginHost.returnData(ifile) 
					pluginHost.zoomToFullExtent()
				}
        	}

			def filledFiles = new ArrayList<String>()
			if (fillDataHoles) {
				for (String f : importedFiles) { 
					String fillFile = f.replace(".dep", "_no_holes.dep")
					filledFiles.add(fillFile)
					args = [f, fillFile, "false"]
					pluginHost.runPlugin("FillMissingDataHoles", args, false, true)
					(new File(f)).delete()
					(new File(f.replace(".dep", ".tas"))).delete()
					if (!mosaicImages || numTiles == 1) {
						if (fixElevationsBool) {
							fixElevations(fillFile)
						}
						pluginHost.returnData(fillFile)
						pluginHost.zoomToFullExtent()
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
			}

			
			if (mosaicImages && numTiles > 1) {
				if (filledFiles.size() == 0) {
					for (String str : importedFiles) {
						filledFiles.add(str)
					}
				}
			
				String inputFiles = filledFiles.get(0)
				for (int k = 1; k < filledFiles.size(); k++) {
					inputFiles = inputFiles + ";" + filledFiles.get(k) 
				}
				def resamplingMethod = "nearest neighbour" 
				args = [inputFiles, outputFile, resamplingMethod] 
				pluginHost.runPlugin("Mosaic", args, false, true)

				// delete the filledFiles
				for (int k = 0; k < filledFiles.size(); k++) {
					(new File(filledFiles.get(k))).delete()
					(new File(filledFiles.get(k).replace(".dep", ".tas"))).delete()
				}
				if (fixElevationsBool) {
					fixElevations(outputFile)
				}
				pluginHost.returnData(outputFile)
				//pluginHost.zoomToFullExtent()
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

	@CompileStatic
    void fixElevations(String fileName) {
    	def wbr = new WhiteboxRaster(fileName, "rw")
    	if (wbr.getMinimumValue() < 0.01) { 
    		// fix ocean values so they are nodata
			int cols = wbr.getNumberColumns()
			int rows = wbr.getNumberRows()
			double nodata = wbr.getNoDataValue()
			double[] data
			double z
			int row, col
			int oldProgress = -1
			int progress
	        for (row = 0; row < rows; row++) {
	        	for (col = 0; col < cols; col++) {
	        		z = wbr.getValue(row, col)
	            	if (z <= 0.01 && z > -0.01 && z != nodata) {
	            		/* due to limits of floating point 
	            		 *  representation you can't simply 
	            		 *  check for zero values. */
	      
	            		wbr.setValue(row, col, nodata)
	            	}
	            }
	            
	            progress = (int)(100f * row / (rows - 1))
				if (progress > oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
	        }
			wbr.findMinAndMaxVals()
    	}
    	wbr.setPreferredPalette("high_relief.pal")
    	wbr.setXYUnits("degrees")
    	wbr.setZUnits("metres")
    	wbr.setProjection("geographic coordinates")
		wbr.close()
    }

	// This method is from NeilMonday, http://stackoverflow.com/questions/981578/how-to-unzip-files-recursively-in-java
	// Modified by John Lindsay
	static String extractFolder(String zipFile) throws ZipException, IOException {
	    String ret = ""
	    int BUFFER = 2048;
	    File file = new File(zipFile);
	
	    ZipFile zip = new ZipFile(file);
	    String newPath = zipFile.substring(0, zipFile.length() - 4);
	
	    new File(newPath).mkdir();
	    Enumeration zipFileEntries = zip.entries();
	
	    // Process each entry
	    while (zipFileEntries.hasMoreElements()) {
	        // grab a zip file entry
	        ZipEntry entry = (ZipEntry) zipFileEntries.nextElement();
	        String currentEntry = entry.getName();
	        File destFile = new File(newPath, currentEntry);
	        ret = destFile
	        File destinationParent = destFile.getParentFile();
	
	        // create the parent directory structure if needed
	        if (!destinationParent.exists()) {
	        	destinationParent.mkdirs();
	        }
			
	        if (!entry.isDirectory()) {
	            BufferedInputStream bis = new BufferedInputStream(zip.getInputStream(entry));
	            int currentByte;
	            // establish buffer for writing file
	            byte[] data = new byte[BUFFER];
	
	            // write the current file to disk
	            FileOutputStream fos = new FileOutputStream(destFile);
	            BufferedOutputStream dest = new BufferedOutputStream(fos,
	            BUFFER);
	
	            // read and write until last byte is encountered
	            while ((currentByte = bis.read(data, 0, BUFFER)) != -1) {
	                dest.write(data, 0, currentByte);
	            }
	            dest.flush();
	            dest.close();
	            bis.close();
	        } else {
	        	destFile.mkdirs();
	        }
	
	        if (currentEntry.endsWith(".zip")) {
	            // found a zip file, try to open
	            extractFolder(destFile.getAbsolutePath());
	        }
	    }
	    return ret
	}
	
	public static void deleteFolder(File folder) {
	    File[] files = folder.listFiles();
	    if(files != null) {
	        for(File f: files) {
	            if(f.isDirectory()) {
	                deleteFolder(f);
	            } else {
	                f.delete();
	            }
	        }
	    }
	    folder.delete();
	}
	
	public class Downloader {
		String site = ""
		String filename = ""
		WhiteboxGui pluginHost;
		
		public Downloader(WhiteboxGui pluginHost, String filename, String site) {
			this.pluginHost = pluginHost
			this.site = site
			this.filename = filename
		}
		
	    public int download(String fn) throws Exception {
	    
		    try {
		    	
		        URL url=new URL(site);
		        HttpURLConnection connection = (HttpURLConnection) url.openConnection();
		        connection.setReadTimeout(1000000)
		        int filesize = connection.getContentLength();
		        float totalDataRead = 0;
	            java.io.BufferedInputStream inStream = new java.io.BufferedInputStream(connection.getInputStream());
	            java.io.FileOutputStream fos = new java.io.FileOutputStream(filename);
	            java.io.BufferedOutputStream bout = new BufferedOutputStream(fos,1024);
	            byte[] data = new byte[1024];
	            int i = 0;
	            while((i = inStream.read(data, 0, 1024)) >= 0) {
		            totalDataRead=totalDataRead+i;
		            bout.write(data,0,i);
		            float percent = (totalDataRead * 100f) / filesize;
		            pluginHost.updateProgress("Downloading ${fn}...", (int) percent)
		            // check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return false
					}
	            }
	            bout.close();
	            inStream.close();
	            pluginHost.updateProgress("Progress:", 0)
	            return 0
		    } catch (FileNotFoundException e) {
			    return 1
		    } catch(SocketTimeoutException e) {
			    return 2
		    } catch(IOException e) {
			    return 3
		    } catch(Exception e) {
		    	pluginHost.logException("Error in " + descriptiveName, e)
			    return 4
		    }
		}
	}
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new RetrieveSRTMData(pluginHost, args, name, descriptiveName)
}
