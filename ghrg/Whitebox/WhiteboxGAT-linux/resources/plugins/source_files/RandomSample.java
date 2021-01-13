/*
 * Copyright (C) 2011-2012 Dr. John Lindsay <jlindsay@uoguelph.ca>
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
package plugins;

import java.util.Date;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Random;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;


/**
 * This tool can be used to create a random sample of grid cells. 
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class RandomSample implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;

    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name
     * containing no spaces.
     *
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "RandomSample";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Random Sample";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Creates an image containing randomly located sample grid cells.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "StatisticalTools" };
    	return ret;
    }

    /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the
     * class that the plugin will send all feedback messages, progress updates,
     * and return objects.
     *
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }

    /**
     * Used to communicate feedback pop-up messages between a plugin tool and
     * the main Whitebox user-interface.
     *
     * @param feedback String containing the text to display.
     */
    private void showFeedback(String message) {
        if (myHost != null) {
            myHost.showFeedback(message);
        } else {
            System.out.println(message);
        }
    }

    /**
     * Used to communicate a return object from a plugin tool to the main
     * Whitebox user-interface.
     *
     * @return Object, such as an output WhiteboxRaster.
     */
    private void returnData(Object ret) {
        if (myHost != null) {
            myHost.returnData(ret);
        }
    }

    private int previousProgress = 0;
    private String previousProgressLabel = "";
   
    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progressLabel A String to use for the progress label.
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(String progressLabel, int progress) {
        if (myHost != null && ((progress != previousProgress) || 
                (!progressLabel.equals(previousProgressLabel)))) {
            myHost.updateProgress(progressLabel, progress);
        }
        previousProgress = progress;
        previousProgressLabel = progressLabel;
    }

    /**
     * Used to communicate a progress update between a plugin tool and the main
     * Whitebox user interface.
     *
     * @param progress Float containing the progress value (between 0 and 100).
     */
    private void updateProgress(int progress) {
        if (myHost != null && progress != previousProgress) {
            myHost.updateProgress(progress);
        }
        previousProgress = progress;
    }
    
     /**
     * Sets the arguments (parameters) used by the plugin.
     *
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    
    private boolean cancelOp = false;
  
    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
     *
     * @param cancel Set to true if the plugin should be canceled.
     */
    @Override
    public void setCancelOp(boolean cancel) {
        cancelOp = cancel;
    }
    
    private void cancelOperation() {
        showFeedback("Operation cancelled.");
        updateProgress("Progress: ", 0);
    }
    
    private boolean amIActive = false;
   
    /**
     * Used by the Whitebox GUI to tell if this plugin is still running.
     *
     * @return a boolean describing whether or not the plugin is actively being
     * used.
     */
    @Override
    public boolean isActive() {
        return amIActive;
    }

    /**
     * Used to execute this plugin tool.
     */
    @Override
    public void run() {
        amIActive = true;
        
        String inputHeader = null;
        String outputHeader = null;
    	
        WhiteboxRaster image;
        WhiteboxRaster output;
        int cols, rows;
        int progress = 0;
        int col, row;
        int i;
        int numSamplePoints = 0;
                
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        inputHeader = args[0];
        outputHeader = args[1];
        numSamplePoints = Integer.parseInt(args[2]);
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader == null) || (outputHeader == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            image = new WhiteboxRaster(inputHeader, "r");
            rows = image.getNumberRows();
            cols = image.getNumberColumns();
            
            if (rows * cols < numSamplePoints) {
                showFeedback("The number of samples cannot exceed the number of cells.");
                return;
            }

            output = new WhiteboxRaster(outputHeader, "rw", inputHeader, WhiteboxRaster.DataType.FLOAT, 0);
            output.setPreferredPalette("qual.pal");
            output.setDataScale(WhiteboxRaster.DataScale.CATEGORICAL);

            image.close();
            
            // A priority gueue is used so that we can later add each of the
            // points to the grid sorted by row value. Otherwise we will end up
            // reading and writing to disc a lot as we enter values to the grid
            // in random locations.
            NonDuplicatingPriorityQueue queue = new NonDuplicatingPriorityQueue(numSamplePoints);
            Random generator = new Random();
            GridCell gc;
            i = 0;
            do {
                row = generator.nextInt(rows);
                col = generator.nextInt(cols);
                gc = new GridCell(row, col);
                if (queue.add(gc)) {
                   i++; 
                   progress = (int)(100f * i / numSamplePoints);
                   updateProgress("Loop 1 of 2:", progress);
                }
            } while (i < numSamplePoints);
            
            Iterator<GridCell> it = queue.iterator();
            i = 1;
            do  {
                // Returns an iterator over the elements in this queue.
                gc = queue.poll();
                output.setValue(gc.row, gc.col, i);
                i++;
                progress = (int)(100f * i / numSamplePoints);
                updateProgress("Loop 2 of 2:", progress);
            } while (i < numSamplePoints);
            
            
            output.addMetadataEntry("Created by the "
                    + getDescriptiveName() + " tool.");
            output.addMetadataEntry("Created on " + new Date());
            
            output.close();

            // returning a header file string displays the image.
            returnData(outputHeader);
            
            
        } catch (OutOfMemoryError oe) {
            myHost.showFeedback("An out-of-memory error has occurred during operation.");
        } catch (Exception e) {
            myHost.showFeedback("An error has occurred during operation. See log file for details.");
            myHost.logException("Error in " + getDescriptiveName(), e);
        } finally {
            updateProgress("Progress: ", 0);
            // tells the main application that this process is completed.
            amIActive = false;
            myHost.pluginComplete();
        }
    }
    
    class NonDuplicatingPriorityQueue extends PriorityQueue<GridCell> {
        public NonDuplicatingPriorityQueue(int initialCapacity) {
            super(initialCapacity);
        }
        
        @Override
        public boolean add(GridCell gc) {
            try {
                // This will only add the gridcell if it's not already in the queue.
                Iterator<GridCell> i = this.iterator();
                while (i.hasNext()) {
                    // Returns an iterator over the elements in this queue.
                    GridCell gc2 = i.next();
                    if (gc.compareTo(gc2) == 0) {
                        return false;
                    }
                }
                super.add(gc);
                return true;
            } catch (Exception e) {
                return false;
            }
        }
    }
    
    class GridCell implements Comparable<GridCell> {

        public int row;
        public int col;
        
        public GridCell(int row, int col) {
            this.row = row;
            this.col = col;
        }

        @Override
        public int compareTo(GridCell cell) {
            final int BEFORE = -1;
            final int EQUAL = 0;
            final int AFTER = 1;

            if (this.row < cell.row) {
                return BEFORE;
            } else if (this.row > cell.row) {
                return AFTER;
            }

            if (this.col < cell.col) {
                return BEFORE;
            } else if (this.col > cell.col) {
                return AFTER;
            }

            return EQUAL;
        }
    }
    
//    // This method is only used during testing.
//    public static void main(String[] args) {
//        args = new String[3];
//        args[0] = "/Users/johnlindsay/Documents/Data/Waterloo streams.dep";
//        args[1] = "/Users/johnlindsay/Documents/Data/tmp2.dep";
//        args[2] = "1000";
//        
//        RandomSample rs = new RandomSample();
//        rs.setArgs(args);
//        rs.run();
//    }
}
