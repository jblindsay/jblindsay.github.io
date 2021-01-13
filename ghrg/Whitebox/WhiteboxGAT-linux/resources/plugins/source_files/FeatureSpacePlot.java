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

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.geom.AffineTransform;
import java.text.DecimalFormat;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.io.*;
import java.awt.print.*;
import javax.print.attribute.*;
import whitebox.structures.ExtensionFileFilter;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import java.util.ArrayList;

/**
 * A feature space plot is similar to a scattergram, it plots the frequency of occurrance of brightness value pairs in the two-dimensional feature space of defined by two images.
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class FeatureSpacePlot implements WhiteboxPlugin {
    
    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "FeatureSpacePlot";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
    	return "Feature Space Plot";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
    	return "Creates a feature space plot for two multispectral bands.";
    }
     /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
    	String[] ret = { "ImageClass", "StatisticalTools" };
    	return ret;
    }
     /**
     * Sets the WhiteboxPluginHost to which the plugin tool is tied. This is the class
     * that the plugin will send all feedback messages, progress updates, and return objects.
     * @param host The WhiteboxPluginHost that called the plugin tool.
     */  
    @Override
    public void setPluginHost(WhiteboxPluginHost host) {
        myHost = host;
    }
    /**
     * Used to communicate feedback pop-up messages between a plugin tool and the main Whitebox user-interface.
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
     * Used to communicate a return object from a plugin tool to the main Whitebox user-interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * Used to communicate a progress update between a plugin tool and the main Whitebox user interface.
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
     * @param args An array of string arguments.
     */
    @Override
    public void setArgs(String[] args) {
        this.args = args.clone();
    }
    
    private boolean cancelOp = false;
    /**
     * Used to communicate a cancel operation from the Whitebox GUI.
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
     * @return a boolean describing whether or not the plugin is actively being used.
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
        
        String inputHeader1 = null;
        String inputHeader2 = null;
        int row, col;
        double z;
        float progress = 0;
        int m, n;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }
        
        inputHeader1 = args[0];
        inputHeader2 = args[1];
        
        
        // check to see that the inputHeader and outputHeader are not null.
        if ((inputHeader1 == null) || (inputHeader2 == null)) {
            showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        try {
            WhiteboxRasterInfo image1 = new WhiteboxRasterInfo(inputHeader1);
            int rows = image1.getNumberRows();
            int cols = image1.getNumberColumns();
            double noData1 = image1.getNoDataValue();

            WhiteboxRasterInfo image2 = new WhiteboxRasterInfo(inputHeader2);
            if (image2.getNumberRows() != rows || image2.getNumberColumns() != cols) {
                showFeedback("The input images must have the same number of rows and columns");
                return;
            }
            double noData2 = image2.getNoDataValue();
            
            int image1Min = (int)image1.getMinimumValue();
            int image2Min = (int)image2.getMinimumValue();
            int image1Max = (int)image1.getMaximumValue();
            int image2Max = (int)image2.getMaximumValue();
            int image1Range = image1Max - image1Min + 1;
            int image2Range = image2Max - image2Min + 1;
            
            int[][] featureSpace = new int[image1Range][image2Range];
            
            double data1[] = null;
            double data2[] = null;
            for (row = 0; row < rows; row++) {
                data1 = image1.getRowValues(row);
                data2 = image2.getRowValues(row);
                for (col = 0; col < cols; col++) {
                    if (data1[col] != noData1 && data2[col] != noData2) {
                        m = (int)(data1[col] - image1Min);
                        n = (int)(data2[col] - image2Min);
                        featureSpace[m][n]++;
                    }

                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (float) (100f * row / (rows - 1));
                updateProgress((int) progress);
            }

            
            Plot plot = new Plot(featureSpace, image1.getShortHeaderFile(), 
                    image2.getShortHeaderFile(), image1Min, image1Max, image2Min,
                    image2Max, myHost.getDefaultFont());

            image1.close();
            image2.close();
            
            // returning a header file string displays the image.
            returnData(plot);

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
    
    
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        FeatureSpacePlot fsp = new FeatureSpacePlot();
//        args = new String[2];
//        args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band1 clipped.dep";
//        args[1] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band2 clipped.dep";
//        
//        fsp.setArgs(args);
//        fsp.run();
//        
//    }
}

class Plot extends JPanel implements ActionListener, Printable, MouseMotionListener, MouseListener {
    private boolean plotCreated = false;
    private String shortName1, shortName2;
    private int bottomMargin = 75;
    private int topMargin = 45;
    private int leftMargin = 90;
    private int rightMargin = 20;
    private int[][] featureSpace;     
    private int[] imageData;
    private int width;
    private int height;
    private int numPix;
    private int image1Min, image2Min, image1Max, image2Max;
    private boolean isHillshaded = false;
    private JPopupMenu myPopup = null;
    private JCheckBoxMenuItem cmi = null;
    private JCheckBoxMenuItem plotPeaksMi = null;
    private JMenuItem backgroundMi = null;
    private JMenuItem gridMi = null;
    private JMenuItem paletteMi = null;
    private int backgroundColor = 255;
    private boolean gridOn = true;
    private boolean spectrumPaletteMode = true;
    private Font myFont;
    // Constructors
    public Plot(int[][] featureSpace, String shortName1, String shortName2, int image1Min,
            int image1Max, int image2Min, int image2Max, Font myFont) {
        this.featureSpace = featureSpace.clone();
        this.shortName1 = shortName1;
        this.shortName2 = shortName2;
        this.image1Min = image1Min;
        this.image2Min = image2Min;
        this.image1Max = image1Max;
        this.image2Max = image2Max;
        this.myFont = myFont;
        
        setMouseMotionListener();
        setMouseListener();
        setUp();
    }
    
    // Methods
    private void setUp() {
        try {
            width = featureSpace[0].length;
            height = featureSpace.length;
            numPix = width * height;
            createPopupMenus();
        } catch (Exception e) {
            
        }
    }
    
    //java.util.ArrayList<Dimension> peaks = new java.util.ArrayList<Dimension>();
    boolean plotPeaks = false;
    boolean peaksFound = false;
    int[][] peaks = null;
    private void findLocalPeaks() {
        // local peaks are higher than all their neighbours, and are relatively high in the landscape
        double z;
        boolean flag = false;
        int rN, cN;
        double peakSize = 0.00001;
        int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
        int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
        double[] N = new double[8];
        int numPeaks = 0;
        
        // apply a 3x3 mean filter to remove insignificant local peaks
        double[][] featureSpace2 = new double[height][width];
        double cellTotal = 0;
        int numNeighbours = 0;
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                if (featureSpace[r][c] > 0) {
                    cellTotal = featureSpace2[r][c];
                    numNeighbours = 1;
                    for (int j = 0; j < 8; j++) {
                        rN = r + Dy[j];
                        cN = c + Dx[j];
                        if (rN >= 0 && rN < height && cN >= 0 && cN < width) {
                            if (featureSpace[rN][cN] != 0) {
                                cellTotal += featureSpace[rN][cN];
                                numNeighbours++;
                            }
                        }
                    }
                    featureSpace2[r][c] = cellTotal / numNeighbours;
                }
            }
        }
        
        long total = 0;
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                total += featureSpace[r][c];
            }
        }
        
        // first scan to find how many peaks there are
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                z = featureSpace2[r][c];
                if (z / total > peakSize) {
                    flag = true;
                    for (int j = 0; j < 8; j++) {
                        rN = r + Dy[j];
                        cN = c + Dx[j];
                        if (rN >= 0 && rN < height && cN >= 0 && cN < width) {
                            if (featureSpace2[rN][cN] > z) {
                                flag = false;
                                break;
                            }
                        }
                    }
                    if (flag) {
                        // it's a peak
                        numPeaks++;
                    }
                }
            }
        }
        
        peaks = new int[numPeaks][3];
        
        int k = 0;
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                z = featureSpace2[r][c];
                if (z / total > peakSize) {
                    flag = true;
                    for (int j = 0; j < 8; j++) {
                        rN = r + Dy[j];
                        cN =c + Dx[j];
                        if (rN >= 0 && rN < height && cN >= 0 && cN < width) {
                            if (featureSpace2[rN][cN] > z) {
                                flag = false;
                                break;
                            }
                        }
                    }
                    if (flag) {
                        // it's a peak
                        peaks[k][0] = r + image1Min;
                        peaks[k][1] = c + image2Min;
                        peaks[k][2] = featureSpace[r][c];
                        k++;
                        //peaks.add(new Dimension(r + image1Min, c + image2Min));
                    }
                }
            }
        }
    }
    
    private void createPlotImageData() {
        // find the maximum point density in featureSpace
        int maxDensityVal = -1;
        long total = 0;
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                if (featureSpace[r][c] > maxDensityVal) {
                    maxDensityVal = featureSpace[r][c];
                }
                total += featureSpace[r][c];
            }
        }
        double maxDensityValLog = Math.log(maxDensityVal);
        imageData = new int[numPix];
        int index = 0;
        int red, green, blue;
        double loc;
        double paletteZoneSize = 1d / 6; 
        if (isHillshaded) {
            double[][] hillshade = new double[height][width];
            int[] Dy = {-1, 0, 1, 1, 1, 0, -1, -1};
            int[] Dx = {1, 1, 1, 0, -1, -1, -1, 0};
            double[] N = new double[8];
            double term1, term2, term3;
            final double radToDeg = 180 / Math.PI;
            final double degToRad = Math.PI / 180;
            double azimuth = 135 * degToRad;
            double altitude = 45 * degToRad;
            double z;
            double fx, fy, aspect;
            double Rad180 = 180 * degToRad;
            double Rad90 = 90 * degToRad;
            int rN, cN;
            double sinTheta;
            double cosTheta;
            double tanSlope;
            sinTheta = Math.sin(altitude);
            cosTheta = Math.cos(altitude);


            double minHill = 99999;
            double maxHill = -99999;
            for (int r = height - 1; r >= 0; r--) {
                for (int c = 0; c < width; c++) {
                    z = (double) (featureSpace[r][c] / (double)total * 1000);
                    for (int j = 0; j < 8; j++) {
                        rN = r + Dy[j];
                        cN = c + Dx[j];
                        if (rN >= 0 && cN >= 0 && rN < height && cN < width) {
                            N[j] = (featureSpace[rN][cN] / (double)total) * 1000;
                        } else {
                            N[j] = z;
                        }
                    }
                    // calculate slope and aspect
                    fy = (N[6] - N[4] + 2 * (N[7] - N[3]) + N[0] - N[2]) / 8;
                    fx = (N[2] - N[4] + 2 * (N[1] - N[5]) + N[0] - N[6]) / 8;
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
                    //if (z > 1) { z = 1; }
                    //if (z < 0.3) { z = 0.3; }
                    if (z < minHill) { minHill = z; }
                    if (z > maxHill) { maxHill = z; }
                    hillshade[r][c] = z;
                }
            }
            
            double range = maxHill - minHill;
            for (int r = height - 1; r >= 0; r--) {
                for (int c = 0; c < width; c++) {
                    z = (hillshade[r][c] - minHill) / (range);
                    if (z < 0.3) { z = 0.3; }
                    hillshade[r][c] = z;
                }
            }
            
            for (int r = height - 1; r >= 0; r--) {
                for (int c = 0; c < width; c++) {
                    //loc = (double) featureSpace[r][c] / (double) maxDensityVal;
                    loc = Math.log(featureSpace[r][c]) / maxDensityValLog;
                    if (spectrumPaletteMode) {
                        if (featureSpace[r][c] == 0) {
                            red = backgroundColor;
                            green = backgroundColor;
                            blue = backgroundColor;
                        } else if (loc < paletteZoneSize) {
                            red = (int) (128 - (loc / paletteZoneSize * 128));
                            green = 0;
                            blue = (int) (128 + loc / paletteZoneSize * 128);
                        } else if (loc < 2 * paletteZoneSize) {
                            red = 0;
                            green = (int) ((loc - paletteZoneSize) / paletteZoneSize * 255);
                            blue = 255;
                        } else if (loc < 3 * paletteZoneSize) {
                            red = 0;
                            green = 255;
                            blue = (int) (255 - ((loc - 2 * paletteZoneSize) / paletteZoneSize * 255));
                        } else if (loc < 4 * paletteZoneSize) {
                            red = (int) ((loc - 3 * paletteZoneSize) / paletteZoneSize * 255);
                            green = 255;
                            blue = 0;
                        } else if (loc < 5 * paletteZoneSize) {
                            red = 255;
                            green = (int) (255 - ((loc - 4 * paletteZoneSize) / paletteZoneSize * 255));
                            blue = 0;
                        } else {
                            red = (int) (255 - ((loc - 5 * paletteZoneSize) / paletteZoneSize * 180));
                            green = 0;
                            blue = 0;
                        }
                        if (featureSpace[r][c] != 0) {
                            red = (int) (red * hillshade[r][c]);
                            green = (int) (green * hillshade[r][c]);
                            blue = (int) (blue * hillshade[r][c]);
                        }
                    } else {
                        if (featureSpace[r][c] == 0) {
                            red = 0;
                            green = 0;
                            blue = 0;
                        } else {
                            red = (int)(255 * loc);
                            green = (int)(255 * loc);
                            blue = (int)(255 * loc);
                        }
                    }
                    imageData[index] = ((255 << 24) | (red << 16) | (green << 8) | blue);
                    index++;
                }
            }
        } else {
            for (int r = height - 1; r >= 0; r--) {
                for (int c = 0; c < width; c++) {
                    loc = Math.log(featureSpace[r][c]) / maxDensityValLog;
                    if (spectrumPaletteMode) {
                        if (featureSpace[r][c] == 0) {
                            red = backgroundColor;
                            green = backgroundColor;
                            blue = backgroundColor;
                        } else if (loc < paletteZoneSize) {
                            red = (int) (128 - (loc / paletteZoneSize * 128));
                            green = 0;
                            blue = (int) (128 + loc / paletteZoneSize * 128);
                        } else if (loc < 2 * paletteZoneSize) {
                            red = 0;
                            green = (int) ((loc - paletteZoneSize) / paletteZoneSize * 255);
                            blue = 255;
                        } else if (loc < 3 * paletteZoneSize) {
                            red = 0;
                            green = 255;
                            blue = (int) (255 - ((loc - 2 * paletteZoneSize) / paletteZoneSize * 255));
                        } else if (loc < 4 * paletteZoneSize) {
                            red = (int) ((loc - 3 * paletteZoneSize) / paletteZoneSize * 255);
                            green = 255;
                            blue = 0;
                        } else if (loc < 5 * paletteZoneSize) {
                            red = 255;
                            green = (int) (255 - ((loc - 4 * paletteZoneSize) / paletteZoneSize * 255));
                            blue = 0;
                        } else {
                            red = (int) (255 - ((loc - 5 * paletteZoneSize) / paletteZoneSize * 128));
                            green = 0;
                            blue = 0;
                        }
                    } else {
                        if (featureSpace[r][c] == 0) {
                            red = 0;
                            green = 0;
                            blue = 0;
                        } else {
                            red = (int)(255 * loc);
                            green = (int)(255 * loc);
                            blue = (int)(255 * loc);
                        }
                    }
                    imageData[index] = ((255 << 24) | (red << 16) | (green << 8) | blue);
                    index++;
                }
            }
        }
        plotCreated = true;
    }
    
    private void setMouseMotionListener() {
        this.addMouseMotionListener(this);
    }
    
    private void setMouseListener() {
        this.addMouseListener(this);
    }
    
    public void refresh() {
        plotCreated = false;
        repaint();
    }
    
    @Override
    public void paint (Graphics g) {
        drawPlot(g);
    }
    
    private void drawPlot(Graphics g) {
        if (!plotCreated) { createPlotImageData(); }
        
        leftMargin = 60;
        Graphics2D g2d = (Graphics2D) g;
        
        g2d.setFont(myFont);
        
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, getWidth(), getHeight());
        
        double activeWidth = getWidth() - leftMargin - rightMargin;
        double activeHeight = getHeight() - topMargin - bottomMargin;
        int bottomY = getHeight() - bottomMargin;
        int rightX = getWidth() - rightMargin;
        int rectLeft, rectWidth, rectTop, rectHeight;
        
        Image image = createImage(new MemoryImageSource(width, height, imageData, 0, width));
        g2d.drawImage(image, leftMargin, topMargin, (int)activeWidth, (int)activeHeight, this);

        // draw axes
        g2d.setColor(Color.black);
        g2d.drawLine(leftMargin, bottomY, rightX, bottomY);
        g2d.drawLine(leftMargin, bottomY, leftMargin, topMargin);
        g2d.drawLine(leftMargin, topMargin, rightX, topMargin);
        g2d.drawLine(rightX, bottomY, rightX, topMargin);
        
        // draw ticks
        int tickSize = 4;
        g2d.drawLine(leftMargin, bottomY, leftMargin, bottomY + tickSize);
        g2d.drawLine(rightX, bottomY, rightX, bottomY + tickSize);
        g2d.drawLine(leftMargin, bottomY, leftMargin - tickSize, bottomY);
        g2d.drawLine(leftMargin, topMargin, leftMargin - tickSize, topMargin);
        
        // histo labels
        DecimalFormat df = new DecimalFormat("#,###,###.###");
        Font font = new Font("SanSerif", Font.PLAIN, 11);
        FontMetrics metrics = g.getFontMetrics(font);
        int hgt, adv;
        hgt = metrics.getHeight();
        
        // x-axis labels
        String label;
        label = df.format(image2Min);
        g2d.drawString(label, leftMargin, bottomY + hgt + 4);
        label = df.format(image2Max);
        adv = metrics.stringWidth(label);
        g2d.drawString(label, rightX - adv, bottomY + hgt + 4);
        label = shortName2;
        adv = metrics.stringWidth(label);
        int xAxisMidPoint = (int)(leftMargin + activeWidth / 2);
        g2d.drawString(label, xAxisMidPoint - adv / 2, bottomY + 2 * hgt + 4);
        
        // y-axis labels
        
        // rotate the font
        Font oldFont = g.getFont();
        Font f = oldFont.deriveFont(AffineTransform.getRotateInstance(-Math.PI / 2.0));
        g2d.setFont(f);
        
        int yAxisMidPoint = (int)(topMargin + activeHeight / 2);
        int offset;
        label = shortName1;
        offset = metrics.stringWidth("0") + 12 + hgt;
        adv = metrics.stringWidth(label);
        g2d.drawString(label, leftMargin - offset, yAxisMidPoint + adv / 2);
        
        // replace the rotated font.
        g2d.setFont(oldFont);

        df = new DecimalFormat("0");
        label = df.format(image1Min);
        adv = metrics.stringWidth(label);
        g2d.drawString(label, leftMargin - adv - 12, bottomY + hgt / 2);
        label = df.format(image1Max);
        adv = metrics.stringWidth(label);
        g2d.drawString(label, leftMargin - adv - 12, topMargin + hgt / 2);
        
        // title
        
        // bold font
        oldFont = g.getFont();
        font = font = new Font("SanSerif", Font.BOLD, 12);
        g2d.setFont(font);
        
        label = "Feature Space Plot: " + shortName1 + " vs. " + shortName2;
        adv = metrics.stringWidth(label);
        g2d.drawString(label, getWidth() / 2 - adv / 2, topMargin - hgt - 5);
        
        g2d.setFont(oldFont);
        
        if (gridOn) {
            // vertical grid lines
            int gridSpacing = 50;
            if (image2Max > 255) {
                gridSpacing = 100;
            }
            int loc = 0;
            if (backgroundColor == 0 || !spectrumPaletteMode) {
                g2d.setColor(Color.white);
            } else {
                g2d.setColor(Color.black);
            }
            for (int i = image2Min; i <= image2Max; i++) {
                if (i % gridSpacing == 0) {
                    loc = (int)(leftMargin + (i - image2Min) / ((double) image2Max - image2Min) * activeWidth);

                    g2d.drawLine(loc, bottomY, loc, topMargin - 1);
                }
            }

            // horizontal grid lines
            if (image1Max > 255) {
                gridSpacing = 100;
            }
            for (int i = image1Min; i <= image1Max; i++) {
                if (i % gridSpacing == 0) {
                    loc = (int)(bottomY - (i - image1Min) / ((double) image1Max - image1Min) * activeHeight);

                    g2d.drawLine(leftMargin + 1, loc, rightX, loc);
                }
            }
        }
        
        if (plotPeaks) {
            if (!peaksFound) {
                findLocalPeaks();
            }
            g2d.setColor(Color.black);
            int crossHairSize = 4;
            int row, col;
            for (int j = 0; j < peaks.length; j++) {
                row = (int)(bottomY - (peaks[j][0] - image1Min) / ((double) image1Max - image1Min) * activeHeight);
                col = (int)(leftMargin + (peaks[j][1] - image2Min) / ((double) image2Max - image2Min) * activeWidth);
                g2d.drawLine(col - crossHairSize, row, col + crossHairSize, row);
                g2d.drawLine(col, row - crossHairSize, col, row + crossHairSize);
            }
         }
        
        if (drawPositionLine) {
            // draw a red vertical line at this position.
            g2d.setColor(Color.red);
            g2d.drawLine(posX, bottomY, posX, topMargin);
            g2d.setColor(Color.white);
            g2d.drawLine(posX + 1, bottomY, posX + 1, topMargin);
            
            // draw a red horizontal line at this position.
            g2d.setColor(Color.red);
            g2d.drawLine(leftMargin, posY, rightX, posY);
            g2d.setColor(Color.white);
            g2d.drawLine(leftMargin, posY + 1, rightX, posY + 1);
            
            // which bin is it and how many are in the bin?
            
            int x = (int)(((double)posX - leftMargin) / activeWidth * (width - 1));
            int y = (int)(((1 - ((double)posY - topMargin) / activeHeight) * (height - 1)));
            df = new DecimalFormat("0");
            int val1 = featureSpace[y][x];
            String val = df.format(val1);
            label = "x: " + (x + image2Min) + " y: " + (y + image1Min) + " Value: " + val;
            adv = metrics.stringWidth(label);
            g2d.setColor(Color.black); 
            g2d.drawString(label, 10, getHeight() - hgt - 10);
        
        }
    }
    
    public boolean saveToImage(String fileName) {
        try {
            int width = (int)this.getWidth();
            int height =(int)this.getHeight();
            // TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
            // into integer pixels
            BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

            Graphics ig = bi.createGraphics();
            drawPlot(ig);
            int i = fileName.lastIndexOf(".");
            String extension = fileName.substring(fileName.lastIndexOf(".") + 1).toUpperCase();
            if (!ImageIO.write(bi, extension, new File(fileName))) {
                return false;
            }
            return true;
        } catch (Exception ex) {
            return false;
        }
    }
    
    @Override
    public int print(Graphics g, PageFormat pf, int page)
            throws PrinterException {
        if (page > 0) {
            return NO_SUCH_PAGE;
        }

        int i = pf.getOrientation();
        
        // get the size of the page
        double pageWidth = pf.getImageableWidth();
        double pageHeight = pf.getImageableHeight();
        double myWidth = this.getWidth();// - borderWidth * 2;
        double myHeight = this.getHeight();// - borderWidth * 2;
        double scaleX = pageWidth / myWidth;
        double scaleY = pageHeight / myHeight;
        double minScale = Math.min(scaleX, scaleY);
                
        Graphics2D g2d = (Graphics2D) g;
        g2d.translate(pf.getImageableX(), pf.getImageableY());
        g2d.scale(minScale, minScale);
        
        drawPlot(g);

        return PAGE_EXISTS;
    }
    
    private boolean drawPositionLine = false;
    private int posX = 0;
    private int posY = 0;
    
    @Override
    public void mouseMoved(MouseEvent e) {
        // is the mouse over the chart area?
        int x = e.getX();
        int y = e.getY();
        int width = getWidth();
        int height = getHeight();
        if (x >= leftMargin && x <= (width - rightMargin) &&
                y >= topMargin && y <= (height - bottomMargin) ) {
            drawPositionLine = true;
            posX = x;
            posY = y;
        } else {
            drawPositionLine = false;
        }
        repaint();
    }
    
    private void createPopupMenus() {
            // menu
            myPopup = new JPopupMenu();

            JMenuItem mi = new JMenuItem("Save");
            mi.addActionListener(this);
            mi.setActionCommand("save");
            myPopup.add(mi);
            
            mi = new JMenuItem("Print");
            mi.addActionListener(this);
            mi.setActionCommand("print");
            myPopup.add(mi);
            
            myPopup.addSeparator();
            
            cmi = new JCheckBoxMenuItem("Hillshade Plot");
            cmi.addActionListener(this);
            cmi.setState(isHillshaded);
            cmi.setActionCommand("hillshade plot");
            myPopup.add(cmi);
            
            backgroundMi = new JMenuItem("Black Background");
            backgroundMi.addActionListener(this);
            backgroundMi.setActionCommand("changeBackground");
            myPopup.add(backgroundMi);
            
            paletteMi = new JMenuItem("Grey-scale Palette");
            paletteMi.addActionListener(this);
            paletteMi.setActionCommand("changePalette");
            myPopup.add(paletteMi);
            
            gridMi = new JMenuItem("Turn Off Grid");
            gridMi.addActionListener(this);
            gridMi.setActionCommand("grid");
            myPopup.add(gridMi);
            
            plotPeaksMi = new JCheckBoxMenuItem("Plot Peaks");
            plotPeaksMi.addActionListener(this);
            plotPeaksMi.setState(plotPeaks);
            plotPeaksMi.setActionCommand("plot peaks");
            myPopup.add(plotPeaksMi);
            
            
            myPopup.setOpaque(true);
            myPopup.setLightWeightPopupEnabled(true);


    }
    
    private void printPlot() {
        PrinterJob job = PrinterJob.getPrinterJob();
        PrintRequestAttributeSet aset = new HashPrintRequestAttributeSet();
        //PageFormat pf = job.pageDialog(aset);
        job.setPrintable(this);
        boolean ok = job.printDialog(aset);
        if (ok) {
            try {
                job.print(aset);
            } catch (PrinterException ex) {
                //showFeedback("An error was encountered while printing." + ex);
                /* The job did not successfully complete */
            }
        }
    }
    
    private void savePlotAsImage() {
        // get the possible image name.
        String imageName = shortName1 + "_vs_" + shortName2;

        // Ask the user to specify a file name for saving the histo.
        String pathSep = File.separator;
        JFileChooser fc = new JFileChooser();
        fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        //fc.setCurrentDirectory(new File(workingDirectory + pathSep + imageName + ".png"));
        fc.setAcceptAllFileFilterUsed(false);

        //File f = new File(workingDirectory + pathSep + imageName + ".png");
        //fc.setSelectedFile(f);
                
        // set the filter.
        ArrayList<ExtensionFileFilter> filters = new ArrayList<ExtensionFileFilter>();
        String[] extensions = ImageIO.getReaderFormatNames(); //{"PNG", "JPEG", "JPG"};
        String filterDescription = "Image Files (" + extensions[0];
        for (int i = 1; i < extensions.length; i++) {
            filterDescription += ", " + extensions[i];
        }
        filterDescription += ")";
        ExtensionFileFilter eff = new ExtensionFileFilter(filterDescription, extensions);
        fc.setFileFilter(eff);

        int result = fc.showSaveDialog(this);
        File file = null;
        if (result == JFileChooser.APPROVE_OPTION) {
            file = fc.getSelectedFile();
            // see if file has an extension.
            if (file.toString().lastIndexOf(".") <= 0) {
                String fileName = file.toString() + ".png";
                file = new File(fileName);
            }
            
            String fileDirectory = file.getParentFile() + pathSep;
//            if (!fileDirectory.equals(workingDirectory)) {
//                workingDirectory = fileDirectory;
//            }

            // see if the file exists already, and if so, should it be overwritten?
            if (file.exists()) {
                Object[] options = {"Yes", "No"};
                int n = JOptionPane.showOptionDialog(this,
                        "The file already exists.\n"
                        + "Would you like to overwrite it?",
                        "Whitebox GAT Message",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null, //do not use a custom Icon
                        options, //the titles of buttons
                        options[0]); //default button title

                if (n == JOptionPane.YES_OPTION) {
                    file.delete();
                } else if (n == JOptionPane.NO_OPTION) {
                    return;
                }
            }
            if (!saveToImage(file.toString())) {
//                showFeedback("An error occurred while saving the map to the image file.");
            }
        }
        
    }
    
    @Override
    public void mouseDragged(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
    
    // ActionListener for events
    @Override
    public void actionPerformed(ActionEvent ae) {
        String actionCommand = ae.getActionCommand().toLowerCase();
        if (actionCommand.equals("save")) {
            savePlotAsImage();
        } else if (actionCommand.equals("print")) {
            printPlot();
        } else if (actionCommand.equals("hillshade plot")) {
            isHillshaded = !isHillshaded;
            refresh();
        } else if (actionCommand.equals("plot peaks")) {
            plotPeaks = !plotPeaks;
            refresh();
        } else if (actionCommand.equals("changebackground")) {
            if (backgroundMi.getText().equals("White Background")) {
                backgroundMi.setText("Black Background");
                backgroundColor = 255;
            } else {
                backgroundMi.setText("White Background");
                backgroundColor = 0;
            }
            refresh();
        } else if (actionCommand.equals("grid")) {
            if (gridOn) {
                gridOn = false;
                gridMi.setText("Turn On Grid");
            } else if (!gridOn) {
                gridOn = true;
                gridMi.setText("Turn Off Grid");
            }
            refresh();
        } else if (actionCommand.equals("changepalette")) {
            spectrumPaletteMode = !spectrumPaletteMode;
            if (spectrumPaletteMode) {
                paletteMi.setText("Grey-scale Palette");
                backgroundMi.setText("Black Background");
                backgroundColor = 255;
                backgroundMi.setEnabled(true);
                cmi.setEnabled(true);
            } else {
                paletteMi.setText("Spectrum Palette");
                isHillshaded = false;
                cmi.setState(false);
                cmi.setEnabled(false);
                backgroundMi.setText("White Background");
                backgroundColor = 0;
                backgroundMi.setEnabled(false);
            }
            refresh();
        }
    }

    @Override
    public void mouseClicked(MouseEvent me) {
        if (me.getButton() == 3 || me.isPopupTrigger()) {
            myPopup.show((Component)me.getSource(), me.getX(), me.getY());
        }
    }

    @Override
    public void mousePressed(MouseEvent me) {

    }

    @Override
    public void mouseReleased(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void mouseEntered(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void mouseExited(MouseEvent me) {
        //throw new UnsupportedOperationException("Not supported yet.");
    }
    
}