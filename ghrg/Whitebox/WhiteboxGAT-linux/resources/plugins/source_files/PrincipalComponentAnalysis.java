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

import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.geospatialfiles.WhiteboxRasterInfo;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.interfaces.WhiteboxPlugin;
import java.util.Date;
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import java.util.SortedSet;
import java.util.TreeSet;
import java.text.DecimalFormat;
import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.geom.AffineTransform;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.io.*;
import java.awt.print.*;
import javax.print.attribute.*;
import whitebox.structures.ExtensionFileFilter;
import java.util.ArrayList;

/**
 * Principal component analysis (PCA) is a common data reduction technique that is used to reduce the dimensionality of multi-dimensional space. 
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class PrincipalComponentAnalysis implements WhiteboxPlugin {

    private WhiteboxPluginHost myHost = null;
    private String[] args;
    /**
     * Used to retrieve the plugin tool's name. This is a short, unique name containing no spaces.
     * @return String containing plugin name.
     */
    @Override
    public String getName() {
        return "PrincipalComponentAnalysis";
    }
    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer name (containing spaces) and is used in the interface to list the tool.
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Principal Component Analysis";
    }
    /**
     * Used to retrieve a short description of what the plugin tool does.
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Performs a principal component analysis (PCA) on a multi-spectral dataset.";
    }
    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"ImageTransformations", "ChangeDetection", "StatisticalTools"};
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
        if (myHost != null && ((progress != previousProgress)
                || (!progressLabel.equals(previousProgressLabel)))) {
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

        String inputFilesString = null;
        String[] imageFiles = null;
        String outputName = null;
        String workingDirectory = null;
        WhiteboxRasterInfo[] images = null;
        int nCols = 0;
        int nRows = 0;
        double z;
        int numImages;
        int progress = 0;
        int col, row;
        int a, i, j;
        double[] imageAverages;
        double[] imageTotals;
        double[] imageNumPixels;
        double[][] data;
        double[] noDataValues;
        String pathSep = File.separator;
        boolean standardizedPCA = false;
        int numberOfComponentImagesToCreate = 0;
        
        if (args.length <= 0) {
            showFeedback("Plugin parameters have not been set.");
            return;
        }

        // read the input parameters
        inputFilesString = args[0];
        outputName = args[1];
        if (outputName.toLowerCase().contains(".dep")) {
            outputName = outputName.replace(".dep", "");
        }
        standardizedPCA = Boolean.parseBoolean(args[2]);
        if (args[3].toLowerCase().contains("not")) { // not specified
            numberOfComponentImagesToCreate = 9999999;
        } else {
            numberOfComponentImagesToCreate = Integer.parseInt(args[3]);
        }
        
        try {
            // deal with the input images
            imageFiles = inputFilesString.split(";");
            numImages = imageFiles.length;
            images = new WhiteboxRasterInfo[numImages];
            
            imageAverages = new double[numImages];
            imageTotals = new double[numImages];
            imageNumPixels = new double[numImages];
            noDataValues = new double[numImages];
            data = new double[numImages][];
            
            for (i = 0; i < numImages; i++) {
                images[i] = new WhiteboxRasterInfo(imageFiles[i]);
                noDataValues[i] = images[i].getNoDataValue();
                if (i == 0) {
                    nCols = images[i].getNumberColumns();
                    nRows = images[i].getNumberRows();
                    File file = new File(imageFiles[i]);
                    workingDirectory = file.getParent();
                } else {
                    if (images[i].getNumberColumns() != nCols
                            || images[i].getNumberRows() != nRows) {
                        showFeedback("All input images must have the same dimensions (rows and columns).");
                        return;
                    }
                }
            }

            // Calculate the means
            for (row = 0; row < nRows; row++) {
                for (i = 0; i < numImages; i++) {
                    data[i] = images[i].getRowValues(row);
                }
                for (col = 0; col < nCols; col++) {
                    for (i = 0; i < numImages; i++) {
                        if (data[i][col] != noDataValues[i]) {
                            imageTotals[i] += data[i][col];
                            imageNumPixels[i]++;
                        }

                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (nRows - 1));
                updateProgress("Calculating image means:", progress);
            }
            
            for (i = 0; i < numImages; i++) {
                imageAverages[i] = imageTotals[i] / imageNumPixels[i];
            }
            
            // Calculate the covariance matrix and total deviations
            double[] imageTotalDeviation = new double[numImages];
            double[][] covariances = new double[numImages][numImages];
            double[][] correlationMatrix = new double[numImages][numImages];
            for (row = 0; row < nRows; row++) {
                for (i = 0; i < numImages; i++) {
                    data[i] = images[i].getRowValues(row);
                }
                for (col = 0; col < nCols; col++) {
                    for (i = 0; i < numImages; i++) {
                        if (data[i][col] != noDataValues[i]) {
                            imageTotalDeviation[i] += (data[i][col] - imageAverages[i])
                                    * (data[i][col] - imageAverages[i]);
                            for (a = 0; a < numImages; a++) {
                                if (data[a][col] != noDataValues[a]) {
                                    covariances[i][a] += (data[i][col] - imageAverages[i]) 
                                            * (data[a][col] - imageAverages[a]);
                                }
                            }
                        }

                    }
                }
                if (cancelOp) {
                    cancelOperation();
                    return;
                }
                progress = (int) (100f * row / (nRows - 1));
                updateProgress("Calculating covariances:", progress);
            }
            
            for (i = 0; i < numImages; i++) {
                for (a = 0; a < numImages; a++) {
                    correlationMatrix[i][a] = covariances[i][a] / (Math.sqrt(imageTotalDeviation[i] * imageTotalDeviation[a]));
                }
            }
            
            for (i = 0; i < numImages; i++) {
                for (a = 0; a < numImages; a++) {
                    covariances[i][a] = covariances[i][a] / (imageNumPixels[i] - 1);
                }
            }
            // Calculate the eigenvalues and eigenvectors
            Matrix cov = null;
            if (!standardizedPCA) {
                cov = new Matrix(covariances);
            } else {
                cov = new Matrix(correlationMatrix);
            }
            EigenvalueDecomposition eigen = cov.eig();
            double[] eigenvalues;
            Matrix eigenvectors;
            SortedSet<PrincipalComponent> principalComponents;
            eigenvalues = eigen.getRealEigenvalues();
            eigenvectors = eigen.getV();
            
            double[][] vecs = eigenvectors.getArray();
            int numComponents = eigenvectors.getColumnDimension(); // same as num rows.
            principalComponents = new TreeSet<PrincipalComponent>();
            for (i = 0; i < numComponents; i++) {
                double[] eigenvector = new double[numComponents];
                for (j = 0; j < numComponents; j++) {
                    eigenvector[j] = vecs[j][i];
                }
                principalComponents.add(new PrincipalComponent(eigenvalues[i], eigenvector));
            }
            
            double totalEigenvalue = 0;
            for (i = 0; i < numComponents; i++) {
                totalEigenvalue += eigenvalues[i];
            }
            
            double[][] explainedVarianceArray = new double[numComponents][2]; // percent and cum. percent
            j = 0;
            for (PrincipalComponent pc: principalComponents) {
                explainedVarianceArray[j][0] = pc.eigenValue / totalEigenvalue * 100.0;
                if (j == 0) {
                    explainedVarianceArray[j][1] = explainedVarianceArray[j][0];
                } else {
                    explainedVarianceArray[j][1] = explainedVarianceArray[j][0] + explainedVarianceArray[j - 1][1];
                }
                j++;
            }
            
            DecimalFormat df1 = new DecimalFormat("0.00");
            DecimalFormat df2 = new DecimalFormat("0.0000");
            DecimalFormat df3 = new DecimalFormat("0.000000");
            String ret = "Principal Component Analysis Report:\n\n";
            ret += "Component\tExplained Var.\tCum. %\tEigenvalue\tEigenvector\n";
            j = 0;
            for (PrincipalComponent pc: principalComponents) {
                
                String explainedVariance = df1.format(explainedVarianceArray[j][0]);
                String explainedCumVariance = df1.format(explainedVarianceArray[j][1]);
                double[] eigenvector = pc.eigenVector.clone();
                ret += (j + 1) + "\t" + explainedVariance + "\t" + explainedCumVariance + "\t" + df2.format(pc.eigenValue) + "\t";
                String eigenvec = "[";
                for (i = 0; i < numComponents; i++) {
                    if (i < numComponents - 1) {
                        eigenvec += df3.format(eigenvector[i]) + ", ";
                    } else {
                        eigenvec += df3.format(eigenvector[i]);
                    } 
                }
                eigenvec += "]";
                ret += eigenvec + "\n";
                
                if (j < numberOfComponentImagesToCreate) {
                    // now set up the output image
                    String outputHeader = workingDirectory + pathSep + outputName + "_comp" + (j + 1) + ".dep";
                    WhiteboxRaster output = new WhiteboxRaster(outputHeader, "rw",
                            imageFiles[0], WhiteboxRaster.DataType.FLOAT, 0);
                    output.setDataScale(DataScale.CONTINUOUS);

                    for (row = 0; row < nRows; row++) {
                        for (i = 0; i < numImages; i++) {
                            data[i] = images[i].getRowValues(row);
                        }
                        for (col = 0; col < nCols; col++) {
                            if (data[0][col] != noDataValues[0]) {
                                z = 0;
                                for (i = 0; i < numImages; i++) {
                                    z += data[i][col] * eigenvector[i];
                                }
                                output.setValue(row, col, z);
                            } else {
                                output.setValue(row, col, noDataValues[0]);
                            }
                        }
                        if (cancelOp) {
                            cancelOperation();
                            return;
                        }
                        progress = (int) (100f * row / (nRows - 1));
                        updateProgress("Creating component images:", progress);
                    }

                    output.addMetadataEntry("Created by the "
                            + getDescriptiveName() + " tool.");
                    output.addMetadataEntry("Created on " + new Date());
                    output.addMetadataEntry("Principal Component Num.: " + (j + 1));
                    output.addMetadataEntry("Eigenvalue: " + pc.eigenValue);
                    eigenvec = "[";
                    for (i = 0; i < numComponents; i++) {
                        if (i < numComponents - 1) {
                            eigenvec += eigenvector[i] + ", ";
                        } else {
                            eigenvec += eigenvector[i];
                        }
                    }
                    eigenvec += "]";
                    output.addMetadataEntry("Eigenvector: " + eigenvec);
                    if (!standardizedPCA) {
                        output.addMetadataEntry("PCA Type: unstandardized");
                    } else {
                        output.addMetadataEntry("PCA Type: standardized");
                    }
                    output.close();
                }
                
                j++;
            }
            
            // calculate the factor loadings.
            ret += "\nFactor Loadings:\n";
            ret += "\t\tComponent\n\t";
            for (i = 0; i < numComponents; i++) {
                ret += (i + 1) + "\t";
            }
            ret += "\n";
            double loading = 0;
            if (!standardizedPCA) {
                for (i = 0; i < numImages; i++) {
                    ret += "band" + (i + 1) + "\t";
                    for (PrincipalComponent pc : principalComponents) {
                        double[] eigenvector = pc.eigenVector.clone();
                        double ev = pc.eigenValue;
                        loading = (eigenvector[i] * Math.sqrt(ev)) / Math.sqrt(covariances[i][i]);
                        ret += df1.format(loading) + "\t";
                    }
                    ret += "\n";
                }
            } else {
                for (i = 0; i < numImages; i++) {
                    ret += "band" + (i + 1) + "\t";
                    for (PrincipalComponent pc : principalComponents) {
                        double[] eigenvector = pc.eigenVector.clone();
                        double ev = pc.eigenValue;
                        loading = (eigenvector[i] * Math.sqrt(ev));
                        ret += df1.format(loading) + "\t";
                    }
                    ret += "\n";
                }
            }
    
            
            for (i = 0; i < numImages; i++) {
                images[i].close();
            }
            
            returnData(ret);
            //System.out.println(ret);
            
            ScreePlot plot = new ScreePlot(explainedVarianceArray);
            returnData(plot);
            
            if (numComponents > 3) {
                for (i = 2; i >= 0; i--){ 
                    if (i < numberOfComponentImagesToCreate) {
                        // returning a header file string displays the image.
                        String outputHeader = workingDirectory + pathSep + outputName + "_comp" + (i + 1) + ".dep";
                        returnData(outputHeader);
                    }
                }
            } else {
                for (i = numComponents - 1; i >= 0; i--){ 
                    if (i < numberOfComponentImagesToCreate) {
                        // returning a header file string displays the image.
                        String outputHeader = workingDirectory + pathSep + outputName + "_comp" + (i + 1) + ".dep";
                        returnData(outputHeader);
                    }
                }
            }

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
    
    public static class PrincipalComponent implements Comparable<PrincipalComponent> {

        public double eigenValue;
        public double[] eigenVector;

        public PrincipalComponent(double eigenValue, double[] eigenVector) {
            this.eigenValue = eigenValue;
            this.eigenVector = eigenVector;
        }

        @Override
        public int compareTo(PrincipalComponent o) {
            int ret = 0;
            if (eigenValue > o.eigenValue) {
                ret = -1;
            } else if (eigenValue < o.eigenValue) {
                ret = 1;
            }
            return ret;
        }

        @Override
        public String toString() {
            String ret = "Principle Component, eigenvalue: " + eigenValue + ", eigenvector: [";
            for (int i = 0; i < eigenVector.length; i++) {
                ret += eigenVector[i] + ", ";
            }
            ret += "]";
            return ret;
        }
    }
    
//    // this is only used for debugging the tool
//    public static void main(String[] args) {
//        PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
//        args = new String[3];
//        //args[0] = "/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band1 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band2 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band3 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band4 clipped.dep;/Users/johnlindsay/Documents/Teaching/GEOG3420/Winter 2012/Labs/Lab1/Data/LE70180302002142EDC00/band5 clipped.dep;";
//        args[0] = "/Users/johnlindsay/Documents/Data/LandsatData/band1.dep;/Users/johnlindsay/Documents/Data/LandsatData/band2_cropped.dep;/Users/johnlindsay/Documents/Data/LandsatData/band3_cropped.dep;/Users/johnlindsay/Documents/Data/LandsatData/band4_cropped.dep";
//        args[1] = "PCA";
//        args[2]= "false";
//        
//        pca.setArgs(args);
//        pca.run();
//        
//    }
}

class ScreePlot extends JPanel implements ActionListener, Printable, MouseMotionListener, MouseListener {
    private int bottomMargin = 60;
    private int topMargin = 45;
    private int leftMargin = 80;
    private int rightMargin = 20;
    private JPopupMenu myPopup = null;
    private JMenuItem gridMi = null;
    private boolean gridOn = true;
    private int numComponents = 0;
    
    private double[][] plotData = null;
    
    // Constructors
    public ScreePlot(double[][] plotData) {
        this.plotData = plotData.clone();
        this.numComponents = plotData.length;
        this.setPreferredSize(new Dimension(350, 350));
        setMouseMotionListener();
        setMouseListener();
        setUp();
    }
    
    // Methods
    private void setUp() {
        try {
            createPopupMenus();
        } catch (Exception e) {
            
        }
    }
    
    private void setMouseMotionListener() {
        this.addMouseMotionListener(this);
    }
    
    private void setMouseListener() {
        this.addMouseListener(this);
    }
    
    public void refresh() {
        repaint();
    }
    
    @Override
    public void paint (Graphics g) {
        drawPlot(g);
    }
    
    int plotIndex = 0;
    private void drawPlot(Graphics g) {
        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g2d.setColor(Color.white);
        g2d.fillRect(0, 0, getWidth(), getHeight());
        
        double activeWidth = getWidth() - leftMargin - rightMargin;
        double activeHeight = getHeight() - topMargin - bottomMargin;
        int bottomY = getHeight() - bottomMargin;
        int rightX = getWidth() - rightMargin;
        
        // draw data line
        int x1, x2, y1, y2;
        g2d.setColor(Color.red);
        for (int i = 1; i < numComponents; i++) {
            x1 = (int)(leftMargin + ((double)(i - 1) / (numComponents - 1)) * activeWidth);
            y1 = (int)(bottomY - (plotData[i - 1][plotIndex] * 10.0) / 1000 * activeHeight);
            x2 = (int)(leftMargin + ((double)(i) / (numComponents - 1)) * activeWidth);
            y2 = (int)(bottomY - (plotData[i][plotIndex] * 10) / 1000 * activeHeight);
            g2d.drawLine(x1, y1, x2, y2);
        }
        
        // draw data points
        int radius = 2;
        for (int i = 0; i < numComponents; i++) {
            x1 = (int)(leftMargin + ((double)(i) / (numComponents - 1)) * activeWidth);
            y1 = (int)(bottomY - (plotData[i][plotIndex] * 10.0) / 1000 * activeHeight);
            g2d.drawOval(x1 - radius - 1, y1 - radius - 1, 2 * radius + 2, 2 * radius + 2);
        }
        
        // draw axes
        g2d.setColor(Color.black);
        g2d.drawLine(leftMargin, bottomY, rightX, bottomY);
        g2d.drawLine(leftMargin, bottomY, leftMargin, topMargin);
        g2d.drawLine(leftMargin, topMargin, rightX, topMargin);
        g2d.drawLine(rightX, bottomY, rightX, topMargin);
        
        // draw ticks
        int tickSize = 4;
        for (int i = 1; i <= numComponents; i++) {
            x1 = (int)(leftMargin + ((double)(i - 1) / (numComponents - 1)) * activeWidth);
            g2d.drawLine(x1, bottomY, x1, bottomY + tickSize);
        }
        
        for (int i = 0; i <= 1000; i += 100) {
            y1 = (int)(bottomY - i / 1000.0 * activeHeight);
            g2d.drawLine(leftMargin, y1, leftMargin - tickSize, y1);
        }
        
        // labels
        DecimalFormat df = new DecimalFormat("#,###,###.###");
        Font font = new Font("SanSerif", Font.PLAIN, 11);
        FontMetrics metrics = g.getFontMetrics(font);
        int hgt, adv;
        hgt = metrics.getHeight();
        
        // x-axis labels
        String label;
        for (int i = 1; i <= numComponents; i++) {
            label = String.valueOf(i);
            x1 = (int)(leftMargin + ((double)(i - 1) / (numComponents - 1)) * activeWidth);
            adv = metrics.stringWidth(label) / 2;
            g2d.drawString(label, x1 - adv, bottomY + hgt + 4);
        }
        label = "Component";
        adv = metrics.stringWidth(label);
        int xAxisMidPoint = (int)(leftMargin + activeWidth / 2);
        g2d.drawString(label, xAxisMidPoint - adv / 2, bottomY + 2 * hgt + 6);
        
        // y-axis labels
        
        // rotate the font
        Font oldFont = g.getFont();
        Font f = oldFont.deriveFont(AffineTransform.getRotateInstance(-Math.PI / 2.0));
        g2d.setFont(f);
        
        int yAxisMidPoint = (int)(topMargin + activeHeight / 2);
        int offset;
        label = "Explained Variance (%)";
        offset = metrics.stringWidth("100.0") + 12 + hgt;
        adv = metrics.stringWidth(label);
        g2d.drawString(label, leftMargin - offset, yAxisMidPoint + adv / 2);
        
        // replace the rotated font.
        g2d.setFont(oldFont);

        df = new DecimalFormat("0.0");
        for (int i = 0; i <= 1000; i += 100) {
            label = df.format(i / 10);
            y1 = (int)(bottomY - i / 1000.0 * activeHeight);
            adv = metrics.stringWidth(label);
            g2d.drawString(label, leftMargin - adv - 12, y1 + hgt / 2);
        }
        
        // title
        
        // bold font
        oldFont = g.getFont();
        font = font = new Font("SanSerif", Font.BOLD, 12);
        g2d.setFont(font);
        
        label = "PCA Scree Plot";
        adv = metrics.stringWidth(label);
        g2d.drawString(label, getWidth() / 2 - adv / 2, topMargin - hgt - 5);
        
        g2d.setFont(oldFont);
        
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
    
//    
    @Override
    public void mouseMoved(MouseEvent e) {
        
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
            
//            myPopup.addSeparator();
//            
//            gridMi = new JMenuItem("Turn Off Grid");
//            gridMi.addActionListener(this);
//            gridMi.setActionCommand("grid");
//            myPopup.add(gridMi);
//             
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
        String imageName = "ScreePlot.png";

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
        } else if (actionCommand.equals("grid")) {
            if (gridOn) {
                gridOn = false;
                gridMi.setText("Turn On Grid");
            } else if (!gridOn) {
                gridOn = true;
                gridMi.setText("Turn Off Grid");
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