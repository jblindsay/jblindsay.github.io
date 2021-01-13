/* PLEASE NOTE THAT THIS TOOL HAS BEEN REPLACED WITH A GROOVY SCRIPT. */

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

import java.util.ArrayList;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
import whitebox.geospatialfiles.WhiteboxRaster;
import whitebox.interfaces.WhiteboxPlugin;
import whitebox.interfaces.WhiteboxPluginHost;
import whitebox.internationalization.WhiteboxInternationalizationTools;
import whitebox.stats.TwoSampleKSTest;

/**
 * This tool can be used to compare two images for statistically significant differences by performing either a paired-sample t-test or a two-sample Kolmogorov-Smirnov (K-S) test on the two images..
 *
 * @author Dr. John Lindsay email: jlindsay@uoguelph.ca
 */
public class CompareImagesForDifferences implements WhiteboxPlugin {

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
        return "CompareImagesForDifferences";
    }

    /**
     * Used to retrieve the plugin tool's descriptive name. This can be a longer
     * name (containing spaces) and is used in the interface to list the tool.
     *
     * @return String containing the plugin descriptive name.
     */
    @Override
    public String getDescriptiveName() {
        return "Compare Images For Signifcant Differences";
    }

    /**
     * Used to retrieve a short description of what the plugin tool does.
     *
     * @return String containing the plugin's description.
     */
    @Override
    public String getToolDescription() {
        return "Test for signficant differences between two rasters using a paired-sample t-test.";
    }

    /**
     * Used to identify which toolboxes this plugin tool should be listed in.
     *
     * @return Array of Strings.
     */
    @Override
    public String[] getToolbox() {
        String[] ret = {"StatisticalTools"};
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
        if (myHost != null && ((progress != previousProgress)
                || (!progressLabel.equals(previousProgressLabel)))) {
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

        String inputHeader1 = null;
        String inputHeader2 = null;
        boolean useSampleBool = false;
        int sampleSize = 0;
        DecimalFormat df = new DecimalFormat("0.000");
        DecimalFormat df2 = new DecimalFormat("###,###,###,###");
        
        try {

            if (args.length <= 0) {
                showFeedback("Plugin parameters have not been set.");
                return;
            }

            inputHeader1 = args[0];
            inputHeader2 = args[1];
            if (args[2].toLowerCase().equals("not specified")) {
                useSampleBool = false;
                sampleSize = 0;
            } else {
                useSampleBool = true;
                sampleSize = Integer.parseInt(args[2]);
            }
            String testType = "ks";

            // check to see that the inputHeader1 and outputHeader are not null.
            if (inputHeader1 == null || inputHeader2 == null) {
                showFeedback("One or more of the input parameters have not been set properly.");
                return;
            }

            int row, col;
            double z1, z2;
            int progress = 0;

            WhiteboxRaster image1 = new WhiteboxRaster(inputHeader1, "r");
            int rows = image1.getNumberRows();
            int cols = image1.getNumberColumns();
            double noData1 = image1.getNoDataValue();

            WhiteboxRaster image2 = new WhiteboxRaster(inputHeader2, "r");
            if (rows != image2.getNumberRows() || cols != image2.getNumberColumns()) {
                showFeedback("The input images must have the same dimensions (rows and columns).");
                return;
            }
            double noData2 = image2.getNoDataValue();

            double total, total1, total2, mean, mean1, mean2,
                    stdDev, stdDev1, stdDev2, t, stdErr;
            double totalSquared = 0;
            double totalSquared1 = 0;
            double totalSquared2 = 0;
            double variance;
            long N = 0;
            total = 0;
            total1 = 0;
            total2 = 0;

            // calculate the mean difference and the standard deviation of differences
            if (!useSampleBool) { //performing the test on the whole image.
                double[] data1, data2;
                for (row = 0; row < rows; row++) {
                    data1 = image1.getRowValues(row);
                    data2 = image2.getRowValues(row);
                    for (col = 0; col < cols; col++) {
                        z1 = data1[col];
                        z2 = data2[col];
                        if (z1 != noData1 && z2 != noData2) {
                            total1 += z1;
                            total2 += z2;
                            total += z1 - z2;
                            totalSquared += (z1 - z2) * (z1 - z2);
                            totalSquared1 += z1 * z1;
                            totalSquared2 += z2 * z2;
                            N++;
                        }
                    }
                    if (cancelOp) {
                        cancelOperation();
                        return;
                    }
                    progress = (int) (100f * row / (rows - 1));
                    updateProgress(progress);
                }

            } else {
                //double[] sample = new double[sampleSize];
                Random generator = new Random();
                int[][] rowsAndColumns = new int[sampleSize][2];
                int sampleNumber = 0;
                while (sampleNumber < sampleSize) {
                    rowsAndColumns[sampleNumber][0] = generator.nextInt(rows);
                    rowsAndColumns[sampleNumber][1] = generator.nextInt(cols);
                    sampleNumber++;
                }

                Arrays.sort(rowsAndColumns, new Comparator<int[]>() {

                    @Override
                    public int compare(final int[] entry1, final int[] entry2) {
                        final int int1 = entry1[0];
                        final int int2 = entry2[0];
                        return Integer.valueOf(int1).compareTo(int2);
                    }
                });

                if (testType.contains("t-test")) {

                    for (int i = 0; i < sampleSize; i++) {
                        row = rowsAndColumns[i][0];
                        col = rowsAndColumns[i][1];
                        z1 = image1.getValue(row, col);
                        z2 = image2.getValue(row, col);
                        if (z1 != noData1 && z2 != noData2) {
                            total1 += z1;
                            total2 += z2;
                            total += (z1 - z2);
                            totalSquared += (z1 - z2) * (z1 - z2);
                            totalSquared1 += z1 * z1;
                            totalSquared2 += z2 * z2;
                            N++;
                        }
                        progress = (int) (100f * sampleNumber / (sampleSize - 1));
                        updateProgress(progress);
                    }

                    mean = total / N;
                    mean1 = total1 / N;
                    mean2 = total2 / N;
                    variance = (N * totalSquared - total * total) / (N * (N - 1));
                    stdDev = Math.sqrt(variance);
                    stdDev1 = Math.sqrt((N * totalSquared1 - total1 * total1) / (N * (N - 1)));
                    stdDev2 = Math.sqrt((N * totalSquared2 - total2 * total2) / (N * (N - 1)));

                    stdErr = stdDev / Math.sqrt(N);
                    t = mean / stdErr;
                    double nu = N - 1;

                    double pValue = 1 - tTest(t, nu);

                    String retstr = "";
                    retstr = "Paired-Samples t-Test Results:\n\n";
                    retstr = retstr + "Input Image1:\t\t" + image1.getShortHeaderFile() + "\n";
                    retstr = retstr + "Image1 Mean:\t\t" + df.format(mean1) + "\n";
                    retstr = retstr + "Image1 SD:\t\t" + df.format(stdDev1) + "\n\n";
                    retstr = retstr + "Input Image2:\t\t" + image2.getShortHeaderFile() + "\n";
                    retstr = retstr + "Image2 Mean:\t\t" + df.format(mean2) + "\n";
                    retstr = retstr + "Image2 SD:\t\t" + df.format(stdDev2) + "\n\n";
                    retstr = retstr + "Sample Size (N):\t" + df2.format(N) + "\n";
                    retstr = retstr + "Test Statistic (t):\t" + df.format(t) + "\n";
                    if (pValue > 0.001) {
                        retstr = retstr + "Significance (p-value):\t" + df.format(pValue) + "\n\n";
                    } else {
                        retstr = retstr + "Significance (p-value):\t<0.001\n\n";
                    }
                    String result;
                    if (pValue < 0.05) {
                        result = "The test REJECTS the null hypothesis that there is no significant difference between the means of the two images or sample pixel values \ndrawn from the two images.\n\n";
                    } else {
                        result = "The test FAILS TO REJECT the null hypothesis that there is no significant difference between the means of the two images or sample pixel values \ndrawn from the two images.\n\n";
                    }
                    String caveat = "Caveat: Given a sufficiently large sample, extremely small and non-notable differences can be found to be statistically significant, \nand statistical significance says nothing about the practical significance of a difference.\n";

                    retstr += result + caveat;

                    returnData(retstr);
                } else if (testType.contains("ks")) {
                    ArrayList<Double> data1 = new ArrayList<>();
                    ArrayList<Double> data2 = new ArrayList<>();

                    N = 0;
                    for (int i = 0; i < sampleSize; i++) {
                        row = rowsAndColumns[i][0];
                        col = rowsAndColumns[i][1];
                        z1 = image1.getValue(row, col);
                        z2 = image2.getValue(row, col);
                        if (z1 != noData1 && z2 != noData2) {
                            data1.add(z1);
                            data2.add(z2);
                            N++;
                        }
                        progress = (int) (100f * sampleNumber / (sampleSize - 1));
                        updateProgress(progress);
                    }
                    double[] data1Dbl = new double[data1.size()];
                    double[] data2Dbl = new double[data2.size()];
                    for (int i = 0; i < data1.size(); i++) {
                        data1Dbl[i] = data1.get(i);
                        data2Dbl[i] = data2.get(i);
                    }
                    TwoSampleKSTest ks = new TwoSampleKSTest(data1Dbl, data2Dbl);
                    double Dmax = ks.getDmax();
                    double pValue = ks.getPvalue();
                    System.out.println("Dmax = " + Dmax);
                    System.out.println("p-value = " + pValue);
                    
                    String retstr = "K-S Test Results:\n\n";
                    retstr = retstr + "Input Image1:\t\t" + image1.getShortHeaderFile() + "\n";
                    retstr = retstr + "Input Image2:\t\t" + image2.getShortHeaderFile() + "\n";
                    retstr = retstr + "Test Statistic (Dmax):\t" + df.format(Dmax) + "\n";
                    retstr = retstr + "Sample Size (N):\t" + df2.format(N) + "\n";
                    if (pValue > 0.001) {
                        retstr = retstr + "Significance (p-value):\t" + df.format(pValue) + "\n\n";
                    } else {
                        retstr = retstr + "Significance (p-value):\t<0.001\n\n";
                    }
                    String result = "";
                    if (pValue < 0.05) {
                        result = "The test REJECTS the null hypothesis that there is no significant difference between the means of the two images or sample pixel values \ndrawn from the two images.\n\n";
                    } else {
                        result = "The test FAILS TO REJECT the null hypothesis that there is no significant difference between the means of the two images or sample pixel values \ndrawn from the two images.\n\n";
                    }
                    String caveat = "Caveat: Given a sufficiently large sample, extremely small and non-notable differences can be found to be statistically significant, \nand statistical significance says nothing about the practical significance of a difference.\n";
                    
                    retstr += result + caveat;
                    
                    returnData(retstr);
                }
            }
            image1.close();
            image2.close();

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

    // I got this code for calculating the p-value of a t-test from:
    // http://pfc-gtg-mspacman.googlecode.com/svn/trunk/pfc-gtg-mspacman/src/stats/StatisticalTests.java
    /**
     *
     * Applies the two-sided t-test given the value of t and nu. To do this it
     * calls betai.
     *
     * @param t t statistic
     * @param nu sample size minus one
     * @return p-value
     */
    public static double tTest(double t, double nu) {
        double a = nu / 2.0;
        double b = 0.5;
        double x = nu / (nu + t * t);
        return 1.0 - betai(a, b, x); 
    }

    protected static double betai(double a, double b, double x) {
        // can be used to find t statistic
        double bt;
        if ((x < 0.0) || (x > 1.0)) {
            System.out.println("Error in betai: " + x);
        }
        if ((x == 0.0) || (x == 1.0)) {
            bt = 0.0;
        } else {
            bt = Math.exp(gammln(a + b) - gammln(a) - gammln(b)
                    + a * Math.log(x)
                    + b * Math.log(1.0 - x));
        }
        if (x < (a + 1.0) / (a + b + 2.0)) {
            return bt * betacf(a, b, x) / a;
        } else {
            return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
        }
    }

    protected static double gammln(double xx) {
        double stp = 2.50662827465;

        double x, tmp, ser;
        x = xx - 1.0;
        tmp = x + 5.5;
        tmp = (x + 0.5) * Math.log(tmp) - tmp;
        ser = 1.0
                + 76.18009173 / (x + 1.0)
                - 86.50532033 / (x + 2.0)
                + 24.01409822 / (x + 3.0)
                - 1.231739516 / (x + 4.0)
                + 0.120858003 / (x + 5.0)
                - 0.536382e-5 / (x + 6.0);

        return tmp + Math.log(stp * ser); // finish
    }

    protected static double betacf(double a, double b, double x) {

        int maxIts = 100;
        double eps = 3.0e-7;

        double tem, qap, qam, qab, em, d;
        double bz, bpp, bp, bm, az, app;
        double am, aold, ap;

        am = 1.0;
        bm = 1.0;
        az = 1.0;
        qab = a + b;
        qap = a + 1.0;
        qam = a - 1.0;
        bz = 1.0 - qab * x / qap;
        for (int m = 1; m <= maxIts; m++) {
            em = m;
            tem = em + em;
            d = em * (b - m) * x / ((qam + tem) * (a + tem));
            ap = az + d * am;
            bp = bz + d * bm;
            d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
            app = ap + d * az;
            bpp = bp + d * bz;
            aold = az;
            am = ap / bpp;
            bm = bp / bpp;
            az = app / bpp;
            bz = 1.0;
            if (Math.abs(az - aold) < eps * Math.abs(az)) {
                return az;
            }
        }
        System.out.println("a or b too big, or maxIts too small");
        return -1;
    }
}
