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
import java.text.DecimalFormat;
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList
import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import whitebox.stats.TwoSampleKSTest
import groovy.transform.CompileStatic

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "CompareImagesForDifferences"
def descriptiveName = "Compare Images For Differences"
def description = "Test for signficant differences between two rasters."
def toolboxes = ["StatisticalTools"]

public class CompareImagesForDifferences implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd;
	private String descriptiveName
	
	public CompareImagesForDifferences(WhiteboxPluginHost pluginHost, 
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
			def helpFile = "CompareImagesForDifferences"
			sd.setHelpFile(helpFile)
		
			// Specifying the source file allows the 'view code' 
			// button on the tool dialog to be displayed.
			def pathSep = File.separator
			def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + "CompareImagesForDifferences.groovy"
			sd.setSourceFile(scriptFile)
			
			// add some components to the dialog
			sd.addDialogFile("Input image 1 file", "Input Image 1:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Input image 2 file", "Input Image 2:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogComboBox("Select a test type", "Test Type:", ["paired-sample t-test", "two-sample K-S test"], 0)
			sd.addDialogDataInput("Enter the desired sample size", "Sample Size (blank for whole image):", "", true, true)
            
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
	  	DecimalFormat df = new DecimalFormat("0.000");
        DecimalFormat df2 = new DecimalFormat("###,###,###,###");
        int progress, oldProgress
	  	int sampleSize = 0
		boolean useSampleBool = false;
		int row, col;
        double z1, z2;
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

		if (args.length != 4) {
			pluginHost.showFeedback("Incorrect number of arguments given to tool.")
			return
		}
		// read the input parameters
		String inputFile1 = args[0]
		String inputFile2 = args[1]
		String testType = "t-test"
		if (args[2].toLowerCase().contains("k-s") || 
		    args[2].toLowerCase().contains("ks")) {
			testType = "k-s"
		}
		if (args[3].toLowerCase().equals("not specified")) {
            useSampleBool = false;
            sampleSize = 0;
        } else {
            useSampleBool = true;
            sampleSize = Integer.parseInt(args[3]);
        }
		
		// check to see that the inputHeader1 and outputHeader are not null.
        if (inputFile1 == null || inputFile2 == null) {
            pluginHost.showFeedback("One or more of the input parameters have not been set properly.");
            return;
        }

        WhiteboxRaster image1 = new WhiteboxRaster(inputFile1, "r");
        int rows = image1.getNumberRows();
        int cols = image1.getNumberColumns();
        double noData1 = image1.getNoDataValue();

        WhiteboxRaster image2 = new WhiteboxRaster(inputFile2, "r");
        if (rows != image2.getNumberRows() || cols != image2.getNumberColumns()) {
            pluginHost.showFeedback("The input images must have the same dimensions (rows and columns).");
            return;
        }
        double noData2 = image2.getNoDataValue();

        
        // calculate the mean difference and the standard deviation of differences
        if (!useSampleBool) { //performing the test on the whole image.
            if (testType.contains("t-test")) {
	            double[] data1, data2;
	            oldProgress = -1
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
	                progress = (int)(100f * row / (rows - 1))
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
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

                pluginHost.returnData(retstr);
            } else if (testType.contains("k-s")) {
            	double[] data1, data2;
            	
                N = 0;
                oldProgress = -1
	            for (row = 0; row < rows; row++) {
	                data1 = image1.getRowValues(row);
	                data2 = image2.getRowValues(row);
	                for (col = 0; col < cols; col++) {
	                    z1 = data1[col];
	                    z2 = data2[col];
	                    if (z1 != noData1 && z2 != noData2) {
	                    	N++;
	                    }
	                }
	                progress = (int)(100f * row / (rows - 1))
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
	            }
				
	            double[] data1Dbl = new double[(int)N];
                double[] data2Dbl = new double[(int)N];
                oldProgress = -1
                int j = 0;
	            for (row = 0; row < rows; row++) {
	                data1 = image1.getRowValues(row);
	                data2 = image2.getRowValues(row);
	                for (col = 0; col < cols; col++) {
	                    z1 = data1[col];
	                    z2 = data2[col];
	                    if (z1 != noData1 && z2 != noData2) {
	                    	data1Dbl[j] = z1
	                    	data2Dbl[j] = z2
	                    	j++
	                    }
	                }
	                progress = (int)(100f * row / (rows - 1))
					if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
	            }
	            
                TwoSampleKSTest ks = new TwoSampleKSTest(data1Dbl, data2Dbl);
                double Dmax = ks.getDmax();
                double pValue = ks.getPvalue();
                
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
                    result = "The test REJECTS the null hypothesis that there is no significant difference between the distributions of the two images or sample pixel values \ndrawn from the two images.\n\n";
                } else {
                    result = "The test FAILS TO REJECT the null hypothesis that there is no significant difference between the distributions of the two images or sample pixel values \ndrawn from the two images.\n\n";
                }
                String caveat = "Caveat: Given a sufficiently large sample, extremely small and non-notable differences can be found to be statistically significant, \nand statistical significance says nothing about the practical significance of a difference.\n";
                
                retstr += result + caveat;
                
                pluginHost.returnData(retstr);
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
                    if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
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

                pluginHost.returnData(retstr);
            } else if (testType.contains("k-s")) {
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
                    if (progress > oldProgress) {
						pluginHost.updateProgress(progress)
						oldProgress = progress
					}
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
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
                    result = "The test REJECTS the null hypothesis that there is no significant difference between the distributions of the two images or sample pixel values \ndrawn from the two images.\n\n";
                } else {
                    result = "The test FAILS TO REJECT the null hypothesis that there is no significant difference between the distributions of the two images or sample pixel values \ndrawn from the two images.\n\n";
                }
                String caveat = "Caveat: Given a sufficiently large sample, extremely small and non-notable differences can be found to be statistically significant, \nand statistical significance says nothing about the practical significance of a difference.\n";
                
                retstr += result + caveat;
                
                pluginHost.returnData(retstr);
            }
        }
        image1.close();
        image2.close();

		// reset the progress bar
		pluginHost.updateProgress(0)
	  } catch (OutOfMemoryError oe) {
            pluginHost.showFeedback("An out-of-memory error has occurred during operation.")
      } catch (Exception e) {
            pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
            pluginHost.logException("Error in " + descriptiveName, e)
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

if (args == null) {
	pluginHost.showFeedback("Plugin arguments not set.")
} else {
	def p = new CompareImagesForDifferences(pluginHost, args, descriptiveName)
}
