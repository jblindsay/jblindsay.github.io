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

 /*
  * This tool performs a conditional evaluation (if-then-else) on 
  * a raster image.
  */

import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.Future
import java.util.concurrent.*
import java.util.Date
import java.util.ArrayList

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineFactory;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import whitebox.interfaces.WhiteboxPluginHost
import whitebox.geospatialfiles.WhiteboxRaster
import whitebox.geospatialfiles.WhiteboxRasterInfo
import whitebox.geospatialfiles.WhiteboxRasterBase.DataType
import whitebox.geospatialfiles.WhiteboxRasterBase.DataScale
import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.*
import whitebox.ui.plugin_dialog.ScriptDialog
import whitebox.utilities.StringUtilities
import groovy.transform.CompileStatic
import groovy.time.TimeDuration
import groovy.time.TimeCategory

// The following four variables are required for this 
// script to be integrated into the tool tree panel. 
// Comment them out if you want to remove the script.
def name = "ConditionalEvaluation"
def descriptiveName = "Conditional Evaluation"
def description = "Performs a conditional evaluaton (if-then-else) operation on a raster."
def toolboxes = ["MathTools", "GISTools"]

public class ConditionalEvaluation implements ActionListener {
	private WhiteboxPluginHost pluginHost
	private ScriptDialog sd
	private String descriptiveName
	private String name
	
	public ConditionalEvaluation(WhiteboxPluginHost pluginHost, 
		String[] args, def name, def descriptiveName) {
		this.pluginHost = pluginHost
		this.name = name
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
			sd.addDialogFile("Input raster file", "Input Raster:", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogDataInput("<html>Conditional statement.\nThe <i>conditional statement</i> can be any valid Groovy statement, e.g. value > 35.0.</html>", "Conditional Statement e.g. value > 35.0:", "", false, false)
            sd.addDialogFile("Value  where TRUE (input raster or constant value)", "Value Where TRUE (Raster File Or Constant Value):", "open", "Raster Files (*.dep), DEP", true, false)
			sd.addDialogFile("Value  where FALSE (input raster or constant value). If this parameter is not specified, NoData will be used.", "Value Where FALSE (Raster File Or Constant Value):", "open", "Raster Files (*.dep), DEP", true, true)
			sd.addDialogFile("Output raster file", "Output Raster File:", "save", "Raster Files (*.dep), DEP", true, false)
            
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
			//Date start = new Date()
			double value, trueValue = 0, falseValue = 0
			boolean ret
			int progress, oldProgress
			String trueValueFile = ""
			String falseValueFile = ""
			boolean trueValueConstant = true
			boolean falseValueConstant = true
			
			if (args.length < 4) {
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return
			}
			
			// read the input parameters

			// the input file
			String inputFile = args[0]
			WhiteboxRaster image = new WhiteboxRaster(inputFile, "r")
			double nodata = image.getNoDataValue()
			int rows = image.getNumberRows()
			int cols = image.getNumberColumns()

			// the conditional statement
			String conStatement = args[1]
			boolean hasConStatement = !(conStatement.isEmpty() | conStatement.toLowerCase().equals("not specified"))

			// the TRUE value
			byte trueType = -1
			String trueExpression = ""
			if (args[2] != null && !args[2].isEmpty() && !(args[2].toLowerCase().equals("not specified"))) {
				// see if it's a file
				def file = new File(args[2]);
				if ((file).exists()) {
					trueValueFile = args[2];
					trueType = 1; // raster file
				} else {
					String str = file.getName().replace(".dep", "")
					if (StringUtilities.isNumeric(str)) {
						trueValue = Double.parseDouble(str);
						trueType = 0; // constant numerical value
					} else {
						trueExpression = str;
						trueType = 2; // expression
					}
				}
			} else {
				trueValue = nodata;
				trueType = 0; // constant numerical value
			}

			// the FALSE value
			byte falseType = -1
			String falseExpression = ""
			if (args[3] != null && !args[3].isEmpty() && !(args[3].toLowerCase().equals("not specified"))) {
				// see if it's a file
				def file = new File(args[3]);
				if ((file).exists()) {
					falseValueFile = args[3];
					falseType = 1; // raster file
				} else {
					String str = file.getName().replace(".dep", "")
					if (StringUtilities.isNumeric(str)) {
						falseValue = Double.parseDouble(str);
						falseType = 0; // constant numerical value
					} else {
						falseExpression = str;
						falseType = 2; // expression
					}
				}
			} else {
				falseValue = nodata;
				falseType = 0; // constant numerical value
			}
			
			// the output file
			String outputFile = args[4]
			
			WhiteboxRaster trueValueRaster
			if (trueType == 1) { //!trueValueConstant) {
				trueValueRaster = new WhiteboxRaster(trueValueFile, "r")
				if (trueValueRaster.getNumberRows() != rows || trueValueRaster.getNumberColumns() != cols) {
					pluginHost.showFeedback("The input images must all have the same dimensions \n(i.e. number of rows and columns and spatial extent)")
					return
				}
			}

			WhiteboxRaster falseValueRaster
			if (falseType == 1) { //!falseValueConstant) {
				falseValueRaster = new WhiteboxRaster(falseValueFile, "r")
				if (falseValueRaster.getNumberRows() != rows || falseValueRaster.getNumberColumns() != cols) {
					pluginHost.showFeedback("The input images must all have the same dimensions \n(i.e. number of rows and columns and spatial extent)")
					return
				}
			}

			WhiteboxRaster output = new WhiteboxRaster(outputFile, "rw", 
  		  	     inputFile, DataType.FLOAT, nodata)
  		  	output.setNoDataValue(nodata)
			if (trueType == 0 && falseType == 0) {
				output.setPreferredPalette("qual.plt")
				output.setDataScale(DataScale.CATEGORICAL)
			} else if (trueType == 1) {
				output.setPreferredPalette(trueValueRaster.getPreferredPalette())
				output.setDataScale(trueValueRaster.getDataScale())
			} else if (falseType == 1) {
				output.setPreferredPalette(falseValueRaster.getPreferredPalette())
				output.setDataScale(falseValueRaster.getDataScale())
			} else {
				output.setPreferredPalette("qual.plt")
				//output.setDataScale(DataScale.CATEGORICAL)
			}

			ScriptEngineManager mgr = new ScriptEngineManager();
   			ScriptEngine engine = mgr.getEngineByName("groovy");
   			pluginHost.showFeedback(engine.NAME);
   			//ScriptEngine engine = mgr.getEngineByName("python");
   			//ScriptEngine engine = mgr.getEngineByName("javascript");
   			String expression = conStatement.replace("\"Value\"", "value").replace("\'Value\'", "value").replace("\"value\"", "value").replace("\'value\'", "value").replace("VALUE", "value").replace("NoData", "nodata").replace("NODATA", "nodata").replace("Null", "nodata").replace("NULL", "nodata").replace("null", "nodata").replace(";", "");
   			boolean expressionEvaluatesNoData = expression.contains("nodata")
   			String myScript = "for (int column = 0; column < $cols; column++) { double value = inData[column]; if (value != nodata || expressionEvaluatesNoData) { if ($expression) { outData[column] = trueData[column]; } else { outData[column] = falseData[column]; } } else { outData[column] = nodata; } }"
   			//String myScript = "for column in xrange(0, $cols):\n\tvalue = inData[column]\n\tif (value != nodata) | expressionEvaluatesNoData:\n\t\tif $expression:\n\t\t\toutData[column] = trueData[column]\n\t\telse:\n\t\t\toutData[column] = falseData[column]\n\telse:\n\t\toutData[column] = nodata\n"
   			//String myScript = "for (column = 0; column < $cols; column++) { value = inData[column]; if (value != nodata || expressionEvaluatesNoData) { if ($expression) { outData[column] = trueData[column]; } else { outData[column] = falseData[column]; } } else { outData[column] = nodata; } }"
   			
   			//pluginHost.showFeedback(myScript)
   			CompiledScript generate_data = ((Compilable) engine).compile(myScript);
   			double oldValue = nodata
   			boolean oldReturn
			Bindings bindings = engine.createBindings();

			double[] outData = new double[cols];
			double[] trueData = new double[cols];
			double[] falseData = new double[cols];
			if (trueType == 0) {
				Arrays.fill(trueData, trueValue)
				bindings.put("trueData", trueData);
			}
			boolean addRowNumberTrue = false;
			boolean addRowYTrue = false;
			boolean addValueTrue = false;
			if (trueType == 2) {
				if (trueExpression.toLowerCase().equals("nodata") ||
				     trueExpression.toLowerCase().equals("null")) {
					Arrays.fill(trueData, nodata);
					bindings.put("trueData", trueData);				
				} else if (trueExpression.toLowerCase().equals("rows")) {
					Arrays.fill(trueData, rows);
					bindings.put("trueData", trueData);			
				} else if (trueExpression.toLowerCase().equals("columns")) {
					Arrays.fill(trueData, cols);
					bindings.put("trueData", trueData);			
				} else if (trueExpression.toLowerCase().equals("minvalue")) {
					Arrays.fill(trueData, image.getMinimumValue());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("maxvalue")) {
					Arrays.fill(trueData, image.getMaximumValue());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("displayminvalue")) {
					Arrays.fill(trueData, image.getDisplayMinimum());
					bindings.put("trueData", trueData);			
				} else if (trueExpression.toLowerCase().equals("displaymaxvalue")) {
					Arrays.fill(trueData, image.getDisplayMaximum());
					bindings.put("trueData", trueData);			
				} else if (trueExpression.toLowerCase().equals("north")) {
					Arrays.fill(trueData, image.getNorth());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("south")) {
					Arrays.fill(trueData, image.getSouth());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("east")) {
					Arrays.fill(trueData, image.getEast());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("west")) {
					Arrays.fill(trueData, image.getWest());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("cellsizex")) {
					Arrays.fill(trueData, image.getCellSizeX());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("cellsizey")) {
					Arrays.fill(trueData, image.getCellSizeY());
					bindings.put("trueData", trueData);		
				} else if (trueExpression.toLowerCase().equals("row")) {
					addRowNumberTrue = true;
				} else if (trueExpression.toLowerCase().equals("column")) {
					for (int col = 0; col < cols; col++) {
						trueData[col] = col
					}
					bindings.put("trueData", trueData);
				} else if (trueExpression.toLowerCase().equals("rowy")) {
					addRowYTrue = true;
				} else if (trueExpression.toLowerCase().equals("columnx")) {
					for (int col = 0; col < cols; col++) {
						trueData[col] = image.getXCoordinateFromColumn(col);
					}
					bindings.put("trueData", trueData);
				} else if (trueExpression.toLowerCase().equals("value")) {
					addValueTrue = true
				} else {
					pluginHost.showFeedback("Invalid TRUE value.");
					return;
				}
			}
			if (falseType == 0) { //falseValueConstant) {
				Arrays.fill(falseData, falseValue)
				bindings.put("falseData", falseData);
			}
			boolean addRowNumberFalse = false;
			boolean addRowYFalse = false;
			boolean addValueFalse = false;
			if (falseType == 2) {
				if (falseExpression.toLowerCase().equals("nodata") ||
				     falseExpression.toLowerCase().equals("null")) {
					Arrays.fill(falseData, nodata);
					bindings.put("falseData", falseData);				
				} else if (falseExpression.toLowerCase().equals("rows")) {
					Arrays.fill(falseData, rows);
					bindings.put("falseData", falseData);			
				} else if (falseExpression.toLowerCase().equals("columns")) {
					Arrays.fill(falseData, cols);
					bindings.put("falseData", falseData);			
				} else if (falseExpression.toLowerCase().equals("minvalue")) {
					Arrays.fill(falseData, image.getMinimumValue());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("maxvalue")) {
					Arrays.fill(falseData, image.getMaximumValue());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("displayminvalue")) {
					Arrays.fill(falseData, image.getDisplayMinimum());
					bindings.put("falseData", falseData);			
				} else if (falseExpression.toLowerCase().equals("displaymaxvalue")) {
					Arrays.fill(falseData, image.getDisplayMaximum());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("north")) {
					Arrays.fill(falseData, image.getNorth());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("south")) {
					Arrays.fill(falseData, image.getSouth());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("east")) {
					Arrays.fill(falseData, image.getEast());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("west")) {
					Arrays.fill(falseData, image.getWest());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("cellsizex")) {
					Arrays.fill(falseData, image.getCellSizeX());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("cellsizey")) {
					Arrays.fill(falseData, image.getCellSizeY());
					bindings.put("falseData", falseData);		
				} else if (falseExpression.toLowerCase().equals("row")) {
					addRowNumberFalse = true;
				} else if (falseExpression.toLowerCase().equals("column")) {
					for (int col = 0; col < cols; col++) {
						falseData[col] = col
					}
					bindings.put("falseData", falseData);
				} else if (falseExpression.toLowerCase().equals("rowy")) {
					addRowYFalse = true;
				} else if (falseExpression.toLowerCase().equals("columnx")) {
					for (int col = 0; col < cols; col++) {
						falseData[col] = image.getXCoordinateFromColumn(col);
					}
					bindings.put("falseData", falseData);
				} else if (falseExpression.toLowerCase().equals("value")) {
					addValueFalse = true
				} else {
					pluginHost.showFeedback("Invalid FALSE value.");
					return;
				}
			}
			bindings.put("outData", outData);
			bindings.put("nodata", nodata);
			bindings.put("expressionEvaluatesNoData", expressionEvaluatesNoData);
			bindings.put("maxvalue", image.getMaximumValue());
			bindings.put("minvalue", image.getMinimumValue());
			bindings.put("displaymaxvalue", image.getDisplayMaximum());
			bindings.put("displayminvalue", image.getDisplayMinimum());
			bindings.put("rows", rows);
			bindings.put("columns", cols);
			bindings.put("north", image.getNorth());
			bindings.put("south", image.getSouth());
			bindings.put("east", image.getEast());
			bindings.put("west", image.getWest());
			bindings.put("cellsizex", image.getCellSizeX());
			bindings.put("cellsizey", image.getCellSizeY());
			bindings.put("raster", image);

			if (addValueTrue || addValueFalse) {
				output.setPreferredPalette(image.getPreferredPalette())
				output.setDataScale(image.getDataScale())
			}
			
			for (int row = 0; row < rows; row++) {
				double[] inData = image.getRowValues(row);
				if (trueType == 1) {
					trueData = trueValueRaster.getRowValues(row)
					bindings.put("trueData", trueData);
				} else if (addValueTrue) {
					trueData = image.getRowValues(row)
					bindings.put("trueData", trueData);
				}
				if (falseType == 1) {
					falseData = falseValueRaster.getRowValues(row)
					bindings.put("falseData", falseData);
				} else if (addValueFalse) {
					falseData = image.getRowValues(row)
					bindings.put("falseData", falseData);
				}
				if (addRowNumberTrue) {
					Arrays.fill(trueData, row);
					bindings.put("trueData", trueData);
				}
				if (addRowNumberFalse) {
					Arrays.fill(falseData, row);
					bindings.put("falseData", falseData);
				}
				if (addRowYTrue) {
					value = image.getYCoordinateFromRow(row)
					Arrays.fill(trueData, value);
					bindings.put("trueData", trueData);
				}
				if (addRowYFalse) {
					value = image.getYCoordinateFromRow(row)
					Arrays.fill(falseData, value);
					bindings.put("falseData", falseData);
				}
				bindings.put("inData", inData);
				bindings.put("row", row);
				generate_data.eval(bindings);
				
				
				output.setRowValues(row, outData)
  				progress = (int)(100f * row / rows)
				if (progress != oldProgress) {
					pluginHost.updateProgress(progress)
					oldProgress = progress
					// check to see if the user has requested a cancellation
					if (pluginHost.isRequestForOperationCancelSet()) {
						pluginHost.showFeedback("Operation cancelled")
						return
					}
				}
				
			}
			
			image.close()
			
			output.addMetadataEntry("Created by the "
	                    + descriptiveName + " tool.")
	        output.addMetadataEntry("Created on " + new Date())
			output.close()

//			Date stop = new Date()
//			TimeDuration td = TimeCategory.minus(stop, start)
//			pluginHost.showFeedback("Elapsed time: $td")
	
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

	class PixelValue {
		double value
		
		public PixelValue(double value) {
			this.value = value
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
	def f = new ConditionalEvaluation(pluginHost, args, name, descriptiveName)
}
