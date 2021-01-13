''' 
Copyright (C) 2014 Dr. John Lindsay <jlindsay@uoguelph.ca>

This program is intended for instructional purposes only. The 
following is an example of how to use Whitebox's scripting 
capabilities to automate a geoprocessing workflow. The scripting 
language is Python; more specifically it is Jython, the Python 
implementation targeting the Java Virtual Machine (JVM).

In this script, we will take a digital elevation model (DEM), 
remove all the topographic depressions from it (i.e. hydrologically 
correct the DEM), calculate a flow direction pointer grid, use 
the pointer file to perform a flow accumulation (i.e. upslope area)
calculation, then threshold the upslope area to derive valley lines 
or streams. This is a fairly common workflow in spatial hydrology.

When you run a script from within Whitebox, a reference to the
Whitebox user interface (UI) will be automatically bound to your
script. It's variable name is 'pluginHost'. This is the primary 
reason why the script must be run from within Whitebox's Scripter.

First we need the directory containing the data, and to set 
the working directory to this. We will use the Vermont DEM contained
within the samples directory.
'''

import os

try:
	separator = os.sep # The system-specific directory separator
	wd = pluginHost.getApplicationDirectory() + separator + "resources" + separator + "samples" + separator + "Vermont DEM" + separator
	pluginHost.setWorkingDirectory(wd)
	
	demFile = wd + "Vermont DEM.dep" 
	# Notice that spaces are allowed in file names. There is also no 
	# restriction on the length of the file name...in fact longer,
	# descriptive names are preferred. Whitebox is friendly!
	
	# A raster or vector file can be displayed by specifying the file
	# name as an argument of the returnData method of the pluginHost
	pluginHost.returnData(demFile)
	
	''' 
	Remove the depressions in the DEM using the 'FillDepressions' tool.
	
	The help file for each tool in Whitebox contains a section detailing 
	the required input parameters needed to run the tool from a script.
	These parameters are always fed to the tool in a String array, in  
	the case below, called 'args'. The tool is then run using the 'runPlugin'
	method of the pluginHost. runPlugin takes the name of the tool (see 
	the tool's help for the proper name), the arguments string array, 
	followed by two Boolean arguments. The first of these Boolean 
	arguments determines whether the plugin will be run on its own
	separate thread. In most scripting applications, this should be set
	to 'False' because the results of this tool are needed as inputs to
	subsequent tools. The second Boolean argument specifies whether the
	data that are returned to the pluginHost after the tool is completed 
	should be suppressed. Many tools will automatically display images 
	or shapefiles or some text report when they've completed. It is often
	the case in a workflow that you only want the final result to be 
	displayed, in which case all of the runPlugins should have this final
	Boolean parameter set to 'True' except for the last operation, for 
	which it should be set to 'False' (i.e. don't suppress the output).
	The data will still be written to disc if the output are supressed, 
	they simply won't be automatically displayed when the tool has 
	completed. If you don't specify this last Boolean parameter, the 
	output will be treated as normal.
	'''
	filledDEMFile = wd + "filled DEM.dep"
	flatIncrement = "0.001" # Notice that although this is a numeric parameter, it is provided to the tool as a string.
	args = [demFile, filledDEMFile, flatIncrement]
	pluginHost.runPlugin("FillDepressions", args, False, True)
	
	# Calculate the D8 pointer (flow direction) file.
	pointerFile = wd + "pointer.dep" 
	args = [filledDEMFile, pointerFile] 
	pluginHost.runPlugin("FlowPointerD8", args, False, True)
	
	# Perform the flow accumulation operation.
	flowAccumFile = wd + "flow accumulation.dep" 
	outputType = "number of upslope grid cells" 
	logTransformOutput = "False" 
	args = [pointerFile, flowAccumFile, outputType, logTransformOutput] 
	pluginHost.runPlugin("FlowAccumD8", args, False, True)
	
	# Extract the streams
	streamsFile = wd + "streams.dep" 
	channelThreshold = "1000.0" 
	backValue = "NoData" 
	args = [flowAccumFile, streamsFile, channelThreshold, backValue] 
	pluginHost.runPlugin("ExtractStreams", args, False, False) # This final result will be displayed
	
	'''
	Note that in each of the examples above, I have created new variables
	to hold each of the input parameters for the plugin tools. I've done
	this more for clarity more than anything else. The script could be
	substantially shorted if the shorter variables were directly entered 
	into the args array. For instance, I could have easily used:
	
	args = [flowAccumFile, streamsFile, "1000.0", "NoData"] 
	
	for the last runPlugin and saved myself declaring the two variables.
	Because the file names are generally used in subsequent operations,
	it is a good idea to dedicate variables to those parameters.
	'''
except Exception, e:
	print e
	pluginHost.showFeedback("Error during script execution.")
	''' alternatively, you many want to send the exception to 
	the pluginHost.logException() method '''
finally:
	print "I'm done!"
#	pluginHost.updateProgress(0)
	
