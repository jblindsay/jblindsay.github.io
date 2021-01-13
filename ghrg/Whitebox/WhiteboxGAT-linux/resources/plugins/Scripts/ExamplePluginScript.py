'''
This script is an example of a plugin script, complete with 
integration within Whitebox's tool tree panel, and a dialog
for retrieving user-specified parameters. The main work of 
the tool is performed on a seperate thread that is able to 
communitcate with the main Whitebox user interface, e.g. by 
updating progress and responding to the user's selection of 
the cancel button.

The script is not intended to be run, so much as it is meant 
for demonstrating how you can maniputlate geographical data 
within Whitebox. Python is an excellent language to use to 
quickly develop a plugin tool in Whitebox. It's particularly 
useful if you simply want to automate a workflow by calling
other existing tools. One issue with the use of Python however
is that it does not have the performance that you get with 
Java and Groovy. A good alternative is to create your script
using the Groovy language with static compilation. This gets
you performance very nearly matching that of Java while still
getting the benefits of a scripting type language's terse 
syntax, ideal for rapid plugin tool development. You'll find
plenty of examples of plugins developed using Groovy.

Notice that for a script to appear in the Whitebox toolbox,
it must be saved within the Scripts folder within the 
Plugins and resources directories. You will also need to
relaunch Whitebox before it appears in the toolbox. Any changes
that you make to the script will be live however, meaning that
you can modify a script and then run it immediately afterwards.
'''

import os
from threading import Thread
from whitebox.ui.plugin_dialog import ScriptDialog
from java.awt.event import ActionListener
from whitebox.geospatialfiles import WhiteboxRaster
from whitebox.geospatialfiles.WhiteboxRasterBase import DataType
from whitebox.geospatialfiles import ShapeFile
from whitebox.geospatialfiles.shapefile import ShapeType
from whitebox.geospatialfiles.shapefile import ShapeFileRecord
from whitebox.geospatialfiles.shapefile import Geometry
from whitebox.geospatialfiles.shapefile.attributes import DBFField
from whitebox.geospatialfiles.shapefile.attributes import AttributeTable

'''The following four variables are required for this 
script to be integrated into the tool tree panel. 
Comment them out if you want to remove the script 
from the tool tree panel, uncomment them if you want 
the script to appear as an integrated plugin tool.
You would need to relaunch Whitebox for the tool to
appear in the tool tree panel.'''
#name = "MyScript" # the name variable is the unique identifier for this tool. It cannot contain any spaces.
descriptiveName = "My Script" # the descriptiveName variable is the name of the tool as it appears in the various tool listings. This name can contain spaces.
description = "This tool processes counts the number of records in a shapefile." # this string appears as a brief description in the status bar when the user clicks on your tool.
toolboxes = ["ConversionTools", "StatisticalTools"] # these are the various toolboxes in the tool tree panel that this tool will appear in.
	
class MyScript(ActionListener):
	def __init__(self, args):
		if len(args) != 0:
			''' This will be followed if the parameters have already 
			been set, e.g. if the tool was called from another script. '''
			self.execute(args)
		else:
			''' The script dialog is created here. The first parameter 
			points the dialog to the Whitebox user interface. The second 
			parameter is the title of the dialog. In most cases this should 
			be set to the descriptive name of the tool. The third parameter
			is a listener for the OK button. When the user presses OK on 
			the dialog, actionPerformed method of this class will be informed.''' 
			self.sd = ScriptDialog(pluginHost, "My Script", self)	
			
			''' Sets the name of the helpfile, which should be an html file 
			contained in the Help directory. Ideally, it has the same name as 
			the tool's name.'''
			helpFile = self.__class__.__name__
			self.sd.setHelpFile(helpFile)
	
			''' This next line will enable the 'View Code' button on the dialog.'''
			self.sd.setSourceFile(os.path.abspath(__file__))
	
			''' Now add some components to the dialog. The following is a 
			description of some of the most commonly used dialog components.'''
	
			''' Use the addDialogFile method to add a File input/output component 
			to the dialog. The method has six parameters: 1) String description, 
			2) String labelText, 3) String dialogMode ("open" or "save"), 
			4) String filter (commonly "Raster Files (*.dep), DEP", 
			"Vector Files (*.shp), SHP", or "Whitebox Files (*.dep; *.shp), DEP, SHP"), 
			5) boolean showButton, 6) boolean makeOptional.'''
			self.sd.addDialogFile("Input raster file", "Input Raster File:", "open", "Raster Files (*.dep), DEP", True, False)
			self.sd.addDialogFile("Output rasterfile", "Output Raster File:", "save", "Raster Files (*.dep), DEP", True, False)
			self.sd.addDialogFile("Input vector file", "Input Vector File:", "open", "Vector Files (*.shp), SHP", True, False)
			self.sd.addDialogFile("Output vector file", "Output Vector File:", "save", "Vector Files (*.shp), SHP", True, False)
			
			''' Use the addDialogMultiFile method to add a multi-file selection 
			component to the dialog. This method has three parameters: 1) String description, 
			2) String labelText, 3) String filter (commonly "Raster Files (*.dep), DEP", 
			"Vector Files (*.shp), SHP", or "Whitebox Files (*.dep; *.shp), DEP, SHP")'''
			#self.sd.addDialogMultiFile("This is a multifile input", "Choose some raster files:", "Raster Files (*.dep), DEP")
			
			''' Use the addDialogDataInput method to add a data input component 
			retrieving text strings or numbers from the user. This method has 
			five parameters: 1) String description, 2) String labelText,
	        3) String initialText, 4) boolean numericalInput (text strings will 
	        not be accepted), 5) boolean makeOptional.'''
			self.sd.addDialogDataInput("Enter data here", "Enter a number", "1.3", True, False)
			
			''' Use the addDialogCheckBox method to add a checkbox to the dialog. 
			This method has three parameters: 1) String description, 2) String 
			labelText, 3) boolean initialState'''
			self.sd.addDialogCheckBox("Do you want it to rain?", "Will it rain?", True)
	
			''' Use the addDialogComboBox method to add a drop-down selection 
			list to the dialog. This method has four parameters: 1) String 
			description, 2) String labelText, 3) String[] (or ArrayList<>) listItems, 
			4) int defaultItem'''
			#self.sd.addDialogComboBox("Weather type", "Weather type:", ["rain", "snow", "sun", "hail"], 1)
	
			''' Use the addDialogOption method to add a radio-button style 
			option on the dialog. This method has four parameters: 1) String 
			description, 2) String labelText, 3) String button1String, 4) String 
			button2String'''
			#self.sd.addDialogOption("This is an option", "Choose one:", "Plants", "Animals")
	
			''' Use the addDialogFieldSelector method to add a field selector component 
			to the dialog. A field selector is a combination file input, specifically 
			for selecting shapefiles, and a list selector for specifying a field 
			within the shapefile's attribute table. This method has three parameters: 
			1) String description, 2) String labelText, 3) boolean allowMultipleSelection'''
			#self.sd.addDialogFieldSelector("A field selector", "Choose a field", True)
	
			''' Use the addDialogLabel method to add a label to the dialog. This 
			method has one parameter: 1) String text'''
			#self.sd.addDialogLabel("Model Parameters:")

			''' Resize the dialog to the standard size and display it '''
			self.sd.setSize(800, 400)
			self.sd.visible = True


	''' This method is the main part of your script. It is were the 
	main processing occurs and will be called on its own thread when 
	the user presses the OK button on the dialog (see actionPerformed
	below).'''
	def execute(self, args):
		try:
			''' Make sure that the args array has the expected number 
			of parameters. If the script is being called by another 
			script, this may not be the case, at which point, the script 
			should inform the user of the error and then end gracefully.'''
			if len(args) != 6:
				pluginHost.showFeedback("Incorrect number of arguments given to tool.")
				return

			''' Read in the parameters. All parameters are provided as
			strings and must be converted to their appropriate form 
			for use, e.g. floats etc. '''
			inputfile = args[0]
			outputfile = args[1]
			inputshapefile = args[2]
			outputshapefile = args[3]
			numericalData = float(args[4])
			booleanValue = bool(args[5])
			
			''' The following is an example of how to call another 
			Whitebox plugin from your script. simply create a list of
			strings to hold the parameters you wish to provide the 
			plugin. Then call the runPlugin method of the pluginHost.
			pluginHost is a reference to the Whitebox user interface.
			The last two boolean variables in the runPlugin method are 
			to indicate 1) whether the tool should be run on a dedicated 
			thread, which is not advisable if you would like to run 
			through a workflow sequentially (i.e. complete first task then 
			move onto second task) and 2) whether upon completion, any 
			data that would be automatically returned by the tool (e.g. 
			a displayed image) should be suppressed.''' 
			# args2 = [inputfile, outputfile, "0.0"]
			# pluginHost.runPlugin("FillDepressions", args2, False, True)

			''' The following code sets up an input Whitebox raster 
			for reading. It then reads the number of rows and columns
			and the nodata value.'''
			inputraster = WhiteboxRaster(inputfile, 'r')
			rows = inputraster.getNumberRows()
			cols = inputraster.getNumberColumns()
			nodata = inputraster.getNoDataValue()

			''' The next line of code will create a new WhiteboxRaster 
			of the name 'ouptutfile'. It will be the same dimensions 
			and extent as the input file. The last parameter is the 
			initial value of the raster grid, here set to nodata.'''
			outputraster = WhiteboxRaster(outputfile, "rw", 
  		  	  inputfile, DataType.FLOAT, nodata)
			
			''' This code can be used to scan through a raster grid. '''
			oldprogress = -1
			for row in xrange(0, rows):
				for col in xrange(0, cols):
					z = inputraster.getValue(row, col)
					if z != nodata:
						outputraster.setValue(row, col, z * 1000)
					else:
						outputraster.setValue(row, col, nodata)
				
				progress = (int)(100.0 * row / (rows - 1))
				if progress > oldprogress:
					oldprogress = progress
					pluginHost.updateProgress(progress)
				if pluginHost.isRequestForOperationCancelSet():
					pluginHost.showFeedback("Operation cancelled")
					return
			
			inputraster.close()
			outputraster.close()

			''' Call the returnData method of pluginHost. You can 
			return the file names of raster files and shapefiles 
			and they will be automatically displayed. If you return 
			a text string (other than a file name) it will be displayed 
			in the textbox at the bottom of the Whitebox user interface. 
			If you return html, it will be rendered in a window. If you 
			return a JPanel, it will be placed in a JDialog and displayed.'''
			pluginHost.returnData(outputfile)


			''' The following code is an example of how to work with 
			vector data.'''

			''' Open and existing shapefile'''
			input = ShapeFile(inputshapefile)
			shapetype = input.getShapeType()
			table = input.getAttributeTable()

			''' reading from the records in a shapefile '''
			r = 0
			for record in input.records:
				shapegeometry = record.getGeometry()
				''' Read the points into an array '''
				points = shapegeometry.getPoints()
				# x = points[pointnum][0]
				# y = points[pointnum][1]

				''' polylines and polygons can contain multiple
				parts. The 'parts' array can be used to identify
				the starting node from the points array for each
				part in the geometry. '''
				parts = shapegeometry.getParts()

				''' the following shows how you read the record's
				attributes and update an entry'''
				# recData = table.getRecord(r)
                # recData[recData.length - 1] = new Double(1.0)
                # table.updateRecord(r, recData)
                # r += 1

            
			''' The following code can be used to create a new shapefile. '''
			''' First, set up the fields within the attribute table. '''
			field1 = DBFField()
			field1.setName("FID")
			field1.setDataType(DBFField.DBFDataType.NUMERIC)
			field1.setFieldLength(10)
			field1.setDecimalCount(0)
			field2 = DBFField()
			field2.setName("NAME")
			field2.setDataType(DBFField.DBFDataType.STRING)
			field2.setFieldLength(20)
			fields = [field1, field2]

			''' Now create the shapefile, feeding the constructor 
			the fields array. '''
			output = ShapeFile(outputshapefile, ShapeType.POLYLINE, fields)

			# outputGeom = Point(x, y)
			# outputGeom = PolyLine(parts, points)
			# outputGeom = Polygon(parts, points)
			''' Note that for ShapeTypes of M or Z dimensions, 
			there are alternative constructors that allow for 
			setting the M and Z data using a double array. '''
			# Object[] rowData = new Object[2]
			# rowData[0] = (float)FID
			# rowData[1] = "New feature"
			# output.addRecord(outputGeom, rowData);
			
  			output.write()
			
		except Exception, e:
			print e
			pluginHost.showFeedback("Error during script execution.")
			''' alternatively, you many want to send the exception to 
			the pluginHost.logException() method '''
			return
		finally:
			pluginHost.updateProgress(0)


	''' This method is used as a listener for the dialog. If 
	the user presses the OK button, it will collect the 
	user-specified parameters then send them to the execute 
	method above. '''
	def actionPerformed(self, event):
		if event.getActionCommand() == "ok":
			args = self.sd.collectParameters()
			t = Thread(target=lambda: self.execute(args))
			t.start()
		#elif event.getActionCommand() == "close":
		#	print "The user cancelled the operation."

if args is None:
	''' Your script will have a variable called 'pluginHost' injected 
	into it. This variable is the Whitebox user interface and you 
	can use it to interact with the interface to do things like 
	show the user feedback, change the progress bar, check to see 
	if the user pressed the cancel button, display a file, etc.'''
	pluginHost.showFeedback("The arguments array has not been set.")
else:		
	MyScript(args)
