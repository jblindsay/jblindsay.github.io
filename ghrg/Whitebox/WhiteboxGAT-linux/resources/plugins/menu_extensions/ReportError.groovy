import java.io.File
import whitebox.utilities.FileUtilities
import whiteboxgis.WhiteboxGui

try {
	// The following lines are necessary for the script 
	// to be recognized as a menu extension.
	def parentMenu = "Help"
	def menuLabel = "Report An Error"

	def versionNum = ((WhiteboxGui)(pluginHost)).versionNumber
	
	def logDirectory = pluginHost.getLogDirectory()
	
	def logFiles = FileUtilities.findAllFilesWithExtension(logDirectory, "xml", false)
	
	StringBuilder errorLogs = new StringBuilder()
	for (int i = logFiles.size() - 1; i >= 0; i--) { //String str : logFiles) {
		def logFile = new File(logFiles.get(i)) // str)
		errorLogs.append("\nLOG FILE: " + logFile.getName() + "\n")
		logFile.eachLine { line -> errorLogs.append(line.replace("<","!LT!").replace(">","!GT!").replace("\"", "!QUOTE!") + "\n") }
		errorLogs.append("\n")
	}
	
	String username = System.getProperty("user.name")
	
	StringBuilder sb = new StringBuilder()
	sb.append("<!DOCTYPE html>")
	sb.append("<html>")
	sb.append("<head><meta charset=\"UTF-8\">")
	sb.append("<title>Error Report</title>")
	sb.append("<style media=\"screen\" type=\"text/css\">")
	sb.append("table {margin-left: 15px;} ")
	sb.append("form {text-align: center;} ")
	sb.append("h1 {font-size: 14pt; margin-left: 15px; margin-right: 15px; text-align: center; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;} ")
	sb.append("p {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
	sb.append("table {font-size: 12pt; font-family: Helvetica, Verdana, Geneva, Arial, sans-serif;}")
	sb.append("table th {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #dedede; }")
	sb.append("table td {border-width: 1px; padding: 8px; border-style: solid; border-color: #666666; background-color: #ffffff; }")
	sb.append("caption {font-family: Helvetica, Verdana, Geneva, Arial, sans-serif; margin-left: 15px; margin-right: 15px;} ")
	sb.append(".numberCell { text-align: right; }")             
	sb.append("</style> </head> <body>")
	sb.append("<h1>Error Report</h1>")
	sb.append("<p>Describe the details of what you were doing when you encountered the error (e.g. What tool were you using? What type of data were you processing? Was there an error message, and if so, what did it say? etc.) in the text box below:</p><br>")
	sb.append("<form onsubmit=\"return true\" action=\"http://www.uoguelph.ca/cgi-bin/FormMail.pl\"")
	sb.append("method=\"post\" name=\"form1\" id=\"form1\">")
	sb.append("<input value=\"${versionNum}\" name=\"version\" type=\"hidden\">")
	sb.append("<textarea name=\"errorDescription\" form=\"form1\" rows=\"15\" cols=\"40\"></textarea><br>")
	sb.append("<p>The contents of your log files will be attached.</p><br>")
	sb.append("<input value=\"").append(errorLogs.toString()).append("\" name=\"errorlogs\" type=\"hidden\">")
	sb.append("<input value=\"Clear\" type=\"reset\"><input value=\"Submit\" type=\"submit\">") 
	sb.append("<input value=\"jlindsay@uoguelph.ca\" name=\"recipient\" type=\"hidden\">")
	sb.append("<input value=\"Whitebox GAT error report\" name=\"subject\" type=\"hidden\">")
	sb.append("<input name=\"realname\" value=\"").append(username).append("\" type=\"hidden\">")
	sb.append("<input value=\"").append(username).append("\" name=\"username\" type=\"hidden\">")
	sb.append("<input name=\"required\" value=\"errorDescription, realname\" type=\"hidden\"> </div>")
	sb.append("</form>")
	sb.append("</body></html>")
	
	pluginHost.returnData(sb.toString())

} catch (Exception e) {
    pluginHost.showFeedback(e.getMessage())
    pluginHost.logException("Error in ReportError", e)
}
