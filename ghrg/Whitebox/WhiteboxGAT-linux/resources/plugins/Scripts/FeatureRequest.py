'''
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
'''
 
from java.io import File
from java.lang import StringBuilder

# The following four variables are required for this 
# script to be integrated into the tool tree panel. 
name = "FeatureRequest"
descriptiveName = "New Feature Request"
description = "Send a request for a new feature (requires internet connection)"
toolboxes = ["topmost"]

# The following lines are necessary for the script 
# to be recognized as a menu extension.
parentMenu = "Help"
menuLabel = "New Feature Request"

try:
	sb = StringBuilder()
	sb.append("<!DOCTYPE html>")
	sb.append("<html>")
	sb.append("<head><meta charset=\"UTF-8\">")
	sb.append("<title>New Feature Request</title>")
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
	sb.append("<h1>New Feature Request</h1>")
	sb.append("<p>You can help to improve Whitebox GAT by suggesting ways of making it ")
	sb.append("even better. Describe the new feature or tool that you would like to ")
	sb.append("see developed in the text-box below. At present, we can only accept ")
	sb.append("feature requests written in English.</p><br>")
	sb.append("<form onsubmit=\"return true\" action=\"http://www.uoguelph.ca/cgi-bin/FormMail.pl\"")
	sb.append("method=\"post\" name=\"form1\" id=\"form1\">")
	sb.append("<div align=\"center\">")
	sb.append("Enter your name: <input name=\"realname\" type=\"text\"><br><br>")
	sb.append("Enter your email: <input name=\"email\" type=\"text\"><br><br>")
	sb.append("<textarea name=\"featureRequest\" form=\"form1\" rows=\"10\" cols=\"40\"></textarea><br><br>")
	sb.append("<input value=\"Clear\" type=\"reset\"><input value=\"Submit\" type=\"submit\">")
	sb.append("<input value=\"jlindsay@uoguelph.ca\" name=\"recipient\" type=\"hidden\">")
	#sb.append("<input name=\"redirect\" value=\"http://www.uoguelph.ca/~hydrogeo/Whitebox/FeatureRequestRedirect.html\" type=\"hidden\">")
	sb.append("<input value=\"Whitebox GAT feature request\" name=\"subject\" type=\"hidden\">")
	sb.append("<input name=\"required\" value=\"email,realname\" type=\"hidden\"></div>")
	sb.append("</form></body></html>")
	sb.append("</form>")
	sb.append("</body></html>")
	
	pluginHost.returnData(sb.toString())

except Exception, value:
	print value
