/* The following lines are necessary for the script 
   to be recognized as a menu extension.
   To see this menu extension, uncomment the lines 
   and then relaunch Whitebox. */
//def parentMenu = "File"
//def menuLabel = "Say Hello"
//def acceleratorKey = "U" // adding an acceleratorKey is optional

/* This is the code that will be executed when the user
   selects the menu. In this case, it simply has Whitebox 
   show a popup window message but the script can be 
   as sophisticated as you desire. */
pluginHost.showFeedback("Hello from a menu extension")
