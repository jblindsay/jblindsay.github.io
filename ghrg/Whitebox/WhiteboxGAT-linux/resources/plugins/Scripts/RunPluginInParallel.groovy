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
 
import java.util.concurrent.Future
import java.util.concurrent.*
import whitebox.interfaces.WhiteboxPluginHost
import whiteboxgis.WhiteboxGui
import whitebox.plugins.ReturnedDataEvent
import whitebox.interfaces.ReturnedDataListener
import whitebox.ui.plugin_dialog.*
import java.awt.event.ActionListener
import java.awt.event.ActionEvent
import java.util.concurrent.atomic.AtomicInteger
import groovy.transform.CompileStatic

def name = "RunPluginInParallel";
def descriptiveName = "Run Plugin In Parallel";
def description = "Runs a plugin tool multiple times in parallel";
def toolboxes = ["topmost"];

public class RunPluginInParallel implements ActionListener {
    private WhiteboxPluginHost pluginHost
    private ScriptDialog sd
    private String descriptiveName
    private String basisFunctionType = ""
	private int numTasks = 1;
	private AtomicInteger numSolvedTasks = new AtomicInteger(0)
	private int oldProgress = -1
	
    public RunPluginInParallel(WhiteboxPluginHost pluginHost, 
        String[] args, def name, def descriptiveName) {
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
            sd.setHelpFile(name)
		
            // Specifying the source file allows the 'view code' 
            // button on the tool dialog to be displayed.
            def pathSep = File.separator
            def scriptFile = pluginHost.getResourcesDirectory() + "plugins" + pathSep + "Scripts" + pathSep + name + ".groovy"
            sd.setSourceFile(scriptFile)

			List plugNames = pluginHost.returnPluginList()
			String[] plugs = new String[plugNames.size()]
			int p = 0
			plugNames.each() {
				plugs[p] = it
				p++
			}
            // add some components to the dialog
            sd.addDialogComboBox("Which plugin tool should be run?", "Plugin:", plugs, 0)
            sd.addDialogFile("Parameter file", "Input Parameter File:", "open", "Text Files (*.txt), TXT", true, false)
            sd.addDialogCheckBox("Suppress return data", "Suppress return data?", false)
            
            // resize the dialog to the standard size and display it
            sd.setSize(800, 400)
            sd.visible = true
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

    @CompileStatic
    private void execute(String[] args) {

		String pluginName = args[0]
		String argsFile = args[1]
		boolean suppressReturns = Boolean.parseBoolean(args[2])

		try {
			ArrayList<String[]> pluginArgs = new ArrayList<>()
			File file = new File(argsFile)
			if (file.exists()) {
				new File(argsFile).eachLine { line -> 
					String str = String.valueOf(line)
					if (str != null && !(str.trim()).isEmpty()) {
						String[] s = str.replace("\"", "").split(",")
						for (int j in 0..(s.length - 1)) {
							s[j] = s[j].trim()
						}
						pluginArgs.add(s)
					}
				}
			} else if (argsFile.contains("\n")) {
				// The second parameter is actually a string
				// containing the parameters
				String[] lines = argsFile.split("\n")
				lines.each() { line ->
					String str = String.valueOf(line)
					if (str != null && !(str.trim()).isEmpty()) {
						String[] s = str.replace("\"", "").split(",")
						for (int j in 0..(s.length - 1)) {
							s[j] = s[j].trim()
						}
						pluginArgs.add(s)
					}
				}
			} else {
				pluginHost.showFeedback("There is something incorrect with the second parameter.")
				return
			}

			numTasks = pluginArgs.size()
			
			pluginHost.updateProgress("Please wait...", 0)
			ArrayList<DoWork> tasks = new ArrayList<>();
			for (int i in 0..(numTasks - 1)) {
				String[] pArgs = pluginArgs.get(i)
				tasks.add(new DoWork(pluginName, pArgs, suppressReturns, pluginHost))
			}
			
			/* If each of the operating plugins are simutaneously updating 
			 their progress, it will make collisions that are bad. This is
			 handled by disabling the progress bar, which must later be undone.
			 */
			((WhiteboxGui)pluginHost).setUpdateProgressEnabled(false)
		
			ReturnListener rl = new ReturnListener()
			((WhiteboxGui)pluginHost).addReturnedDataEventListener(rl)
		
			((WhiteboxGui)pluginHost).isPluginReturnDataSuppressed(true)
			
			ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
			// the only reason for the getExecutorResults method 
	  	    // is that Groovy throws a compilation type mis-match
	  	    // error when compiled statically. I think it's a bug.
	  	    List<Future<Boolean>> results = getExecutorResults(executor, tasks);
			executor.shutdown();
			
//			int i = 0
//			int progress
//			int oldProgress = -1
			int numSuccessful = 0
			for (Future<Boolean> result : results) {
				Boolean data = result.get()
				if (data) { numSuccessful ++ }
//				i++
//				// update progress bar
//				progress = (int)(100f * i / numTasks)
//				if (progress > oldProgress) {
//					pluginHost.updateProgress("Progress:", progress)
//					oldProgress = progress
//				}
//				// check to see if the user has requested a cancellation
//				if (pluginHost.isRequestForOperationCancelSet()) {
//					pluginHost.showFeedback("Operation cancelled")
//					return
//				}
			}
		
			((WhiteboxGui)pluginHost).removeReturnedDataEventListener(rl)


			// re-enable returned data.
			((WhiteboxGui)pluginHost).isPluginReturnDataSuppressed(false)
				
			if (!suppressReturns) { 
				ArrayList returns = rl.getReturns()
				for (Object obj : returns) {
					pluginHost.returnData(obj)
				}
				((WhiteboxGui)pluginHost).removeReturnedDataEventListener(rl)
			}
			
			if (numSuccessful == numTasks) {
				pluginHost.showFeedback("Operations Complete. All of the operations were successful.")
			} else {
				pluginHost.showFeedback("Operations Complete. ${numSuccessful} of the ${numTasks} operations were successful.")
			}
		
		} catch (Exception e) {
			((WhiteboxGui)pluginHost).isPluginReturnDataSuppressed(false)
			pluginHost.showFeedback("An error has occurred during operation. See log file for details.")
			pluginHost.logException("Error in " + descriptiveName, e)
		} finally {
			// re-enable returned data.
			((WhiteboxGui)pluginHost).isPluginReturnDataSuppressed(false)
			
			// re-enable the progress bar.
			((WhiteboxGui)pluginHost).setUpdateProgressEnabled(true)
		
			pluginHost.updateProgress("Progress:", 0)
			
		}
    }

    public List<Future<Boolean>> getExecutorResults(ExecutorService executor, ArrayList<DoWork> tasks) {
    	List<Future<Boolean>> results = executor.invokeAll(tasks);
		return results
    }


	class DoWork implements Callable<Boolean> {
		private String[] args
		private String pluginName
		private WhiteboxPluginHost pluginHost
		private boolean suppressReturns = true
		
	    public DoWork(String pluginName, String[] args, boolean suppressReturns, WhiteboxPluginHost pluginHost) {
	    	this.args = args
	    	this.pluginName = pluginName
	        this.pluginHost = pluginHost
	        this.suppressReturns = suppressReturns
	   	}
	    	
	    @Override
	    public Boolean call() {
	    	try {
	    		if (pluginHost.isRequestForOperationCancelSet()) { return }
	    		pluginHost.runPlugin(pluginName, args, false)
		    	int solved = numSolvedTasks.incrementAndGet()
				int progress = (int) (100f * solved / numTasks)
				if (progress != oldProgress) {
					// re-enable the progress bar.
					((WhiteboxGui)pluginHost).setUpdateProgressEnabled(true)
		
					pluginHost.updateProgress("Solved ${solved} tasks of ${numTasks}:", progress)

					// disable the progress bar.
					((WhiteboxGui)pluginHost).setUpdateProgressEnabled(false)
		
				}
		    	return Boolean.TRUE
	    	} catch (Exception e) {
	    		return false
	    	}
	    }                
	}
	
	class ReturnListener implements ReturnedDataListener {
		private ArrayList returns = new ArrayList()
		
		public void dataReturned(ReturnedDataEvent evt) {
			returns.add(evt.getData())
		}
	
		public ArrayList getReturns() {
			return returns
		}
	}
}

if (args == null) {
    pluginHost.showFeedback("Plugin arguments not set.")
} else {
    def f = new RunPluginInParallel(pluginHost, args, name, descriptiveName)
}
