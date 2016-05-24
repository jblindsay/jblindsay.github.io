Whitebox GAT Introduction
=========================
Download
--------
The latest version of Whitebox GAT can be downloaded [here](http://www.uoguelph.ca/~hydrogeo/Whitebox/download.shtml "Whitebox GAT 'Iguazu' v. 3.2.2 (60 MB)")

Installation
------------
Whitebox will be downloaded as a compressed zip file. In order to run Whitebox, the folder must first be unzipped (decompressed). Once the file is unzipped, to launch the program simply locate and run the WhiteboxGIS.jar file.

In order to run Whitebox you must also ensure that you have installed and enabled an up to date version of [Java](https://java.com/en/download/ "Java download") (version 8.0 or higher).

A Brief Overview
----------------
![Whitebox Breakdown](/Images/Image1.png "Whitebox Breakdown")

The Whitebox display is divided up into 7 main sections, the Menu Bar, the Toolbar, the Sidebar, the Search Bar, the Progress Bar, the Cartographic Layout Toolbar, and the Map Display Area.

**Menu Bar**<br>
You will notice that the main menu bar contains 6 drop down menus, *File*, *Data Layers*, *View*, *Cartographic*, *Tools*, and *Help*. The *File* menu contains basic functions for dealing with recent data layers, maps, and working directories. The *Data Layers* menu allows you to add/remove data layers, modify the order in which they appear, and many other basic functions dealing with data layers. The *view* menu allows you to change the way in which you view your map through functions like zooming or panning. The *Cartographic* menu allows you to open a map to work on, save a current map, as well as add various visuals to your map such as a title, scale bar, north arrow, etc. The *Tools* menu includes tools such as the measurement tool, raster calculator, as well as tools to write your own script and help files for Whitebox. The *Help* menu contains an Index which outlines what each tool does as well as provides an example of the script used by the tool.

**Toolbar**<br>
The toolbar located directly below the Menu Bar contains many of the same tools as the Menu Bar. The tools located within the toolbar are more easily accessible through their representative icons. Icons within the toolbar are broken up into 5 sections, representing tools from 5 of the 6 drop down menus of the main menu. The first section on the left contains tool from the Data Layers menu, the second from the View menu, the third from the Cartographic menu, the fourth form the Tools menu, and the fifth from the Help menu.

**Sidebar**<br>
On the left-hand side you will notice a sidebar made up of 3 tabs. The Tools tab is made up of various folders each pertaining to a broad lens of tools. Within each folder are subfolders each lading to a more specific niche of tools. Included in every tool is a description of what the tool is used for, how the tool works, and an example of the scripting used in the tool. The Layers tab allows you to toggle the order in which your layers will be displayed, as well all allows you to switch between active maps being displayed. The Features tab displays vector features as selected by the Select Features tool.

**Search Bar**<br>
Below the sidebar you will also find a search bar which allows you to quickly perform a search within Whitebox for any tool you that you may require. Below the search bar you fill find 3 tabs. The All tab lists every tool in alphabetical order, Most Used lists some of your frequently used tools, and Recent lists some of your recently used tools.

**Progress Bar**<br>
The right-hand side of progress bar shows you the progress of whatever tool function you are performing. If you have a map displayed, you will notice as you move your cursor around the map display area that the left-hand side of the bar will display the spatial coordinates, row and column numbers, as well as z value as you move your cursor around your map.

**Cartographic Layout Toolbar**<br>
Within this toolbar you will find tools to adjust the layout of your map within the screen

**Map Display Area** <br>
Upon opening Whitebox this area will be blank. Upon adding data layers however, this is where your active map will be displayed

Getting Started
---------------

**Adding Data**<br>
Adding data layers to your map can be done in two ways. Firstly, you can select the Data Layers dropdown menu from the top menu bar and select "Add Layers to Map". Similarly, the first symbol on the left-hand side of the Toolbar can also be selected to add data layers to your map. Your Whitebox download should have included some  sample files, try adding the Vermont DEM.dep data layer to your map.

Importing from a different file format? No problem. Whitebox's import tool can convert many file types to those compatible to run in Whitebox. To simplify this process further, there is no need to even run the tool manually from the sidebar. Simply add the data as you would normally through the Add Layers to Map function, and select the file format and file to be imported and the corresponding tool be run automatically. Note data ca be added/imported in batches, rather than individually.

**Selecting the Active Map**<br>
It is possible to have multiple maps open at the same time, however only one map can be the active map. The active map is the map to which new layers and any other cartographic elements will be added to. To set the active map, simply navigate to the Layers tab of the sidebar, right click on whichever map you would like set as the active map, and click Set as Active Map.

**Navigating the Map Display Area**<br>
In the toolbar you will notice that there are 2 sets of tools which will be of use in navigating through your map display, zoom tools and the pan tool. With the Zoom-In or Zoom-Out tool selected zooming can be achieved either by clicking on the map, by using the scroll wheel on your mouse, or by clicking and holding down the right or left mouse button and dragging the cursor across the screen, which will create a black window which will be zoomed on once the mouse button is let go. The Zoom to Full Extent tool will return you to the original full view of the map extent. When zoomed in on a data layer, clicking and dragging with the Pan tool enabled allows you to navigate across your zoomed image.

![Whitebox Map Display Area Options](/Images/Image2.png "Whitebox Map Display Area Options")

**Viewing Vector Feature Attributes**<br>
The attributes of vector features are an integral component of many GIS analysis's. Viewing vector feature attributes can be achieved in two ways. The first and broader way in which this can be done is by navigating to the Layers tab of the sidebar, right clicking on the vector layer of interest, and click View Attribute Table. This will display a table of each point in the vector layer and any attribute data attributed to it. Another method for viewing vector feature attributes involves the Select Feature tool from the toolbar. When this tool is selected you will notice that the sidebar automatically switches to the Features tab. As you click on numerous polylines/polygons you will notice that the corresponding feature REC# will be displayed in the Selected Features box displayed in the sidebar. You can then click on these features individually from within the Selected Features box to view their corresponding attributes.

**Display Properties**<br>
Accessing the properties of you data layer can be integral to your analysis. With the Vermont DEM still loaded click on the Layers tab of the sidebar, then right click on Vermont DEM and select Layer Display Properties. A window with 2 tabs will be displayed. The Display tab will allow you to modify things like the palette and the opacity of the layer. The File tab displays more technical information such as the units in which x,y,z are displayed in (metres), the cell size, as well as any metadata that may be included in the data layer.

Additionally, you will notice that this popup menu contains much more than just the Layer Display Properties option. From here you can also adjust the layer's order, toggle its visibility/visibility in the legend, adjust the colour palette, export the layer to another file format, and view a histogram of the data. Histograms can be viewed in two forms, Probability Density Function (PDF), or Cumulative Distribution Function (CDF).

![Whitebox Layer Right-Click](/Images/Image3.png "Whitebox Layer Right-Click")


**Working Directories and Saving Your Map**<br>
Lets take it back to the basics for a second. When running nearly any tool at the bare minimum it will require both an input and output file. After you run the tool, a new file will be created (the output). Unlike other GIS programs, Whitebox does not require that you manually select a working directory for your new output files to be saved to, rather it  automatically saves any new outputs to the last folder from which data was added from. That being said, under the File tab in the top menu bar, Recent Working Directories can be accessed from the dropdown menu, from which you can toggle through various past directories to which you can have your outputs saved to.

Additionally, when saving multiple layers as one map, select "Save Map" under the Cartographic tab in the top menu bar to save them as a .wmap file. It should however noted that this file alone does not contain the actual map layers. A .wmap file simply contains the description and layout of the map, meaning that without all of the actual layers saved on your computer a well, nothing will be loaded from the .wmap file.

**Printing and Exporting your Map**<br>
Printing your map can be done easily through the Cartographic dropdown menu. Through this same menu your completed map can also be exported to an image file so that it can be more easily shared. Individual layers can also be exported by right-clicking on the layer in the Layers tab and selecting Export Layer. Layers can be exported to file formats compatible with other GIS software such as ArcGIS, GRASS, Idrisi, Surfer, and Saga.

**Formating your Map**<br>
When formating the layout of your map the first thing that you are going to want to do is make sure that your map is going to fit correctly onto the page that you are printing on. This can be achieved through the Draw the Page Function, which can be easily accessed by clicking the top icon in the Cartographic Layout Toolbar on the right-hand side of the screen. Now that the page that your map will be printed on is displayed, you can begin adding additional layout elements to you map.

Many features can be added to your map from the Cartographic dropdown menu in the top menu bar such as titles, north arrow, scale bar, legend, neatline, text areas, images, and additional map areas. Once one or multiple of these elements have been added to your map you will now be ale to utilize the Select Map Elements tool form the toolbar to move around and resize any of the added elements.

Once you have added some of these features to your map you will notice that they are still not attached to your map in any way, rather they are just overlaid ontop of it. Meaning that if you were to click and drag your map to reposition it, all of the features would stay in their exact same spot on the page rather than moving cohesively with the map. This can be solved by using the Group Elements tool.

When adding annotations or other elements to your map the Group Elements tool found under the Cartographic menu can be quite useful. To group multiple elements together simply hold down the shift key and click on each element that you would like to grouped together and then select the Group Elements tool from the Cartographic dropdown menu. This tool can be quite useful for labeling your map, as it will append the labels to their specific point on the map. Notice: Once your map has been grouped, it is no longer considered to be a layer, and is now purely visual. Once you have grouped your map layer in with another object you will no longer be able to zoom in on and pan around your map. Any zooming done now will simply zoom into the entire page itself as if zooming into a web page. Therefor, grouping should only be done once all other analysis are completed. However, if you group your map and realize that you still have further analysis to perform, fear not, as by once again individually selecting each grouped element and selecting Ungroup Elements from the Cartographic menu, you regain full access to your map and all of its layers.

If you just need to quickly reposition your map on the page and don't want to have to group everything together, then the Select All Map Elements is the tool for you. This will allow you to move your map and all of the features on it is unison, and will disengage once you click leaving all of your features in their new location.

Modifying text, font, borders, background colour and more can all be accessed by right clicking the map element. To then modify these elements, click on the corresponding boxes on the right-hand side and make any necessary adjustments. Additionally, double clicking anywhere in the map display area will bring up the Map Properties menu, which will allow you to scroll through all displayed map elements and make any necessary changes to them.

![Whitebox Formating](/Images/Image4.png "Whitebox Formating")
