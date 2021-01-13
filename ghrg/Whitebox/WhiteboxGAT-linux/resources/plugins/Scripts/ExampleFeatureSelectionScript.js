var ShapeFile = Java.type('whitebox.geospatialfiles.ShapeFile');
var AttributeTable = Java.type('whitebox.geospatialfiles.shapefile.attributes.AttributeTable');
var MapLayer = Java.type('whitebox.interfaces.MapLayer');
var VectorLayerInfo = Java.type('whitebox.geospatialfiles.VectorLayerInfo');

/* Notes: 
 *  1. 'pluginHost' in the script below is the name of the WhiteboxGui
 *  2. I've displayed the vector layer 'TM_WORLD_BORDERS_SIMPL-0.2.shp'
 *     from the samples directory for this example.
 */

var myVectorsName = "TM_WORLD_BORDERS_SIMPL-0.2";

// get a list of all the displayed MapLayers
var openLayers = pluginHost.getAllMapLayers();

var myVli;

// find the vector layer named 'TM_WORLD_BORDERS_SIMPL-0.2'
for (l = 0 ; l < openLayers.size(); l++) {
	var layer = openLayers[l]
	if (layer instanceof VectorLayerInfo) {
		if (layer.getShapefile().getShortName().equals(myVectorsName)) {
			myVli = layer;
		}
	}
}

if (myVli != null) {
	/* Make sure that the vli is set as the active map layer.
	 * Selections will only render for the active map layer.
	 */
	pluginHost.setActiveMapLayer(myVli.getOverlayNumber());
	
	// go through the attributes and select all countries with names starting with 'c'
	var myShapeFile = myVli.getShapefile();
	var myTable = myShapeFile.getAttributeTable();
	var myFieldName = "NAME";
	myVli.clearSelectedFeatures();
	for (r = 0; r < myShapeFile.getNumberOfRecords(); r++) {
		var countryName = myTable.getValue(r, myFieldName);
		if (countryName.toLowerCase().trim().startsWith("c")) {
			myVli.setSelectedFeature(r + 1); // the selected feature is a reference to the shapefile feature (base 1)
		}
	}
	
	// refresh the map
	pluginHost.refreshMap(true);

	// Now retrieve a list of selected features (work in the other direction)
	// Note that the indices are feature numbers and are therefore base 1
	var listOfSelectedIndices = myVli.getSelectedFeatureNumbers();
	for (i = 0; i < listOfSelectedIndices.size(); i++) {
		var index = listOfSelectedIndices.get(i);
		var countryName = myTable.getValue(index - 1, myFieldName);
		print("Feature Number: " + index + ", Country: " + countryName);
	}
}

print("Done");
