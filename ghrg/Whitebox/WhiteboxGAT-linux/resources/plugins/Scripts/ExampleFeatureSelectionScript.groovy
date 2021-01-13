import whitebox.geospatialfiles.ShapeFile
import whitebox.geospatialfiles.shapefile.attributes.AttributeTable
import whitebox.interfaces.MapLayer
import whitebox.geospatialfiles.VectorLayerInfo

/* Notes: 
 *  1. 'pluginHost' in the script below is the name of the WhiteboxGui
 *  2. I've displayed the vector layer 'TM_WORLD_BORDERS_SIMPL-0.2.shp'
 *     from the samples directory for this example.
 */

String myVectorsName = "TM_WORLD_BORDERS_SIMPL-0.2";

// get a list of all the displayed MapLayers
ArrayList<MapLayer> openLayers = pluginHost.getAllMapLayers();

VectorLayerInfo myVli;

// find the vector layer named 'TM_WORLD_BORDERS_SIMPL-0.2'
for (MapLayer layer : openLayers) {
	if (layer instanceof VectorLayerInfo) {
		VectorLayerInfo vli = (VectorLayerInfo) layer;
		if (vli.getShapefile().getShortName().equals(myVectorsName)) {
			myVli = vli;
		}
	}
}

if (myVli != null) {
	/* Make sure that the vli is set as the active map layer.
	 * Selections will only render for the active map layer.
	 */
	pluginHost.setActiveMapLayer(myVli.getOverlayNumber());
	
	// go through the attributes and select all countries with names starting with 'c'
	ShapeFile myShapeFile = myVli.getShapefile();
	AttributeTable myTable = myShapeFile.getAttributeTable();
	String myFieldName = "NAME";
	myVli.clearSelectedFeatures();
	for (int r = 0; r < myShapeFile.getNumberOfRecords(); r++) {
		String countryName = (String)myTable.getValue(r, myFieldName);
		if (countryName.toLowerCase().trim().startsWith("c")) {
			myVli.setSelectedFeature(r + 1); // the selected feature is a reference to the shapefile feature (base 1)
		}
	}
	
	// refresh the map
	pluginHost.refreshMap(true);

	// Now retrieve a list of selected features (work in the other direction)
	// Note that the indices are feature numbers and are therefore base 1
	ArrayList<Integer> listOfSelectedIndices = myVli.getSelectedFeatureNumbers();
	for (Integer i : listOfSelectedIndices) {
		String countryName = (String)myTable.getValue(i - 1, myFieldName);
		println("Feature Number: $i, Country: $countryName");
	}
}

println("Done");
