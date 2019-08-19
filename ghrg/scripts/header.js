function insertHeader(page) {
  var imageFolder = "img";
  if (page == "whitebox") {
    imageFolder = "../img";
  }
  if (page == "whiteboxtools") {
    imageFolder = "../img";
  }
  var el = document.getElementById("titleContainer");
  if (el != null) {
    el.innerHTML = "";
    s = "<a href=\"http://www.uoguelph.ca\"><img id=\"UoGLogo\" src=\"" + imageFolder + "/UGlogo.jpg\" alt=\"UofG conerstone\" width=\"72\" height=\"72\"/></a>" + eol +
        "<img id=\"titleRunner\" src=\"" + imageFolder + "/headerRunner.png\" alt=\"title runner image\"/>" + eol +
        "<div id=\"headerImageBox\"></div>" + eol +
        "<span id=\"johnText\">John Lindsay<small>, PhD</small></span>" + eol +
        "<span id=\"johnSubtext\"><a href=\"http://www.uoguelph.ca/geography/\">Geography, Environment &amp; Geomatics</a>&#8729;<a href=\"http://www.uoguelph.ca/\">University of Guelph</a></span>" + eol +
        "<span id=\"labText\">Geomorphometry &amp; Hydrogeomatics Research Group</span>";
    el.innerHTML += s;
  }
}

function insertFooter() {
  var el = document.getElementById("footer");
  if (el != null) {
    el.innerHTML = "";
    s = "<span id=\"copyrightText\" class=\"footerText\">Copyright &copy; John Lindsay 2015-2019</span>" + eol +
        "<span id=\"footerLinkBar\" class=\"footerText\"><a href=\"http://www.uoguelph.ca/geography/\">Department of Geography, Environment &amp; Geomatics</a> &#8729; <a href=\"http://www.uoguelph.ca/\">The University of Guelph</a> &#8729; <a href=\"https://whiteboxgeospatial.wordpress.com\">Whitebox Blog</a></span>";
    el.innerHTML += s;
  }
}

function getNavBar(page) {
  var el = document.getElementById("navBar");
  if (el != null) {
    el.innerHTML = "";
    var homeClass = "";
    var publicationsClass = "";
    var softwareClass = "";
    var whiteboxClass = "";
    var whiteboxtoolsClass = "";
    var researchGroupClass = "";
    var opportunitiesClass = "";
    var newsClass = "";
    var pageFolder = "";
    var whiteboxFolder = "./Whitebox/";
    var whiteboxtoolsFolder = "./WhiteboxTools/";
    switch (page) {
      case "home":
        homeClass = "class=\"activated\"";
        break;
      case "publications":
        publicationsClass = "class=\"activated\"";
        break;
      case "software":
        softwareClass = "class=\"activated\"";
        break;
      case "whitebox":
        softwareClass = "class=\"activated\"";
        whiteboxClass = "class=\"activated\"";
        pageFolder = "../";
        whiteboxFolder = "";
        whiteboxtoolsFolder = "../WhiteboxTools/";
        break;
      case "whiteboxtools":
        softwareClass = "class=\"activated\"";
        whiteboxtoolsClass = "class=\"activated\"";
        pageFolder = "../";
        whiteboxFolder = "../Whitebox/";
        whiteboxtoolsFolder = "";
        break;
      case "researchGroup":
        researchGroupClass = "class=\"activated\"";
        break;
      case "opportunities":
        opportunitiesClass = "class=\"activated\"";
        break;
      case "news":
        newsClass = "class=\"activated\"";
        break;
      // default:
      //   console.log("Could not create navigation bar due to unrecognized page.");
    }
    s = "<ul>" + eol +
        "  <li><a href=\"" + pageFolder + "index.html\"" + homeClass + ">HOME</a></li>" + eol  +
        "  <li><a href=\"" + pageFolder + "publications.html\"" + publicationsClass + ">PUBLICATIONS</a></li>" + eol  +
        "  <li><a href=\"" + pageFolder + "software.shtml\"" + softwareClass + ">SOFTWARE</a>" + eol  +
        "  <ul>" + eol +
        "    <li " + whiteboxClass + "><a href=\"" + whiteboxFolder + "index.html\">Whitebox GAT</a></li>" + eol  +
        "    <li " + whiteboxtoolsClass + "><a href=\"" + whiteboxtoolsFolder + "index.html\">WhiteboxTools</a></li>" + eol  +
        "    <li><a href=\"" + pageFolder + "software.shtml#GoSpatial\">GoSpatial</a></li>" + eol  +
        "  </ul>" + eol  +
        "  </li>" + eol  +
        "  <li><a href=\"" + pageFolder + "research_group.html\"" + researchGroupClass + ">RESEARCH GROUP</a></li>" + eol  +
        "  <li><a href=\"" + pageFolder + "opportunities.html\"" + opportunitiesClass + ">OPPORTUNITIES</a></li>" + eol  +
        "  <li><a href=\"" + pageFolder + "news.html\"" + newsClass + ">NEWS</a></li>" + eol  +
        "</ul>";
    el.innerHTML += s;
  }
}

var OSName = "Unknown";
if (window.navigator.userAgent.indexOf("Windows") != -1) OSName="Windows";
if (window.navigator.userAgent.indexOf("Mac")!=-1) OSName="Mac/iOS";
if (window.navigator.userAgent.indexOf("X11")!=-1) OSName="UNIX";
if (window.navigator.userAgent.indexOf("Linux")!=-1) OSName="Linux";
var eol = "\n";
if (OSName === "Windows") {
  eol = "\r\n";
}

// News
var newsJson = {
  newsItems:
  [
    {
      date: "2019-08-19",
      headline: "Anthony Francioni Publishes in Remote Sensing",
      fullStory: "Congratulations to now graduated GHRG member Anthony Francioni for his co-authorship on a new \
      article in the journal <em>Remote Sensing</em>. The article, titled \"LiDAR DEM smoothing and the preservation \
      of drainage features\" is now out in the latest issue of the journal. "
    },
    {
      date: "2019-08-12",
      headline: "New lab publication in Remote Sensing",
      fullStory: "Kevin Roberts has had his GHRG-based Masters research published as a new article in \
      Remote Sensing. His article, titled \"An Analysis of Ground-Point Classifiers for Terrestrial LiDAR\" \
      has been accepted for publication in the 10th Anniversary Issue of Remote Sensing. Well done Kevin!"
    },
    {
      date: "2019-06-23",
      headline: "Successful Defences in the GHRG",
      fullStory: "Congratulations to Simon Gudim, Kevin Roberts, and Anthony Francioni for each \
      successfully defending their Masters research. Each of you did an incredible job writing your \
      theses and presenting your fantastic research. I'm excited to see what comes next for each \
      of you. Well done guys!"
    },
    {
      date: "2019-05-24",
      headline: "WhiteboxTools 0.16.0 now available",
      fullStory: "The WhiteboxTools development team is proud to announce the \
      the WhiteboxTools 0.16.0. This release largely focusses on bug-fixes \
      rather than the addition of new features. The following tools were added \
      to the platform:<br> \
        MergeLineSegments \
        SphericalStdDevOfNormals tools. \
      <li>Fixed a bug with reading LAS files with point records with extra bytes. Previously, the LAS decoder \
        assumed the Point Record Length matched that of the LAS specifications (with the variable of the \
        optional intensity and user data). Some LAS files in the wild (particularly those created using \
        LASTools and of LAS version 1.2) have larger Point Record Lengths, which presumably carry extra \
        bytes of information. These extra byes are ignored, but they no longer throw off the decoding.</li> \
      <li>Fixed a bug with writing Big-Ending GeoTIFF files. The 'MM' file header was not correct previously.</li> \
      <li>Significantly reduced the memory requirements of the StochasticDepressionAnalysis tool. The tool \
        may be somewhat slower as a result, but it should be applicable to larger DEMs than was previously \
        possible.</li> \
      <li>Fixed bugs with the Union and SplitWithLines tools.</li> \
      <li>WhiteboxTools can now read and write Shapefiles of MultiPointZ, PolyLineZ, and PolygonZ ShapeTypes \
        missing the optional 'M' values (i.e. measures).</li> \
      <li>SelectTilesByPolygon and LidarTileFootprint are now compatible with LAZ file inputs. Both of these \
        tools only rely on information in the input LiDAR file's header, which is the same for a LAZ file \
        as a LAS file.</li> \
      <li>Fixed a bug with writing Saga GIS files (*.sdat) that inverted rasters.</li>"
    },
    {
      date: "2019-05-15",
      headline: "2018 Best Paper Award",
      fullStory: "Congratulations to former GHRG member Melanie Chabot for being awarded the Canadian Remote Sensing Society's 2018 Best Paper Award, for her Canadian Journal of Remote Sensing paper  <a href=\"https://www.tandfonline.com/doi/abs/10.1080/07038992.2018.1461559\"><em>Comparing the Use of Terrestrial LiDAR Scanners and Pin Profilers for Deriving Agricultural Roughness Statistics</em></a>. Melanie will receive the award at the 40th Canadian Symposium on Remote Sensing and Geomatics Atlantic 2019 in Fredericton, NB in June. Well done Melanie!"
    },
    {
      date: "2019-05-15",
      headline: "Dan Newman Awarded OGS",
      fullStory: "Congratulations to GHRG member Dan Newman for being awarded an Ontario Graduate Scholarship (OGS). \
      Dan's doctoral research into hyperscale geomorphometric landscape characterization is state-of-the-science and \
      this prestigious award is well-deserved recognition. "
    },
    {
      date: "2019-01-08",
      headline: "WhiteboxTools 0.13.0 now available",
      fullStory: "The WhiteboxTools development team is proud to announce the \
      the first release of 2019. This release largely focusses on bug-fixes \
      rather than the addition of new features. The following tools were added \
      to the project:<br> \
          MosaicWithFeathering \
      <ul><li>Support was added for GeoTIFF MODELTRANSFORMATIONTAG (Tag 33920).</li> \
      <li>Support was added for reading GeoTIFFs that have coordinate transformations \
        defined by multiple tiepoints contained with the ModelTiepointTag (Tag 33922). \
        These rasters have their raster-to-model transform defined by a 2D polynomial \
        regression of the 3rd order.</li> \
      <li>The initialize_using_file function in the abstract Raster model now transfers \
        information contained in an input GeoTIFF's ModelTiePoint, ModelPixelScale, \
        ModelTransform, GeoKeyDirectory, GeoDoubleParms, and GeoAsciiParams tags to \
        the output raster. This means that if a GeoTIFF file is input to a Whitebox  \
        tool, and the output raster is specified to be of GeoTIFF format as well, \
        all of the coordinate information contain in the input raster will now be \
        contained in the output raster.</li> \
      <li>The FeaturePreservingDenoise and DrainagePreservingSmoothing tools, both of \
        which are used for DEM generalization, now represent surface normal vectors  \
        using 32-bit floats instead of the original double-precision values. This  \
        does not alter the results of these tools significantly, but does reduce the  \
        memory requirements and run-times of these tools substantially.</li> \
      <li>The LidarKappa tool now outputs a raster displaying the spatial distribution  \
        of the overall accuracy per grid cell (i.e. percent agreement).</li> \
      <li>Fixed a bug with the RasterStreamsToVector tool that resulted in overlapping \
        traced streams.</li> \
      <li>The D8FlowAccumulation tool has been modifed to use a fixed flow-width to  \
        calculate specific contributing area, equal to the average grid cell resolution.  \
        The tool previously used a variable flow-width for SCA calculations, however, \
        1. this differs from the constant value used in Whitebox GAT, and 2. a  \
        variable flow-width means that flow accumulation values do not increase  \
        continuously in a downstream direction. This last issue was causing problems \
        with applications involving stream network extraction. This change does not \
        affect the 'cells' nor 'catchment area' outputs of the tool.</li> \
      <li>Fixed a bug with the GeoTIFF NoData tag.</li> \
      <li>Fixed a bug with the SetNodataValue tool.</li></ul> \
      "
    },
    {
      date: "2018-11-22",
      headline: "WhiteboxTools 0.12.0 now available",
      fullStory: "This release of the WhiteboxTools library is marked by the \
      addition of several useful vector overlay and data processing tools. Most \
      notably, this includes support for TINing and TIN based gridding (vector \
      and LiDAR), as well as several vector patch shape indicies. The following \
      tools were added to the project:<br>\
      BlockMaximumGridding <br>\
      BlockMinimumGridding <br>\
      Clip <br>\
      Dissolve <br>\
      Erase <br>\
      JoinTables <br>\
      Intersect <br>\
      LasToShapefile <br>\
      LidarClassifySubset <br>\
      LinearityIndex <br>\
      LineIntersections <br>\
      LongestFlowpath <br>\
      MergeTableWithCsv <br>\
      MergeVectors <br>\
      NearestNeighbourGridding <br>\
      PatchOrientation <br>\
      Polygonize <br>\
      RasterToVectorLines <br>\
      SplitWithLines <br>\
      SymmetricalDifference <br>\
      Union <br>\
      VoronoiDiagram <br>\
      <p>Other changes include:</p>\
      <ul><li>Modified the algorithm used by the CostDistance tool from an iterative method of \
      finding the minimum cost surface to one that uses a priority-flood approach. This \
      is far more efficient. Also, there was a bug in the original code that was the \
      result of a mis-match between the neighbouring cell distances and the back-link\
      direction. In some cases this resulted in an infinite loop, which is now resolved.\
      <li>Improvements have been made to the WhiteboxTools GeoTIFF reader. A bug has been \
        fixed that prevented tile-oriented (in contrast to the more common strip-oriented) \
        TIFF files from being read properly. Support has been added for reading rasters \
        that have been compressed using the DEFLATE algorithm. Lastly, the WhiteboxTools \
        GeoTIFF reader now supports sparse rasters, as implemented by GDAL's GeoTIFF driver.</li> \
      <li>An issue in the SAGA raster format reader has been fixed.</li></ul> \
      <p>Please see the readme.txt associated with the release for further details.</p>"
    },

    {
      date: "2018-10-13",
      headline: "Open-access book \"Geocomputation with R\" now available",
      fullStory: "Researcher at the Leeds Institute for Transport Studies and \
      friend of the GHRG, Robin Lovelace, and his co-authors Jakub \
      Nowosad, and Jannes Muenchow, have recently released a wonderful new resource for \
      geocomputation. The book, titled, \"<em>Geocomputation with R</em>\" is open-access \
      and can be found <a href=\"https://geocompr.robinlovelace.net/\">on Robin's \
      blog site</a>."
    },

    {
      date: "2018-10-01",
      headline: "WhiteboxTools 0.11.0 now available",
      fullStory: "This release of the WhiteboxTools library is marked by the \
      addition of several useful vector data processing capabilities. Most \
      notably, this includes support for TINing and TIN based gridding (vector \
      and LiDAR), as well as several vector patch shape indicies. The following \
      tools were added to the project:<br>\
      AddPointCoordinatesToTable<br> \
      CentroidVector<br> \
      CompactnessRatio<br> \
      ConstructVectorTIN<br> \
      ElongationRatio<br> \
      ExtendVectorLines<br> \
      HoleProportion<br> \
      LayerFootprint<br> \
      LidarConstructVectorTIN<br> \
      LidarTINGridding<br> \
      LinesToPolygons<br> \
      Medoid<br> \
      MinimumBoundingCircle<br> \
      MinimumBoundingEnvelope<br> \
      MultiPartToSinglePart<br> \
      PerimeterAreaRatio<br> \
      PolygonArea<br> \
      PolygonPerimeter<br> \
      RasterStreamsToVector<br> \
      RasterToVectorPoints<br> \
      RelatedCircumscribingCircle<br> \
      RemovePolygonHoles<br> \
      ShapeComplexityIndex<br> \
      SinglePartToMultiPart<br> \
      SmoothVectors<br> \
      SumOverlay<br> \
      TINGridding<br> \
      <br>Please see the readme.txt associated with the release for further details."
    },

    {
      date: "2018-09-18",
      headline: "WhiteboxTools 0.10.0 now available",
      fullStory: "We are pleased to announce the release of the latest version of \
      WhiteboxTools. This release sees the addition of several vector and LiDAR \
      analysis tools and numerous bug fixes and enhancement. The following tools \
      have been added to WhiteboxTools in this release:<br> \
      CreateHexagonalVectorGrid<br> \
      CreateRectangularVectorGrid<br> \
      DrainagePreservingSmoothing<br> \
      EliminateCoincidentPoints<br> \
      ExtractNodes<br> \
      HighPassMedianFilter<br> \
      LasToMultipointShapefile<br> \
      LidarHexBinning and VectorHexBinning<br> \
      LidarTileFootprint<br> \
      MaxDifferenceFromMean<br> \
      MinimumBoundingBox<br> \
      MinimumConvexHull<br> \
      PolygonLongAxis and PolygonShortAxis<br> \
      PolygonsToLines<br> \
      ReinitializeAttributeTable<br><br> \
      Refactoring of some data related to Point2D, and common algorithms (e.g. \
        point-in-poly, convex hull).<br> \
      Added unit tests to BoundingBox, point_in_poly, convex_hull, and elsewhere. \
      Fixed a bug in LiDAR join related to tiles with fewer than two points. LAS files \
        now issue a warning upon saving when they contain less than two points. \
      The default callback can now be modified in whitebox_tools.py, such that \
        a single custom callback can be used without having to specify it for each \
        tool function call.<br> \
      Added initial support for getting projection ESPG and WKT info from LAS files  \
        and GeoTiff data. This is the start of a more fullsome approach to handling \
        spatial reference system information in the library.<br> \
      Fixed a bug in saving Shapefile m and z data.<br> \
      Fixed a bug that wouldn't allow the LidarIdwInterpolation and  \
        LidarNearestNeighbourGridding tool to interpolate point classification data.<br> \
      LidarGroundPointFilter now has the ability to output a classified LAS file rather  \
        than merely filtering non-ground points. Ground points are assigned classification \
        values of 2 while non-ground points are classified as 1.<br> \
      Updated the LidarKappaIndex tool to use a NN-search to find matching points between \
        the compared point clouds.<br> \
      Modified the FixedRadiusSearch structure to use 64-bit floats for storing coordinates. \
        This impacts performance efficiency but is needed for the fine precision of  \
        positional information found in terrestrial LiDAR data. FixedRadiusSearch structures \
        have also had approximate kNN search methods added.<br>"
    },

    {
      date: "2018-08-22",
      headline: "Geomorphometry 2018",
      fullStory: "GHRG research on <a href=\"http://2018.geomorphometry.org/geomorphometery_2018_program.html\"> \
      hyper-scale geomorphometric analysis</a> has \
      just been presented at <a href=\"http://2018.geomorphometry.org/\">Geomorphometry 2018</a> \
      in Boulder, CO, USA.  The paper can be found on <a href=\"https://peerj.com/preprints/27110/\">PeerJ</a>. \
      Overall, the conference was a tremendous success, thanks in large part to \
      the tireless efforts of Peter Guth, Scott Peckham, his wife, Lucy, and the \
      rest of the organizing committee. The community is looking forward to \
      Geomorphometry 2020!"
    },

    {
      date: "2018-07-28",
      headline: "GHRG member Dan Newman publishes research",
      fullStory: "Congratulations to <a href=\"research_group.html#DanNewman\">Dan \
      Newman</a> for publishing the second paper from his geomorphometry-focused MSc research, \
      <a href=\"publications.html#NewmanGeosciencePaper\">\
      Measuring hyperscale topographic anisotropy as a continuous landscape property</a> \
      in <em>Geomorphology</em>."
    },

    {
      date: "2018-05-30",
      headline: "WhiteboxTools 0.8.0 now available",
      fullStory: "The following tools have been added to WhiteboxTools in this release:<br> \
    CornerDetection<br> \
    FastAlmostGaussianFilter<br> \
    GaussianContrastStretch<br> \
    IdwInterpolation<br> \
    ImpoundmentIndex<br> \
    LidarThin<br> \
    StochasticDepressionAnalysis<br> \
    UnsharpMasking<br> \
    WeightedOverlay<br><br> \
Modified some filters to take RGB inputs by operating on the intensity value.  \
  These include AdaptiveFilter, BilateralFilter, ConservativeSmoothingFilter,  \
  DiffOfGaussianFilter, EdgePreservingMeanFilter, EmbossFilter, GaussianFilter,  \
  HighPassFilter, KNearestMeanFilter, LaplacianFilter, LaplacianOfGaussianFilter,  \
  LeeFilter, MaximumFilter, MeanFilter, MedianFilter, MinimumFilter, OlympicFilter,  \
  PrewittFilter, RangeFilter, RobertsCrossFilter, ScharrFilter, SobelFilter, and \
  UserDefinedWeightsFilter.<br><br> \
Fixed a bug with reading/writing Whitebox Raster files containing RGB data.<br><br> \
Modified the MajorityFilter tool to improve efficiency substantially. Also fixed \
  a bug in it and the DiversityFilter tools."
    },

    {
      date: "2018-05-01",
      headline: "WhiteboxTools 0.7.0 now available",
      fullStory: "This release of WhiteboxTools includes the addition of several \
      important tools, including:<br> \
        AttributeCorrelation<br> \
        ChangeVectorAnalysis<br> \
        ClassifyOverlapPoints<br> \
        ClipLidarToPolygon<br> \
        ClipRasterToPolygon<br> \
        CorrectVignetting<br> \
        ErasePolygonFromLidar<br> \
        ExportTableToCsv<br> \
        RaiseWalls<br> \
        TrendSurface<br> \
        TrendSurfaceVectorPoints<br> \
        UnnestBasins<br> \
        UserDefinedWeightsFilter<br><br> \
      Additionally, the TraceDownslopeFlowpaths has been updated to take vector \
      seed point inputs."
    },

    {
      date: "2018-04-22",
      headline: "WhiteboxTools 0.6.0 now available",
      fullStory: "This release of WhiteboxTools includes a couple of significant \
      feature additions, including the ability to read Shapefile attribute data \
      (.dbf files), and support for reading LZW compressed GeoTIFFs. The LZW \
      decoder can also handle the use of a horizontal predictor (TIFF Tag 317). \
      The following tools have been added in this release:<br> \
          AttributeHistogram<br> \
          AttributeScattergram<br> \
          CountIf<br> \
          ListUniqueValues<br> \
          VectorLinesToRaster<br> \
          VectorPointsToRaster<br> \
          VectorPolygonsToRaster"
    },

    {
      date: "2018-04-04",
      headline: "WhiteboxTools 0.5.0 now available",
      fullStory: "This release of the WhiteboxTools library includes new tools, \
      improvements to the User Manual, and other minor bug fixes. The new tools \
      include the following:<br> \
      EdgePreservingMeanFilter<br> \
      ElevationAboveStreamEuclidean<br> \
      ErasePolygonFromRaster<br> \
      FillBurn<br> \
      FlattenLakes<br> \
      ImageStackProfile<br> \
      InPlaceAdd<br> \
      InPlaceDivide<br> \
      InPlaceMultiply<br> \
      InPlaceSubtract<br> \
      MaxAnisotropyDevSignature<br> \
      PrincipalComponentAnalysis<br> \
      RasterizeStreams"
    },

    {
      date: "2018-03-30",
      headline: "GHRG member Dan Newman publishes research",
      fullStory: "Congratulations to <a href=\"research_group.html#DanNewman\">Dan \
      Newman</a> for publishing the first paper from his geomorphometry-focused MSc research, \
      <a href=\"publications.html#NewmanGeomorphologyPaper\">\
      Evaluating metrics of local topographic position for multiscale geomorphometric \
      analysis</a> in <em>Geomorphology</em>. Dan's Masters research into multiscale \
      geomorphometric analysis is truly at the forefront of the field. This is remarkably well done Dan!"
    },

    {
      date: "2018-03-11",
      headline: "Former GHRG member Melanie Chabot publishes research",
      fullStory: "Congratulations to <a href=\"research_group.html#MelChabot\">Melanie \
      Chabot</a> for publishing her MSc research <a href=\"publications.html#ChabotCJRSPaper\">\
      LiDAR derived roughness statistics</a> in the Canadian Journal of Remote Sensing. \
      Well done Mel!"
    },

    {
      date: "2018-03-04",
      headline: "WhiteboxTools 0.4.0 now available",
      fullStory: "This release of the WhiteboxTools library includes significant \
      improvments to the Python scripting interface, the User Manual, and of course, \
      several new and interesting tools, including:<br> \
        LidarColourize, <br> \
        LidarPointStats,<br> \
        LidarRemoveDuplicates,<br> \
        LongProfile,<br> \
        LongProfileFromPoints,<br> \
        MaxElevDevSignature,<br> \
        MultiscaleRoughness,<br> \
        MultiscaleRoughnessSignature,<br> \
        PrintGeoTiffTags,<br> \
        Profile"
    },
    {
      date: "2018-02-15",
      headline: "WhiteboxTools 0.3.1 now available",
      fullStory: "This release of the WhiteboxTools library includes many refinements to LiDAR data \
      processing and a fix for a file-selection bug on Windows. Pre-compiled \
      binary files for the WhiteboxTools geospatial analysis engine are now available \
      from the Geomorphometry and Hydrogeomatics Research Group \
      <a href=\"http://www.uoguelph.ca/~hydrogeo/software.shtml#WhiteboxTools\">Software</a> \
      web site."
    },
    {
      date: "2018-01-12",
      headline: "WhiteboxTools 0.2.0 now available",
      fullStory: "The second public release of the WhiteboxTools library is now available. Pre-compiled \
      binary files for the WhiteboxTools geospatial analysis engine are now available \
      from the Geomorphometry and Hydrogeomatics Research Group \
      <a href=\"http://www.uoguelph.ca/~hydrogeo/software.shtml#WhiteboxTools\">Software</a> \
      web site."
    },
    {
      date: "2017-12-12",
      headline: "WhiteboxTools compiled binary files now available",
      fullStory: "Pre-compiled binary files for the WhiteboxTools geospatial analysis engine are now available \
      from the Geomorphometry and Hydrogeomatics Research Group <a href=\"http://www.uoguelph.ca/~hydrogeo/software.shtml#WhiteboxTools\">Software</a> \
      web site."
    },
    {
      date: "2017-07-12",
      headline: "Announcing the WhiteboxTools library",
      fullStory: "<a href=\"https://github.com/jblindsay/whitebox-geospatial-analysis-tools/tree/master/whitebox_tools\">WhiteboxTools</a> \
      is a sub-project of the Whitebox GAT open-source GIS project. Over the years, there have been a large number of requests\
      to call Whitebox GAT tools and functionality from outside of the user-interface (e.g. from Python automation scripts).\
      WhiteboxTools is intended to meet these usage requirements. WhiteboxTools is an advanced geospatial data analysis engine.\
      Eventually most of the tools contained within Whitebox GAT will be ported to the new WhiteboxTools library, which is being\
      developed using the Rust programing language and is compiled to native code. In addition to separating the user-interface\
      and data processing components of Whitebox GAT, this migration should significantly improve data processing efficiency. This is\
      the result of the highly efficient Rust code, which is generally faster than Java, and because many of the WhiteboxTools\
      functions are designed to process data concurrently (in parallel) wherever possible. The official public announcement can be\
      found <a href=\"https://whiteboxgeospatial.wordpress.com/2017/07/10/announcing-the-whiteboxtools-library/\">here</a> at the Whitebox blog. You may track the progress of tool\
      porting efforts <a href=\"https://github.com/jblindsay/whitebox-geospatial-analysis-tools/blob/master/whitebox_tools/tool_porting.md\">\
      here</a> at the Github code repository. This announcement marks a major milestone in the history of the Whitebox GAT project\
      and represents a very significant commitment to future developments in GIS software."
    },
    {
      date: "2017-04-28",
      headline: "The Whitebox project featured in recent article",
      fullStory: "The Whitebox project has been featured in a recent article on the\
      <a href=\"http://www.gisresources.com/university-guelph-prof-john-lindsay-develops-whitebox-geospatial-analysis-tools-processing-geospatial-data/\">\
      <em>GIS Resources</em></a> geospatial news site. <em>GIS Resources</em> is a global platform for\
      the latest information for the geospatial industry, focusing on recent developments\
      in geospatial science and technology."
    },
    {
      date: "2017-04-28",
      headline: "New GHRG undergraduate independent studies (W17)",
      fullStory: "Congratulations to Anthony Francioni (Assessing the techniques for removal\
      of small-scale noise in LiDAR DEMs) and and Libby George (<em>The effects of data\
      reduction on digital elevation models</em>) for\
      successfully completing geomorphometry related undergratuate independent studies\
      this past winter semester. We are delighted to say that Anthony will be joining\
      the GHRG this fall as well to complete an MSc."
    },
    {
      date: "2017-03-02",
      headline: "How are you using Whitebox GAT?",
      fullStory: "The Whitebox GAT user community is fantastic and we're always receiving feedback\
      from users about their experiences using the software. We're always keenly interested in hearing\
      about all the wonderful ways that people are applying Whitebox GAT either for education, research,\
      or in industry and government settings. If Whitebox GAT has helped you in some interesting way please send\
      your comments to Prof. John Lindsay (jlindsay@uoguelph.ca) in an email with the subject heading\
      '<em>My Whitebox story</em>'. We're looking forward to hearing from you."
    },
    {
      date: "2017-02-10",
      headline: "Whitebox GAT survey shows global distribution of usage",
      fullStory: "A post has recently been published on the\
      <a href=\"https://whiteboxgeospatial.wordpress.com\">Whitebox Blog</a> describing\
      an analysis of the Whitebox GAT user community by analyzing download data. You can\
      read the full blog post <a href=\"https://whiteboxgeospatial.wordpress.com/2017/02/09/whitebox-gat-usage/\">\
      here</a>. The post included a version of the map below, which shows each of the 5149 cities within\
      178 countries in which Whitebox GAT has been downloaded. We're definitely international!<br><br>\
      <a href=\"img/WhiteboxDownloads2.png\">\
      <img src=\"img/WhiteboxDownloads2.png\" width=\"596\" height=\"auto\" alt=\"Whitebox downloads\"/></a>\
      (click image to enlarge)"
    },
    {
      date: "2017-02-04",
      headline: "Kevin Roberts wins best poster prize at CGU Eastern Section 2017 Student Meeting",
      fullStory: "Congratulations to UofG Geography undergrad Kevin Roberts for winning the\
      best poster prize at the Canadian Geophysical Union (CGU) Eastern Section 2017 Student\
      Meeting. Kevin's poster, entitled,\
      <a href=\"http://www.uoguelph.ca/~hydrogeo/pubs/RobertsAndLindsayCGU.pdf\"><em>A comparison\
      of ground-point separation methods applied to terrestrial laser scanner mapped LiDAR point\
      clouds</em></a>, was based on his undergraduate independent study (GEOG4990), which was carried\
      out during the Fall 2016 semester.<br><br>\
      <img src=\"img/KevinCGUPoster.png\" width=\"596\" height=\"auto\" alt=\"Kevin Roberts wins CGU\
      conference best poster\"/>"
    },

    {
      date: "2017-01-30",
      headline: "Prof. Lindsay gives talk to the Waterloo Water Institute",
      fullStory: "Prof. Lindsay will be giving a talk on Jan. 31st to the University of Waterloo\
      Water Institute and Department of Civil and Evironmental Engineering entitled,\
      <a href=\"https://uwaterloo.ca/water-institute/events/using-open-access-gis-address-issues-spatial-hydrological\">\
      <em>Using Open-Access GIS to Address Issues in Spatial Hydrological Modelling</em></a>. A pdf of\
      the presentation can be found\
      <a href=\"http://jblindsay.github.io/GuestTalks/WaterlooWaterInstitute2017/OAGISandIssuesInSpatialHydrology.pdf\">here</a>."
    },

    {
      date: "2017-01-25",
      headline: "Whitebox GAT version 3.4.0 Montreal released",
      fullStory: "Whitebox GAT version 3.4.0 (codename Montreal) has been\
      <a href=\"http://www.uoguelph.ca/~hydrogeo/Whitebox/index.html\">released</a> for Windows,\
      MacOS, and Linux. This release adds several new features, tools, and bug fixes. Several new\
      LiDAR data processing tools have been added. Additionally, native-compiled plugin tools\
      are now supported.\
      <br><br><img src=\"img/screenshots/WbGAT_Montreal.png\" width=\"596\" height=\"auto\" alt=\"Whitebox 3.4.0 Montreal\"/>"
    },

    {
      date: "2016-12-16",
      headline: "New GHRG undergraduate independent studies",
      fullStory: "Well done to Kevin Roberts and Jordan McDonald, for completing geomorphometry based undergraduate independent studies (GEOG4990)."
    },

    {
      date: "2016-07-12",
      headline: "Former GHRG member Kat Woodrow publishes research",
      fullStory: "Congratulations to <a href=\"research_group.html#KatWoodrow\">Kat Woodrow</a> for publishing her MSc research on <a href=\"publications.html#WoodrowPaper\">DEM conditioning techniques</a> in the <a href=\"http://www.sciencedirect.com/science/article/pii/S0022169416304486\">Journal of Hydrology</a>. Kat now works for Geomorphix. Well done!"
    },

    {
      date: "2016-07-11",
      headline: "New Whitebox GAT journal article",
      fullStory: "A <a href=\"publications.html#WhiteboxPaper\">journal article</a> describing the Whitebox GAT project has been accepted for publication in <a href=\"http://www.sciencedirect.com/science/article/pii/S0098300416301820\">Computers &amp; Geosciences</a>."
    },

    {
      date: "2016-06-22",
      headline: "Former GHRG member Colleen Fuss publishes research",
      fullStory: "Congratulations to <a href=\"research_group.html#ColleenFuss\">Colleen Fuss</a> for publishing her MSc research on <a href=\"publications.html#FussPaper\">DEM fusion methods</a> in the <a href=\"http://www.tandfonline.com/doi/abs/10.1080/17538947.2016.1208685\">International Journal of Digital Earth</a>. Colleen also started a new job with the Canada Centre for Mapping and Earth Observation."
    },

    {
      date: "2016-05-20",
      headline: "Whitebox GAT version 3.3.0 Glasgow released",
      fullStory: "Whitebox GAT version 3.3.0 (codename Glasgow) is now available for <a href=\"http://www.uoguelph.ca/~hydrogeo/Whitebox/index.html\">download</a>. See <a href=\"https://whiteboxgeospatial.wordpress.com/2016/05/20/whitebox-gat-3-3-0-has-been-released/\">here</a> for the release announcement."
    }
  ]
}


// function getNews() {
//     var result = "";
//     $.ajax({
//       type: "GET",
//       url:"./data/news",
//       async: false,
//       success:function(data) {
//          result = JSON && JSON.parse(data) || $.parseJSON(data);
//       },
//       error: function(XMLHttpRequest, textStatus, errorThrown) {
//           // alert("Status: " + textStatus); alert("Error: " + errorThrown);
//           result = "No recent updates.";
//       }
//    });
//    return result;
// }

function insertNewsItems(num_stories_to_show = 5, whitebox_site = false) {
  // var el = document.getElementById("news");
  // if (el != null) {
  //   el.innerHTML = "";
  //   var s = "";
  //   var num_stories_to_show = 5;
  //   var n = newsItems.length < num_stories_to_show ? newsItems.length : num_stories_to_show;
  //   for (var i = 0; i < n; i++) {
  //       s += "  <p>" + newsItems[i] + "</p>";
  //   }
  //   s += "  <p>Read older news posts <a href=\"news.html\">here</a>.</p>";
  //   el.innerHTML += s;
  // }

  // createFullNews();

  // $.ajax({
  //   url:"./data/news",
  //   // async: false,
  //   success:function(data) {
  //     //  newsJson = JSON && JSON.parse(data) || $.parseJSON(data);
  //      var el = document.getElementById("newsStories");
  //      if (el != null) {
  //        var s = "";
  //        el.innerHTML = "";
  //       //  var num_stories_to_show = 5;
  //       //  var n = newsJson.newsItems.length < num_stories_to_show ? newsJson.newsItems.length : num_stories_to_show;
  //       //  for (var i = 0; i < n; i++) {
  //       //    s += "<p>" + newsJson.newsItems[i].headline + "</p>";
  //       //  }
  //        s += "  <p>Read older news posts <a href=\"news.html\">here</a>.</p>";
  //        el.innerHTML += s;
  //      }
  //   },
  //   error: function(XMLHttpRequest, textStatus, errorThrown) {
  //       // alert("Status: " + textStatus); alert("Error: " + errorThrown);
  //       var el = document.getElementById("newsStories");
  //       if (el != null) {
  //         var s = "";
  //         el.innerHTML = "";
  //         s += "  <p>Read older news posts <a href=\"news.html\">here</a>.</p>";
  //         el.innerHTML += s;
  //       }
  //   }
  // });

  var s = "";
  var el = document.getElementById("news");
  if (el != null) {
    // var s = "";
    el.innerHTML = "";
    // var num_stories_to_show = 5;
    var n = newsJson.newsItems.length < num_stories_to_show ? newsJson.newsItems.length : num_stories_to_show;
    for (var i = 0; i < n; i++) {
      var ds = new Date(newsJson.newsItems[i].date);
      if (!newsJson.newsItems[i].headline.endsWith('?') && !newsJson.newsItems[i].headline.endsWith('?')) {
        if (!whitebox_site) {
          s += "<p>" + newsJson.newsItems[i].headline + ". (" + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + " Full story <a href=\"news.html#story" + (i+1) + "\">here</a>)</p>";
        } else {
          s += "<p>" + newsJson.newsItems[i].headline + ". (" + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + " Full story <a href=\"../news.html#story" + (i+1) + "\">here</a>)</p>";
        }
      } else {
        if (!whitebox_site) {
          s += "<p>" + newsJson.newsItems[i].headline + " (" + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + " Full story <a href=\"news.html#story" + (i+1) + "\">here</a>)</p>";
        } else {
          s += "<p>" + newsJson.newsItems[i].headline + " (" + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + " Full story <a href=\"../news.html#story" + (i+1) + "\">here</a>)</p>";
        }
      }
    }
    if (!whitebox_site) {
      s += "  <p>Read older news posts <a href=\"news.html\">here</a>.</p>";
    } else {
      s += "  <p>Read older news posts <a href=\"../news.html\">here</a>.</p>";
    }
    el.innerHTML += s;
  }
}

// function createFullNews() {
//   var el = document.getElementById("newsStories");
//   if (el != null) {
//     var s = "";
//     el.innerHTML = "";
//     var d = new Date();
//     var year = d.getFullYear();
//     for (var j = year; j > 2015; j -= 1) {
//       var num_stories = 0;
//       for (var i = 0; i < newsItems.length; i++) {
//         if (newsItems[i].indexOf("/" + String(j) + ")") > 0) {
//           num_stories += 1
//         }
//       }
//       if (num_stories >= 0) {
//         s += "<h3>" + j + "</h3>";
//         s += "<ul>";
//         for (var i = 0; i < newsItems.length; i++) {
//           if (newsItems[i].indexOf("/" + String(j) + ")") > 0) {
//             s += "  <li>" + newsItems[i] + "</li>";
//           }
//         }
//         s += "</ul>";
//       }
//     }
//
//     el.innerHTML += s;
//   }
//   // console.log("I'm here.");
// }

// function addJQuerySource() {
//   var script = document.createElement('script');
//   script.src = 'https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js';
//   script.type = 'text/javascript';
//   document.getElementsByTagName('head')[0].appendChild(script);
//   document.write("Hey there");
// }

// (function(){
//   var newscript = document.createElement('script');
//     //  newscript.type = 'text/javascript';
//     //  newscript.async = false;
//     //  newscript.src = 'https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js';
//     //  newscript.content = "document.write(\"I'm here now\");"
//   // (document.getElementsByTagName('head')[0]||document.getElementsByTagName('body')[0]).appendChild(newscript);
//   // document.getElementsByTagName('head')[0].appendChild(newscript);
//   document.getElementsByTagName('head')[0].appendChild("<script src=\'https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js\'></script>")
// })();

function scrollTo(hash) {
    location.hash = "#" + hash;
}

function createFullNews() {

  var el = document.getElementById("newsStories");
  if (el != null) {
    var d = new Date();
    var year = d.getFullYear();
    var s = "";
    el.innerHTML = "";
    var storyNum = 1;
    s += "<h2>" + year + "</h2>";

    for (var i = 0; i < newsJson.newsItems.length; i++) {
      var ds = new Date(newsJson.newsItems[i].date);
      var storyYear = ds.getFullYear();
      if (storyYear != year) {
        year = storyYear;
        s += "<h2>" + year + "</h2>";
      }
      if (storyNum % 2 == 1) {
        s += "<div class=\"newsLight\">"
        if (newsJson.newsItems[i].headline != null) {
          s += "<a id=\"story" + storyNum + "\" href=\"\"></a>"  + "<h3><em>" + newsJson.newsItems[i].headline + "</em></h3>";
        }
        storyNum++;
        s += "<p>" + newsJson.newsItems[i].fullStory + "<br>Date: " + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + "</p>";
        s += "</div>"
      } else {
        s += "<div class=\"newsDark\">"
        if (newsJson.newsItems[i].headline != null) {
          s += "<a id=\"story" + storyNum + "\" href=\"\"></a>"  + "<h3><em>" + newsJson.newsItems[i].headline + "</em></h3>";
        }
        storyNum++;
        s += "<p>" + newsJson.newsItems[i].fullStory + "<br>Date: " + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + (ds.getDate() + 1) + "</p>";
        s += "</div>"
      }
    }
    el.innerHTML += s;
  }

  // if (!newsJson) {
  //   newsJson = getNews();
  // }
  // document.write(newsJson.newsItems[0].fullStory);
  // var el = document.getElementById("newsStories");
  // if (el != null) {
  //   var d = new Date();
  //   var year = d.getFullYear();
  //   var s = "";
  //   el.innerHTML = "";
  //   var storyNum = 1;
  //   s += "<h2>" + year + " Stories </h2>";
  //
  //   for (var i = 0; i < newsJson.newsItems.length; i++) {
  //     var ds = new Date(newsJson.newsItems[i].date);
  //     var storyYear = ds.getFullYear();
  //     if (storyYear != year) {
  //       year = storyYear;
  //       s += "<h2>" + year + " Stories </h2>";
  //     }
  //     if (newsJson.newsItems[i].headline != null) {
  //       s += "<h3 name=\"story" + storyNum + "\">" + newsJson.newsItems[i].headline + "</h3>";
  //     }
  //     storyNum++;
  //     s += "<p>" + newsJson.newsItems[i].fullStory + "<br>Date: " + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + ds.getDate() + "</p>";
  //   }
  //   el.innerHTML += s;
  // }

  // $.ajax({
  //   url:'./data/news',
  //   success: function (data) {
  //     var el = document.getElementById("newsStories");
  //     if (el != null) {
  //       var d = new Date();
  //       var year = d.getFullYear();
  //       var s = "";
  //       el.innerHTML = "";
  //       var storyNum = 1;
  //       s += "<h2>" + year + " Stories </h2>";
  //       obj = JSON && JSON.parse(data) || $.parseJSON(data);
  //
  //       for (var i = 0; i < obj.newsItems.length; i++) {
  //         var ds = new Date(obj.newsItems[i].date);
  //         var storyYear = ds.getFullYear();
  //         if (storyYear != year) {
  //           year = storyYear;
  //           s += "<h2>" + year + " Stories </h2>";
  //         }
  //         if (obj.newsItems[i].headline != null) {
  //           s += "<h3 name=\"story" + storyNum + "\">" + obj.newsItems[i].headline + "</h3>";
  //         }
  //         storyNum++;
  //         s += "<p>" + obj.newsItems[i].fullStory + "<br>Date: " + ds.getFullYear() + "-" + (ds.getMonth() + 1) + "-" + ds.getDate() + "</p>";
  //       }
  //       el.innerHTML += s;
  //     }
  //   },
  //   error: function(XMLHttpRequest, textStatus, errorThrown) {
  //       alert("Status: " + textStatus); alert("Error: " + errorThrown);
  //   }
  // });

}
