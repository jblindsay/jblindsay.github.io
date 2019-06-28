var currentImageIndex = -1,
  maxImageIndex = 0,
  images = [],
  changeInterval = 8000;

  // prepares relevant variables to cycle through images
  var setUp = function() {
    var images2 = document.images;
    for (i = 0; i < images2.length; i += 1) {
      if (images2[i].id.indexOf("asideImage") > -1) {
        images[images.length] = images2[i];
      }
    }
    maxImageIndex = images.length;
    currentImageIndex = 0;
  };

// changes the banner currently being displayed
var changeBanner = function() {
  var i;
  currentImageIndex = (currentImageIndex >= maxImageIndex - 1) ? 0 : currentImageIndex += 1;

  for (i = 0; i < maxImageIndex; i += 1) {
    images[i].hidden = (i !== currentImageIndex);
  }
};

// a function which is triggered following the `load` event
window.onload = function() {
  setUp();
  images[currentImageIndex].hidden = false; // show the first banner image;
  setInterval(changeBanner, changeInterval); // following a delay, keep changing the banner image by the specified interval
};
