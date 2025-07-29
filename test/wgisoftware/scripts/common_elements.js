// var OSName = "Unknown";
// if (window.navigator.userAgent.indexOf("Windows") != -1) OSName="Windows";
// if (window.navigator.userAgent.indexOf("Mac")!=-1) OSName="Mac/iOS";
// if (window.navigator.userAgent.indexOf("X11")!=-1) OSName="UNIX";
// if (window.navigator.userAgent.indexOf("Linux")!=-1) OSName="Linux";
// var eol = "\n";
// if (OSName === "Windows") {
//   eol = "\r\n";
// }

function insertNavbar() {
  var el = document.getElementById("navbar");
  if (el != null) {
    el.innerHTML = "";
    s = `<div class="navbar" id="myNavbar">
            <div class="company-info">
                <img src="img/logo-muted.svg" alt="FPC logo" style="width: auto; height: 35px;"></img>
                <span>&nbsp;</span>
                <a class="company-name" href="./index.html">Flash Point Classifier</a>
            </div>
            <div class="navbar-links">
                <a href="./download.html">Download</a>
                <a href="./purchase.html">Purchase</a>
                <a href="./support.html">Support</a>
                <a href="./about.html">About Us</a>
                <div><span>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</span></div>
            </div>
        </div>`;
    el.innerHTML += s;
  }
}

function insertFooter() {
  var el = document.getElementById("footer");
  if (el != null) {
    el.innerHTML = "";
    s = `<div class="footer">
            <div style="display: flex; flex-wrap: wrap; justify-content: space-around; gap: 20px;">
                <div>
                    <p style="font-weight: bold;">Quick Links</p>
                    <div style="display: flex; flex-wrap: wrap; gap: 10px;">
                        <div>
                            <p><a href="index.html">Home</a></p>
                            <p><a href="download.html">Download</a></p>
                            <p><a href="purchase.html">Purchase</a></p>
                            <p><a href="support.html">Support</a></p>
                        </div>
                        <div>
                            <p><a href="about.html">About</a></p>
                            <p><a href="license.html">Terms & Conditions</a></p>
                            <p><a href="license.html#refund-policy">Refund Policy</a></p>
                        </div>
                    </div>
                </div>
                <div>
                    <p style="font-weight: bold;">Contact Us</p>
                    <p><img src="img/icons/location.svg" alt="Location Icon" style="width:18px;height:auto;vertical-align: middle;"> WG Inc. Guelph, Canada</p>
                    <p><img src="img/icons/email.svg" alt="Email Icon" style="width:16px;height:auto;vertical-align: middle;"><a href="mailto:contact@wgisoftware.com"> contact@wgisoftware.com</a></p>
                </div>
                <div style="justify-content: center; align-items: center;">
                    <img src="img/logo-muted.svg" alt="Logo" style="width:48px;height:auto;vertical-align:middle;padding-bottom:5px;">
                    <p>&copy; 2025 WG Inc. All rights reserved.</p>
                </div>
            </div>
        </div>`;
    el.innerHTML += s;
  }
}

function redirectToPage(page) {
    window.location.href = `./${page}.html`;
}