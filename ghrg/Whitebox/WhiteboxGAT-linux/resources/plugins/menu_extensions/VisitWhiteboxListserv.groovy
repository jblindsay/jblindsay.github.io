import java.net.URI

def parentMenu = "Help"
def menuLabel = "Visit the Whitebox Listserv"

java.awt.Desktop.getDesktop().browse(new URI("http://www.lsoft.com/scripts/wl.exe?SL1=WHITEBOX-GAT&H=LISTSERV.UOGUELPH.CA"))
