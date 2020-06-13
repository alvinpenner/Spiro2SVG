
package zeeman;

// see : C:\APP\Java\CoreJava\v2ch06\EditorPaneTest\EditorPaneTest.java
// see : C:\APP\Java\CoreJava\v2ch07\DesktopAppTest\DesktopAppTest.java

import java.awt.Desktop;
import java.awt.Dimension;
import java.net.*;
import javax.swing.*;
import javax.swing.event.*;

public final class ZeemanAbout
{
    public static void showDialog(JFrame parent)
    {
        JDialog help = new JDialog(parent, " Zeeman Catastrophe Machine v" + main.VERSION_NO + " - by Alvin Penner", false);
        JEditorPane editorPane = new JEditorPane();
        editorPane.setEditable(false);
        URL helpURL = ZeemanAbout.class.getResource("images/aboutZCM.html");
        try {
            editorPane.setPage(helpURL);
        } catch (java.io.IOException e) {
            System.err.println("Attempted to read a bad URL: " + helpURL);
        }
        editorPane.addHyperlinkListener(new HyperlinkListener() {
            public void hyperlinkUpdate(HyperlinkEvent event)
            {
                if (event.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
                {
                    if (event.getURL().getProtocol().equals("http"))
                    {
                        try {Desktop.getDesktop().browse(new URI(event.getURL().toString()));}
                        catch (URISyntaxException ex) {System.err.println(ex);}
                        catch (java.io.IOException ex) {System.err.println(ex);}
                    }
                    else if (event.getURL().getProtocol().equals("mailto"))
                    {
                        try {Desktop.getDesktop().mail(new URI(event.getURL().toString()));}
                        catch (URISyntaxException ex) {System.err.println(ex);}
                        catch (java.io.IOException ex) {System.err.println(ex);}
                    }
                }
            }
        });

        help.add(editorPane);
        help.setMinimumSize(new Dimension(600, 600));
        help.setResizable(false);
        help.setLocationByPlatform(true);
        help.setVisible(true);
    }
}
