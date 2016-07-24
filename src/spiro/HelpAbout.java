
package spiro;

// see : C:\APP\Java\CoreJava\v2ch06\EditorPaneTest\EditorPaneTest.java
// see : C:\APP\Java\CoreJava\v2ch07\DesktopAppTest\DesktopAppTest.java

import java.awt.Desktop;
import java.awt.Dimension;
import java.net.*;
import javax.swing.*;
import javax.swing.event.*;

public final class HelpAbout
{
    public static void showDialog(JDialog parent)
    {
        JDialog help = new JDialog(parent, " Spiro2SVG v" + main.VERSION_NO + " - by Alvin Penner", true);
        JEditorPane editorPane = new JEditorPane();
        editorPane.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 10));
        editorPane.setBackground(UIManager.getColor("TabbedPane.contentAreaColor"));
        editorPane.setEditable(false);
        URL helpURL = HelpAbout.class.getResource("images/about.html");
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

        help.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        help.add(editorPane);
        help.setMinimumSize(new Dimension(500, 620));
        help.setResizable(false);
        help.setLocationByPlatform(true);
        help.setVisible(true);
    }
}
