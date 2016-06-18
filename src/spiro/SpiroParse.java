
package spiro;

import java.io.File;
import javax.swing.JOptionPane;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.*;

public final class SpiroParse
{
    private static int spiro_count = 0;
    private static String draw_style;

    public static void parse_spiro_file(String fname, String m_style)
    {
        if (fname.isEmpty())
        {
            JOptionPane.showMessageDialog(null, "No input file name was specified.\nPlease try again." , " No input location ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        final File file = new File(fname);
        if (!file.exists())
        {
            JOptionPane.showMessageDialog(null, "The file '" + file.getAbsolutePath() + "' does not exist.\nPlease try again." , " File not found ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        draw_style = m_style;
        DocumentBuilder builder;
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        try
        {
            builder = factory.newDocumentBuilder();
            Document doc = builder.parse(file);
            Element root = doc.getDocumentElement();
            if (!root.getTagName().equals("Spiro"))
            {
                JOptionPane.showMessageDialog(null, "This does not appear to be a valid '.spiro' file.\nElement 'Spiro' not found.", " parse_spiro error ", JOptionPane.WARNING_MESSAGE);
                return;
            }
            NodeList children = root.getChildNodes();
            for (int i = 0; i < children.getLength(); i++)
            {
                Node child = children.item(i);
                if (child instanceof Element)
                {
                    Element childElement = (Element) child;
                    if (childElement.getTagName().equals("Argb"))
                        main.CanvasColor = Integer.parseInt(childElement.getFirstChild().getNodeValue());
                    else if (childElement.getTagName().equals("SpiroCurves"))
                    {
                        main.rowData = new String[main.rowNames.length][count_Children(childElement, "Curve")];
                        for (int ii = 0; ii < main.rowData[0].length; ii++)
                        {
                            main.rowData[11][ii] = "100";              // default Zoom
                            main.rowData[13][ii] = "0";                // default FillArgb
                            main.rowData[14][ii] = "0";                // default FillMode
                        }
                        parse_SpiroCurves(childElement);
                    }
                }
            }
        }
        catch (Exception e) {JOptionPane.showMessageDialog(null, "DocumentBuilderFactory : " + e, " parse_spiro error ", JOptionPane.WARNING_MESSAGE);}
    }

    private static int count_Children(Element element, String tag)
    {
        int count = 0;
        NodeList children = element.getChildNodes();
        for (int i = 0; i < children.getLength(); i++)
        {
            Node child = children.item(i);
            if (child instanceof Element)
            {
                Element childElement = (Element) child;
                if (childElement.getTagName().equals(tag))
                    count++;
            }
        }
        return count;
    }

    private static void parse_SpiroCurves(Element element)
    {
        NodeList children = element.getChildNodes();
        for (int i = 0; i < children.getLength(); i++)
        {
            Node child = children.item(i);
            if (child instanceof Element)
            {
                Element childElement = (Element) child;
                if (childElement.getTagName().equals("Curve"))
                    parse_SpiroCurves(childElement);
                else if (childElement.getTagName().equals("CurveDO"))
                    parse_CurveDO(childElement);
            }
        }
    }

    private static void parse_CurveDO(Element element)
    {
        NodeList children = element.getChildNodes();
        for (int i = 0; i < children.getLength(); i++)
        {
            Node child = children.item(i);
            if (child instanceof Element)
            {
                Element childElement = (Element) child;
                for (int ii = 0; ii < main.rowNames.length; ii++)
                    if (main.rowNames[ii][0].equals(childElement.getTagName()))
                        main.rowData[ii][spiro_count] = childElement.getFirstChild().getNodeValue().trim();
            }
        }
        main.rowData[main.rowNames.length - 1][spiro_count] = draw_style;   // Edit Drawing Style
        if (draw_style.equals(main.STYLE_AUTO))
            if (Float.parseFloat(main.rowData[3][spiro_count])*Float.parseFloat(main.rowData[1][spiro_count])
            >   20*Float.parseFloat(main.rowData[0][spiro_count]))
                main.rowData[main.rowNames.length - 1][spiro_count] = main.STYLE_BEZIER;
            else
                main.rowData[main.rowNames.length - 1][spiro_count] = main.STYLE_LINES;
        spiro_count++;
    }
}
