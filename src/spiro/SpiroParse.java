
package spiro;

import java.io.File;
import javax.swing.JOptionPane;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.*;

//  this will parse two types of spirograph files:
//  .spiro files produced by Spirograph_1.0.2.1 - http://mathiversity.com/online-spirograph
//  .xml   files produced by SpiroJ_1.0.2       - http://sourceforge.net/projects/spiroj

public final class SpiroParse
{
    private static int spiro_count = 0;
    private static String draw_style;

    public static void parse_spiro_file(String fname, String m_style)
    {
        final File file = new File(fname);
        if (!file.exists())
        {
            JOptionPane.showMessageDialog(null, "The file '" + file.getAbsolutePath() + "' does not exist.\nPlease try again." , " File not found ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        main.rowNames = main.spiroNames;
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

    public static void parse_SpiroJ_file(String fname, String m_style)
    {
        final File file = new File(fname);
        if (!file.exists())
        {
            JOptionPane.showMessageDialog(null, "The file '" + file.getAbsolutePath() + "' does not exist.\nPlease try again." , " File not found ", JOptionPane.WARNING_MESSAGE);
            return;
        }

        // SpiroJ file  has 2 rotors, 4 parameters per rotor (2 size, 2 frequency)
        // Farris Wheel has 3 rotors, 3 parameters per rotor (size, frequency, phase)
        // both types use the same tags: shape, generator, operator

        draw_style = m_style;
        DocumentBuilder builder;
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        try
        {
            builder = factory.newDocumentBuilder();
            Document doc = builder.parse(file);
            Element root = doc.getDocumentElement();
            if (!root.getTagName().equals("spiro"))
            {
                JOptionPane.showMessageDialog(null, "This does not appear to be a valid SpiroJ file.\nElement 'spiro' not found.", " parse_SpiroJ error ", JOptionPane.WARNING_MESSAGE);
                return;
            }
            int nrotor = count_Children(root, "operator");
            if (nrotor == 2)                                            // SpiroJ format
                main.rowNames = main.SpiroJNames;
            else if (nrotor == 3)                                       // my version of Farris Wheel format
                main.rowNames = main.FarrisNames;
            else
            {
                JOptionPane.showMessageDialog(null, "This does not appear to be a valid SpiroJ file or a Farris Wheel file.\nNeed two or three 'operator' elements.", " parse_SpiroJ error ", JOptionPane.WARNING_MESSAGE);
                return;
            }
            main.rowData = new String[main.rowNames.length][1];
            NodeList children = root.getChildNodes();
            int irotor = 0;
            for (int i = 0; i < children.getLength(); i++)
            {
                Node child = children.item(i);
                if (child instanceof Element)
                {
                    Element childElement = (Element) child;
                    if (childElement.getTagName().equals("shape"))
                    {
                        main.rowData[7 + nrotor][0]  = childElement.getAttribute("linewidth");  // line_width
                        main.rowData[8 + nrotor][0] = childElement.getAttribute("linecolor");   // line_color
                        main.rowData[9 + nrotor][0] = childElement.getAttribute("fillcolor");   // fill_color
                    }
                    else if (childElement.getTagName().equals("generator"))
                    {
                        main.rowData[6 + nrotor][0] = childElement.getAttribute("steps");   // generator_steps
                    }
                    else if (childElement.getTagName().equals("operator"))
                    {
                        NodeList opchildren = childElement.getChildNodes();
                        for (int j = 0; j < opchildren.getLength(); j++)
                        {
                            Node opchild = opchildren.item(j);
                            if (opchild instanceof Element)
                            {
                                Element opchildElement = (Element) opchild;
                                if (opchildElement.getTagName().equals("radius"))
                                {
                                    if (nrotor == 2)
                                    {
                                        main.rowData[4*irotor][0] = opchildElement.getAttribute("x");       // Radius_x[i]
                                        main.rowData[4*irotor + 1][0] = opchildElement.getAttribute("y");   // Radius_y[i]
                                    }
                                    else if (nrotor == 3)
                                        main.rowData[3*irotor][0] = opchildElement.getAttribute("r");       // Radius[i]
                                }
                                else if (opchildElement.getTagName().equals("frequency"))
                                {
                                    if (nrotor == 2)
                                    {
                                        main.rowData[4*irotor + 2][0] = opchildElement.getAttribute("x");   // Frequency_x[i]
                                        main.rowData[4*irotor + 3][0] = opchildElement.getAttribute("y");   // Frequency_y[i]
                                    }
                                    else if (nrotor == 3)
                                        main.rowData[3*irotor + 1][0] = opchildElement.getAttribute("w");   // Frequency[i]
                                }
                                else if (opchildElement.getTagName().equals("phase"))
                                    if (nrotor == 3)
                                        main.rowData[3*irotor + 2][0] = opchildElement.getAttribute("phi"); // Phase[i]
                            }
                        }
                        irotor++;
                    }
                }
            }
            main.rowData[main.rowNames.length - 1][0] = draw_style;     // Edit Drawing Style
            if (draw_style.equals(main.STYLE_AUTO))
                main.rowData[main.rowNames.length - 1][0] = main.STYLE_BEZIER;
        }
        catch (Exception e) {JOptionPane.showMessageDialog(null, "DocumentBuilderFactory : " + e, " parse_SpiroJ error ", JOptionPane.WARNING_MESSAGE);}
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
