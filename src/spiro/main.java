
// Spiro2SVG - Nov. 2015 - Alvin Penner - penner@vaxxine.com

// source : C:\APP\Java\CoreJava\v2ch02\DOMTreeTest\DOMTreeTest.java
// source : C:\Program Files\Java\jdk1.6.0_16\demo\jfc\TableExample\src\TableExample4.java
// source : C:\APP\Java\jswing2\ch15\ColorTable.java

// Nov 28, 2015 - first build - NetBeans 6.9.1 - jdk1.6.0_16
// Feb 14, 2016 - completed lines and points output to svg

package spiro;

import java.io.*;
import java.util.Properties;
import javax.swing.JOptionPane;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

public class main
{
    public static final String VERSION_NO = "0.9";
    public static final String PAGE_UNITS = "mm";
    public static final float PAGE_WIDTH = 210;
    public static final float PAGE_HEIGHT = 297;
    public static final float mm2px = 96F/25.4F;
    public static final float MITER_LIMIT = 20;
    public static final String STYLE_AUTO = "Auto";
    public static final String STYLE_POINTS = "Points";
    public static final String STYLE_LINES = "Lines";
    public static final String STYLE_BEZIER = "Bezier";
    public static final String[][] rowNames = new String[][] {{"StatorRadius"}, {"RotorRadius"}, {"NumRotations"}, {"AnglesPerCycle"},
                                                              {"RotorSlide"}, {"OriginX"}, {"OriginY"}, {"InitialAngle"},
                                                              {"PenDistance"}, {"Lock"}, {"CurvePenWidth"}, {"Zoom"},
                                                              {"Argb"}, {"FillArgb"}, {"FillMode"}, {"Edit Drawing Style"}};
    public static String[][] rowData;
    public static int CanvasColor = 0;
    private static String draw_style = STYLE_AUTO;

    public static final String strHdr = "<?xml version='1.0' encoding='UTF-8' standalone='no'?>\n"
                          + "<!-- Created with Spiro2SVG v" + VERSION_NO + " (http://launchpad.net/tobedetermined) -->\n\n"
                          + "<svg\n"
                          + "    xmlns:dc='http://purl.org/dc/elements/1.1/'\n"
                          + "    xmlns:cc='http://creativecommons.org/ns#'\n"
                          + "    xmlns:rdf='http://www.w3.org/1999/02/22-rdf-syntax-ns#'\n"
                          + "    xmlns:svg='http://www.w3.org/2000/svg'\n"
                          + "    xmlns:inkscape='http://www.inkscape.org/namespaces/inkscape'\n"
                          + "    xmlns='http://www.w3.org/2000/svg'\n"
                          + "    version='1.1'\n"
                          + "    height='" + PAGE_HEIGHT + PAGE_UNITS + "'\n"
                          + "    width='" + PAGE_WIDTH + PAGE_UNITS + "'>\n"
                          + "    <metadata>\n"
                          + "        <rdf:RDF>\n"
                          + "            <cc:Work rdf:about=''>\n"
                          + "                <dc:format>image/svg+xml</dc:format>\n"
                          + "                <dc:type rdf:resource='http://purl.org/dc/dcmitype/StillImage' />\n"
                          + "            </cc:Work>\n"
                          + "        </rdf:RDF>\n"
                          + "    </metadata>\n";
    public static final String strFtr = "</svg>\n";
    public static double t1 = -1, t2 = -1;                      // test parms for testing bezier fit

//  The parameters fall into four roughly defined classes

//  Class 1 - style parameters
//  overall Argb = CanvasColor : background color only
//  12 - Argb = stroke color and opacity
//  13 - FillArgb = fill color and opacity
//  10 - CurvePenWidth = stroke width in px (not scaled by Zoom)
//  14 - FillMode = evenodd/nonzero
//  15 - Edit Drawing Style = user-defined: Points/Lines/Bezier

//  Class 2 - scaling applied after path is created
//  05 - OriginX - in px relative to page center
//  06 - OriginY - in px relative to page center
//  07 - InitialAngle - clockwise rotation of object: ~500*radians
//  11 - Zoom - object size in %

//  Class 3 - Spirograph 1.0.2.1 rendering parameters
//  see : http://mathiversity.com/online-spirograph
//  or  : http://sourceforge.net/projects/spirograph
//  00 - StatorRadius - static radius in px
//  01 - RotorRadius  - rolling radius in px
//  02 - NumRotations - integer number of full rotations
//  03 - AnglesPerCycle - number of line segments per rotation
//  04 - RotorSlide - modifies Stator/Rotor
//  08 - PenDistance - from center of rolling object in px
//  09 - Lock - forces cycloid behavior

//  Class 4 - standard epi/hypoTrochoid rendering parameters
//  see : http://turnbull.mcs.st-and.ac.uk/~history/Curves/Epitrochoid.html
//  a = stator, b = rotor, c = pen distance
//  x = (a + b) cos(t) - c cos((a/b + 1)t)
//  y = (a + b) sin(t) - c sin((a/b + 1)t)

    public static void main(String[] args)
    {
        String fileName = "";
        String fileNameOnly = "";
        String exportName = "";
        Properties pgmProp = new Properties();

        if (args.length > 0)                                // run from command line
        {
            if (args[0].equals("-?") || args[0].equals("-h") || args[0].equals("--help"))
            {
                System.out.println("\nSpiro2SVG v" + main.VERSION_NO);
                System.out.println("\nUsage: java -jar Spiro2SVG.jar [OPTIONS] [FILEIN (.spiro)] [FILEOUT (.svg)]");
                System.out.println("\nAvailable Options: -?, -h, --help       This message");
                System.out.println("                 : -p                   Points only");
                System.out.println("                 : -l                   Lines only");
                System.out.println("                 : -b                   Bezier curves only");
                System.out.println("default rendering is 'Lines' if less than 20 points/lobe, 'Bezier' otherwise.");
                return;
            }
            if (args[0].equals("-p"))
                draw_style = STYLE_POINTS;
            else if(args[0].equals("-l"))
                draw_style = STYLE_LINES;
            else if(args[0].equals("-b"))
                draw_style = STYLE_BEZIER;
            int iIn = -1;                                   // index of input file name
            if (!args[0].startsWith("-"))
                iIn = 0;
            else if ((args.length > 1) && !args[1].startsWith("-"))
                iIn = 1;
            if (iIn > -1)
            {
                fileName = args[iIn];
                if ((args.length > iIn + 1) && !args[iIn + 1].startsWith("-"))
                    exportName = args[iIn + 1];
                if ((args.length > iIn + 3) && !args[iIn + 1].startsWith("-"))
                {
                    t1 = Double.parseDouble(args[iIn + 2]);   // for testing only
                    t2 = Double.parseDouble(args[iIn + 3]);
                    t1 *= 2*Math.PI;
                    t2 *= 2*Math.PI;
                }
            }
        }
        else                                                // run with gui
        {
            try                                             // recall program properties
            {
                if (new File(System.getProperty("user.home"), "Spiro2SVGPrefs.ini").exists())
                    pgmProp.load(new FileInputStream(new File(System.getProperty("user.home"), "Spiro2SVGPrefs.ini")));
                else
                    {HelpAbout.showDialog(null);}
            }
            catch (IOException e)
                {System.out.println("error reading Spiro2SVGPrefs.ini : " + e);}

            JFileChooser chooser = new JFileChooser(pgmProp.getProperty("pathname", System.getProperty("user.home")));
            chooser.setFileFilter(new FileNameExtensionFilter("spiro files (*.spiro)", "spiro"));
            if (chooser.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) return;
            pgmProp.setProperty("pathname", chooser.getSelectedFile().getAbsoluteFile().getParent());
            fileName = chooser.getSelectedFile().getAbsolutePath();
            if (chooser.getFileFilter().getDescription().equals("spiro files")  && !fileName.endsWith(".spiro"))
                fileName += ".spiro";
            fileNameOnly = chooser.getSelectedFile().getName();
            if (fileNameOnly.endsWith(".spiro"))
                fileNameOnly = fileNameOnly.substring(0, fileNameOnly.length() - 6);
            fileNameOnly += ".svg";
        }

        SpiroParse.parse_spiro_file(fileName, draw_style);
        if (rowData == null) return;

        if (args.length == 0)                                           // run with gui
        {
            if (!SpiroConfig.showDialog(null, fileName)) return;
            JFileChooser chooser = new JFileChooser(pgmProp.getProperty("exportpath", pgmProp.getProperty("pathname", System.getProperty("user.home"))));
            chooser.setSelectedFile(new File(fileNameOnly));
            chooser.setFileFilter(new FileNameExtensionFilter("Scalable Vector Graphics (*.svg)", "svg"));
            if (chooser.showSaveDialog(null) != JFileChooser.APPROVE_OPTION) return;
            pgmProp.setProperty("exportpath", chooser.getSelectedFile().getAbsoluteFile().getParent());
            exportName = chooser.getSelectedFile().getAbsolutePath();
            if (!exportName.endsWith(".svg"))
                exportName += ".svg";
            try
                {pgmProp.store(new FileOutputStream(System.getProperty("user.home") + System.getProperty("file.separator") + "Spiro2SVGPrefs.ini"), "Spiro2SVG v" + VERSION_NO + " Prefs");}
            catch (IOException e)
                {JOptionPane.showMessageDialog(null, e.getMessage(), " Could not save Preferences file", JOptionPane.WARNING_MESSAGE);}
        }
        SpiroWrite.write_svg_file(exportName);

        System.out.println("\n" + fileName + " : rowData = " + rowData.length + ", " + rowData[0].length);
        for (int i = 0; i < rowData.length; i++)
        {
            System.out.print(rowNames[i][0] + "\t");
            for (int j = 0; j < rowData[0].length; j++)
                System.out.print(", " + rowData[i][j]);
            System.out.println();
        }
    }
}
