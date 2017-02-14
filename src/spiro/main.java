
// Spiro2SVG - Nov. 2015 - Alvin Penner - penner@vaxxine.com

// source : C:\APP\Java\CoreJava\v2ch02\DOMTreeTest\DOMTreeTest.java
// source : C:\Program Files\Java\jdk1.6.0_16\demo\jfc\TableExample\src\TableExample4.java
// source : C:\APP\Java\jswing2\ch15\ColorTable.java

// Nov 28, 2015 - first build - NetBeans 6.9.1 - jdk1.6.0_16
// Feb 14, 2016 - completed lines and points output to svg
// June  , 2016 - completed Bezier output
// Jun 18, 2016 - first commit to http://github.com/alvinpenner/Spiro2SVG
// Jul 11, 2016 - rev 0.91, support for Lissajous figures from SpiroJ
// Jul 24, 2016 - rev 0.92, support for all figures from SpiroJ
// Aug  8, 2016 - rev 0.93, support for Farris Wheels, with no inflection points
// Aug 10, 2016 - fix Farris Wheel crash caused by multiple simultaneous roots
// Aug 11, 2016 - support for inflection points in Farris Wheels
// Nov  9, 2016 - for Farris Wheel, use extrema of curvature, plus perpendicular slope
// Nov 25, 2016 - rev 0.94, support for cusps in Farris Wheels
// Dec  5, 2016 - for standard spirograph, Class 4, switch from -c to c

package spiro;

import java.awt.geom.Point2D;
import java.awt.geom.CubicCurve2D;
import java.io.*;
import java.util.Properties;
import javax.swing.JOptionPane;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

public class main
{
    public static final String VERSION_NO = "0.94";
    public static final String PAGE_UNITS = "mm";
    public static final float PAGE_WIDTH = 210;                 // A4 page in mm
    public static final float PAGE_HEIGHT = 297;
//    public static final float PAGE_WIDTH = 150*25.4F/96;        // 150 px
//    public static final float PAGE_HEIGHT = 170*25.4F/96;
    public static final float mm2px = 96F/25.4F;                // SVG resolution is 96 dpi
    public static final float MITER_LIMIT = 20;
    public static final String STYLE_AUTO = "Auto";
    public static final String STYLE_POINTS = "Points";
    public static final String STYLE_LINES = "Lines";
    public static final String STYLE_BEZIER = "Bezier";
    public static final boolean FIT_POINTS_ONLY = false;         // for debugging only
    public static final boolean IS_DEBUG = true;                // for debugging only
    public static final String[][] spiroNames = new String[][] {{"StatorRadius"}, {"RotorRadius"}, {"NumRotations"}, {"AnglesPerCycle"},
                                                                {"RotorSlide"}, {"OriginX"}, {"OriginY"}, {"InitialAngle"},
                                                                {"PenDistance"}, {"Lock"}, {"CurvePenWidth"}, {"Zoom"},
                                                                {"Argb"}, {"FillArgb"}, {"FillMode"}, {"Edit Drawing Style"}};
    public static final String[][] SpiroJNames = new String[][] {{"Radius_x1"}, {"Radius_y1"}, {"Frequency_x1"}, {"Frequency_y1"},
                                                                 {"Radius_x2"}, {"Radius_y2"}, {"Frequency_x2"}, {"Frequency_y2"},
                                                                 {"generator_steps"}, {"line_width"}, {"line_color"}, {"fill_color"},
                                                                 {"Edit Drawing Style"}};
    public static final String[][] FarrisNames = new String[][] {{"Radius_1"}, {"Frequency_1"}, {"Phase_1"}, {"Radius_2"},
                                                                 {"Frequency_2"}, {"Phase_2"}, {"Radius_3"}, {"Frequency_3"}, {"Phase_3"},
                                                                 {"generator_steps"}, {"line_width"}, {"line_color"}, {"fill_color"},
                                                                 {"Edit Drawing Style"}};
    public static String[][] rowNames;
    public static String[][] rowData;
    public static int CanvasColor = 0;
    private static String draw_style = STYLE_AUTO;

    public static final String strHdr = "<?xml version='1.0' encoding='UTF-8' standalone='no'?>\n"
                          + "<!-- Created with Spiro2SVG v" + VERSION_NO + " (http://github.com/alvinpenner/Spiro2SVG) -->\n\n"
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
    public static double t1 = -1, t2 = -1;                              // test parms for testing bezier fit

    public static double[] theta = new double[2];                       // Spiro velocity angle[t = (0,1)]
    public static double[] Cu = new double[2];                          // Spiro curvature[t = (0,1)]
    private static final double TOL = 0.00001;
    private static Point2D.Double[] ptBez = new Point2D.Double[4];      // Point[point index (0-3)]
    private static double xrot0, yrot0, xrot3, yrot3, mrot0;

//  The parameters fall into five roughly defined classes
//  (arbitrarily made up classification)

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
//  x = (a + b) cos(t) + c cos((a/b + 1)t)
//  y = (a + b) sin(t) + c sin((a/b + 1)t)

//  Class 5 - SpiroJ roulette parameters
//  see : http://sourceforge.net/projects/spiroj/
//  x = rx1*cos(wx1*t) + rx2*cos(wx2*t)
//  y = ry1*sin(wy1*t) + ry2*sin(wy2*t)

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
            chooser.setFileFilter(new FileNameExtensionFilter("spirograph files (.spiro) (.xml)", "spiro", "xml"));
            if (chooser.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) return;
            pgmProp.setProperty("pathname", chooser.getSelectedFile().getAbsoluteFile().getParent());
            fileName = chooser.getSelectedFile().getAbsolutePath();
            if (chooser.getFileFilter().getDescription().equals("spirograph files (.spiro) (.xml)")
            && !fileName.endsWith(".spiro") && !fileName.endsWith(".xml"))
                if (new File(fileName + ".spiro").exists())
                    fileName += ".spiro";
                else if (new File(fileName + ".xml").exists())
                    fileName += ".xml";
            fileNameOnly = chooser.getSelectedFile().getName();
            if (fileNameOnly.endsWith(".spiro"))
                fileNameOnly = fileNameOnly.substring(0, fileNameOnly.length() - 6);
            else if (fileNameOnly.endsWith(".xml"))
                fileNameOnly = fileNameOnly.substring(0, fileNameOnly.length() - 4);
            fileNameOnly += ".svg";
        }

        if (fileName.isEmpty())
        {
            JOptionPane.showMessageDialog(null, "No input file name was specified.\nPlease try again." , " No input location ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (fileName.endsWith(".spiro"))
            SpiroParse.parse_spiro_file(fileName, draw_style);
        else if (fileName.endsWith(".xml"))
            SpiroParse.parse_SpiroJ_file(fileName, draw_style);
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
        SpiroWrite.write_svg_file(exportName, fileName.endsWith(".spiro"));

        if (IS_DEBUG)
        {
            System.out.println("\n" + fileName + " : ");
            for (int i = 0; i < rowData.length; i++)
            {
                System.out.print(rowNames[i][0] + "\t");
                for (int j = 0; j < rowData[0].length; j++)
                    System.out.print(", " + rowData[i][j]);
                System.out.println();
            }
        }
    }

    protected static CubicCurve2D.Float calcBezier(Point2D.Double[][] ptSpiro, double t1, double t2, double max_v)
    {
        double delxrot0 = 0, delyrot0, delxrot3 = 0, delyrot3;
        double[] dirx = new double[2];
        double[] scaled_v = new double[2];                                  // check for stationary point

        ptBez[0] = new Point2D.Double(ptSpiro[0][0].x, ptSpiro[0][0].y);
        ptBez[1] = new Point2D.Double(ptBez[0].x, ptBez[0].y);
        ptBez[3] = new Point2D.Double(ptSpiro[0][1].x, ptSpiro[0][1].y);
        ptBez[2] = new Point2D.Double(ptBez[3].x, ptBez[3].y);

        xrot0 = getrotX(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        yrot0 = getrotY(ptBez[0].x, ptBez[0].y, (theta[0] + theta[1])/2);
        xrot3 = getrotX(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        yrot3 = getrotY(ptBez[3].x, ptBez[3].y, (theta[0] + theta[1])/2);
        mrot0 = Math.tan((theta[0] - theta[1])/2);

        for (int i = 0; i < ptSpiro[0].length; i++)                         // effective 'rotated' curvature
        {
            scaled_v[i] = Math.sqrt(ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y)/max_v;
            if (scaled_v[i] < TOL)                                          // stationary point
            {
                if (IS_DEBUG)
                    System.out.println("calcBezier stationary point at i = " + i + " : " + t1 + ", " + t2 + ", " + Math.sqrt(ptSpiro[1][i].x*ptSpiro[1][i].x + ptSpiro[1][i].y*ptSpiro[1][i].y));
                dirx[i] = getrotX(ptSpiro[2][i].x, ptSpiro[2][i].y, (theta[0] + theta[1])/2);   // use x″ & y″
            }
            else                                                            // moving point
                dirx[i] = getrotX(ptSpiro[1][i].x, ptSpiro[1][i].y, (theta[0] + theta[1])/2);
            Cu[i] = Math.signum(dirx[i])*Cu[i]*Math.pow(1 + mrot0*mrot0, 1.5);
//            System.out.println("slopes " + i + " : " + ptSpiro[1][i].x + ", " + ptSpiro[1][i].y + " : " + ptSpiro[2][i].x + ", " + ptSpiro[2][i].y);
//            System.out.println(i + " : " + Cu[i] + ", " + dirx[i] + ", " + ptSpiro[2][i].x + ", " + ptSpiro[2][i].y);
        }
        for (int i = 0; i < ptSpiro[0].length; i++)
            if (scaled_v[i] < TOL)                                          // stationary point
            {
                if (Math.signum(dirx[i]) != Math.signum(dirx[1 - i]))
                    dirx[i] *= -1;                  // force same sign of dirx
                if (Math.signum(Cu[i]) != Math.signum(Cu[1 - i]))
                    Cu[i] *= -1;                    // force same sign of Cu
            }

        if (IS_DEBUG)
        {
//          System.out.println("angles = " + Math.atan2(ptSpiro[1][0].y, ptSpiro[1][0].x)*180/Math.PI + ", " + Math.atan2(ptSpiro[1][1].y, ptSpiro[1][1].x)*180/Math.PI + ", " +  mrot0);
            System.out.println("angles = " + (float)t1 + ", " + (float)ptSpiro[0][0].x + ", " + (float)ptSpiro[0][0].y + ", " + (float)(theta[0]*180/Math.PI) + ", " + (float)Cu[0]);
            System.out.println("angles = " + (float)t2 + ", " + (float)ptSpiro[0][1].x + ", " + (float)ptSpiro[0][1].y + ", " + (float)(theta[1]*180/Math.PI) + ", " + (float)Cu[1]);
            System.out.println("data   = " + (float)xrot0 + ", " + (float)xrot3 + ", " + (float)yrot0 + ", " + (float)yrot3 + ", " + (float)mrot0 + ", " + (float)dirx[0] + ", " + (float)dirx[1]);
        }

        if (Math.abs(theta[1] - theta[0]) < TOL)                // parallel finite slopes
        {
            // the parallel case can be solved as two decoupled quadratic equations
            if (IS_DEBUG)
                System.out.println("parallel finite : " + (float)(t1/2/Math.PI) + ", " + (float)(t2/2/Math.PI) + ", " + (float)mrot0 + ", " + (float)Cu[0] + ", " + (float)Cu[1] + ", " + (float)yrot0 + ", " + (float)yrot3);
            if (Math.abs(yrot3 - yrot0) > TOL)
                delxrot0 =  2*(yrot3 - yrot0)/Cu[0]/3;
            if (delxrot0 < 0)
            {
                System.out.println("parallel slope, wrong curvature at t = 0 : " + delxrot0 + " : abort");
                delxrot0 = Double.NaN;
                return new CubicCurve2D.Float();
            }
            else
                delxrot0 = Math.signum(dirx[0]*(t2 - t1))*Math.sqrt(delxrot0);
            if (Math.abs(yrot3 - yrot0) > TOL)
                delxrot3 = -2*(yrot3 - yrot0)/Cu[1]/3;
            if (delxrot3 < 0)
            {
                System.out.println("parallel slope, wrong curvature at t = 1 : " + delxrot3 + " : abort");
                delxrot3 = Double.NaN;
                return new CubicCurve2D.Float();
            }
            else
                delxrot3 = Math.signum(dirx[1]*(t2 - t1))*Math.sqrt(delxrot3);
        }
        else if ((Math.abs(Cu[0]) < TOL)                        // zero curvature at t = 0
             &&  (scaled_v[1] > TOL))                           // not stationary at t = 1
        {
            if (IS_DEBUG)
                System.out.println("zero curvature at t = 0 : " + Math.signum(ptSpiro[1][0].y) + ", " + Math.signum(t2 - t1) + ", " + Cu[0]);
            delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2/mrot0;
            if (scaled_v[0] < TOL)
                delxrot0 = 0;                                   // stationary point, clamp it
            else if (Math.signum(yrot3 - yrot0 + mrot0*(xrot3 - xrot0)) != Math.signum(yrot3 - yrot0 + mrot0*(xrot3 - xrot0) + 3*Cu[1]*delxrot3*delxrot3/2))
                delxrot0 = 0;                                   // unwanted sign reversal
            else
                delxrot0 =  (yrot3 - yrot0 + mrot0*(xrot3 - xrot0) + 3*Cu[1]*delxrot3*delxrot3/2)/2/mrot0;
            if (IS_DEBUG)
                System.out.println("delxrot0/3 = " + delxrot0 + ", " + delxrot3);
        }
        else if (Math.abs(Cu[1]) < TOL)                         // zero curvature at t = 1
        {
            if (IS_DEBUG)
                System.out.println("zero curvature at t = 1 : " + Math.signum(ptSpiro[1][0].y) + ", " + Math.signum(t2 - t1) + ", " + Cu[1]);
            delxrot0 = (yrot3 - yrot0 + mrot0*(xrot3 - xrot0))/2/mrot0;
            if (scaled_v[1] < TOL)
                delxrot3 = 0;                                   // stationary point, clamp it
            else if (Math.signum(yrot3 - yrot0 - mrot0*(xrot3 - xrot0)) != Math.signum(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delxrot0*delxrot0/2))
                delxrot3 = 0;                                   // unwanted sign reversal
            else
                delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delxrot0*delxrot0/2)/2/mrot0;
            if (IS_DEBUG)
                System.out.println("delxrot0/3 = " + delxrot0 + ", " + delxrot3);
        }
        else                                                    // general slopes
        {
            if (IS_DEBUG)
            {
                System.out.println("general slope = " + (float)t1 + ", " + (float)t2 + ", " + (float)(theta[0]*180/Math.PI) + ", " + (float)(theta[1]*180/Math.PI) + ", " + (float)((theta[0] + theta[1])*90/Math.PI) + ", " + (float)mrot0 + ", " + (float)Cu[0] + ", " + (float)Cu[1]);
                System.out.println("general slope = " + (float)t1 + ", " + (float)t2 + ", " + (float)xrot0 + ", " + (float)yrot0 + ", " + (float)xrot3 + ", " + (float)yrot3);
            }
            delxrot0 = solve_quartic(-27*Cu[1]*Cu[0]*Cu[0]/8,
                                      0,
                                      9*Cu[1]*Cu[0]*(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2,
                                      8*mrot0*mrot0*mrot0,
                                     -(yrot3 - yrot0 + mrot0*(xrot3 - xrot0))*4*mrot0*mrot0
                                     -3*Cu[1]*(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))
                                             *(yrot3 - yrot0 - mrot0*(xrot3 - xrot0))/2,
                                      Math.signum(dirx[0]*(t2 - t1)),
                                      Math.signum(dirx[1]*(t2 - t1)));
            delxrot3 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delxrot0*delxrot0/2)/2/mrot0;
        }
        delyrot0 = mrot0*delxrot0;
        delyrot3 = -mrot0*delxrot3;
        ptBez[1].x += getrotX(delxrot0, delyrot0, -(theta[0] + theta[1])/2);
        ptBez[1].y += getrotY(delxrot0, delyrot0, -(theta[0] + theta[1])/2);
        ptBez[2].x -= getrotX(delxrot3, delyrot3, -(theta[0] + theta[1])/2);
        ptBez[2].y -= getrotY(delxrot3, delyrot3, -(theta[0] + theta[1])/2);
/*
        System.out.println("\nptSpiro[3][2]");
        for (int i = 0; i < ptSpiro.length; i++)
            for (int j = 0; j < ptSpiro[0].length; j++)
                System.out.println(i + " " + j + " : " + ptSpiro[i][j]);
        System.out.println("\nptBez[4]");
        for (int i = 0; i < ptBez.length; i++)
            System.out.println(i + " : " + ptBez[i]);
        System.out.println("\ntheta[2]");
        for (int i = 0; i < theta.length; i++)
            System.out.println(i + " : " + theta[i]);
        System.out.println("\nCu[2]");
        for (int i = 0; i < Cu.length; i++)
            System.out.println(i + " : " + Cu[i]);
        System.out.println("\ndelxrotx : " + delxrot0 + ", " + delxrot3);
*/
        return new CubicCurve2D.Float((float)ptBez[0].x, (float)ptBez[0].y, (float)ptBez[1].x, (float)ptBez[1].y, (float)ptBez[2].x, (float)ptBez[2].y, (float)ptBez[3].x, (float)ptBez[3].y);
    }

    protected static int insert_t_value(int N, int index, double[] t, double new_t)
    {
        // push t[index] up by one, insert at location index
        if (N < index) return 0;
        if (N > index)
            for (int i = N; i > index; i--)
                t[i] = t[i - 1];
        t[index] = new_t;
        return N + 1;
    }

    protected static void write_test_cubic(FileWriter out)
    {
        // this will write a test cubic bezier, just for testing purposes

        System.out.printf("\ncubic test parms = %f, %f\n\n", t1, t2);
        SpiroCalc.getBezier(t1, t2);                                    // refresh ptBez[]
        String strPath;
        for (int i = 0; i < ptBez.length; i++)                          // re-define origin
        {
            ptBez[i].x += main.PAGE_WIDTH*main.mm2px/2;
            ptBez[i].y += main.PAGE_WIDTH*main.mm2px/2;
        }
        try
        {
            strPath = "    <path\n"
                    + "        d='";
            out.write(strPath);
            strPath =  "M " + (float)ptBez[0].x + ", " + (float)ptBez[0].y +
                      " C " + (float)ptBez[1].x + ", " + (float)ptBez[1].y + " " + (float)ptBez[2].x + ", " + (float)ptBez[2].y + " " + (float)ptBez[3].x + ", " + (float)ptBez[3].y;
            out.write(strPath);
            strPath = "'\n"
                    + "        id='test_cubic'\n"
                    + "        style='fill:none;stroke:#0000ff;stroke-width:1px' />\n";
            out.write(strPath);
        }
        catch (IOException e)
            {System.out.println("save test cubic Bezier error = " + e);}
    }

    private static double solve_quartic(double lead, double qua, double qub, double quc, double qud, double sgn0, double sgn1)
    {
        double sol, R, D, E;
        double delx0, delx1;

        qua /= lead;
        qub /= lead;
        quc /= lead;
        qud /= lead;
        if (IS_DEBUG)
            System.out.println("quartic      a,b,c,d = " + (float)qua + ", " + (float)qub + ", " + (float)quc + ", " + (float)qud + ", " + (float)sgn0 + ", " + (float)sgn1);
        sol = solve_cubic(-qub, qua*quc - 4*qud, -qua*qua*qud + 4*qub*qud - quc*quc);
        System.out.println("slope y = " + sol + ", " + (2*sol*sol - 4*qub*sol/3 - 8*qud/3));
        System.out.println("quartic = " + sol + ", " + ((sol - qub)*(sol + qub)*(sol + qub) - 4*quc*quc));
//        System.out.println("test y = " + sol + ", " + (sol*(12*qud + qub*qub) - qub*(4*qud + 3*qub*qub) - 9*quc*quc));
//        System.out.println("cubic sol = " + sol + ", " + (sol*sol*sol + qub*sol*sol - qub*qub*sol - qub*qub*qub - 4*quc*quc));
//        System.out.print("coalesced quad : ");
//        solve_cubic(qub, -qub*qub, -qub*qub*qub - 4*quc*quc); // fix fix dummy call temporary code
//        System.out.print("derivative of quad : ");
//        solve_cubic(3*qua/4, 2*qub/4, quc/4);                 // fix fix dummy call temporary code
//        double qa = 3;
//        double qb = -2*qub;
//        double qc = qua*quc - 4*qud;
//        System.out.println("derivative of resolvent cubic : " + (-qb - Math.sqrt(qb*qb - 4*qa*qc))/2/qa); // fix fix dummy call temporary code
        R = Math.sqrt(qua*qua/4 - qub + sol);
        D = Math.sqrt(3*qua*qua/4 - R*R - 2*qub + (4*qua*qub - 8*quc - qua*qua*qua)/4/R);
        E = Math.sqrt(3*qua*qua/4 - R*R - 2*qub - (4*qua*qub - 8*quc - qua*qua*qua)/4/R);
        if (IS_DEBUG)
        {
//          System.out.println("\ncubic sol = " + sol + ", " + R + ", " + D + ", " + E);
//          System.out.println(t1/2/Math.PI + ", " + t2/2/Math.PI);
            System.out.print("root 1 = "); check_quartic(-qua/4 + R/2 + D/2);
            System.out.print("root 2 = "); check_quartic(-qua/4 + R/2 - D/2);
            System.out.print("root 3 = "); check_quartic(-qua/4 - R/2 + E/2);
            System.out.print("root 4 = "); check_quartic(-qua/4 - R/2 - E/2);
        }
        if (!Double.isNaN(D) && (sgn0 == Math.signum(-qua/4 + R/2 + D/2)))
        {
            if (IS_DEBUG)
                System.out.println("using root 1 = " + sgn0 + ", " + (-qua/4 + R/2 + D/2));
            return (-qua/4 + R/2 + D/2);
        }
        else if (!Double.isNaN(E) && (sgn0 == Math.signum(-qua/4 - R/2 - E/2)))
        {
            if (IS_DEBUG)
                System.out.println("using root 4 = " + sgn0 + ", " + (-qua/4 - R/2 - E/2));
            return (-qua/4 - R/2 - E/2);
        }
        else if (!Double.isNaN(E) && (sgn0 == Math.signum(-qua/4 - R/2 + E/2)))
        {
            if (IS_DEBUG)
                System.out.println("using root 3 = " + sgn0 + ", " + (-qua/4 - R/2 + E/2));
            return (-qua/4 - R/2 + E/2);
        }
        else if (!Double.isNaN(D) && (sgn0 == Math.signum(-qua/4 + R/2 - D/2)))
        {
            if (IS_DEBUG)
                System.out.println("using root 2 = " + sgn0 + ", " + (-qua/4 + R/2 - D/2));
            return (-qua/4 + R/2 - D/2);
        }

        delx0 = -qua/4 + R/2 + D/2;                     // alternate check of root 1
        delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delx0*delx0/2)/2/mrot0;
        if (!Double.isNaN(D) && (sgn1 == Math.signum(delx1)))
        {
            if (IS_DEBUG)
                System.out.println("alternate root 1 = " + sgn1 + ", " + delx0);
            return delx0;
        }
        delx0 = -qua/4 - R/2 - E/2;                     // alternate check of root 4
        delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delx0*delx0/2)/2/mrot0;
        if (!Double.isNaN(E) && (sgn1 == Math.signum(delx1)))
        {
            if (IS_DEBUG)
                System.out.println("alternate root 4 = " + sgn1 + ", " + delx0);
            return delx0;
        }
        delx0 = -qua/4 - R/2 + E/2;                     // alternate check of root 3
        delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delx0*delx0/2)/2/mrot0;
        if (!Double.isNaN(E) && (sgn1 == Math.signum(delx1)))
        {
            if (IS_DEBUG)
                System.out.println("alternate root 3 = " + sgn1 + ", " + delx0);
            return delx0;
        }
        delx0 = -qua/4 + R/2 - D/2;                     // alternate check of root 2
        delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delx0*delx0/2)/2/mrot0;
        if (!Double.isNaN(D) && (sgn1 == Math.signum(delx1)))
        {
            if (IS_DEBUG)
                System.out.println("alternate root 2 = " + sgn1 + ", " + delx0);
            return delx0;
        }
        System.out.println("general quartic : Bad solution = " + R + ", " + D + ", " + E);
        return Double.NaN;
    }

    private static void check_quartic(double delx0)
    {
        // calculate other quartic root (delx1), and cross-check
        double delx1 = -(yrot3 - yrot0 - mrot0*(xrot3 - xrot0) - 3*Cu[0]*delx0*delx0/2)/2/mrot0;
        System.out.print((float)delx0 + ", " + (float)delx1);
        System.out.println(", " + (float)(-yrot3 + yrot0 - mrot0*(xrot3 - xrot0) + delx0*2*mrot0 - 3*Cu[1]*delx1*delx1/2));
    }

    private static double solve_cubic(double p, double q, double r)
    {
        // see Math CRC book, page 392

        double cua = (3*q - p*p)/3;
        double cub = (2*p*p*p - 9*p*q + 27*r)/27;
        double cud = cub*cub/4 + cua*cua*cua/27;

//        System.out.println("\ncubic p,q,r = " + p + ", " + q + ", " + r);
//        System.out.println("\ncubic a,b,d = " + cua + ", " + cub + ", " + cud);
        if (cud < 0)
        {
            double phi = Math.acos(-cub/2/Math.sqrt(-cua*cua*cua/27));
            System.out.println("3 cubic d < 0 : " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3 + 2*Math.PI/3) - p/3) + ", " + (2*Math.sqrt(-cua/3)*Math.cos(phi/3 + 4*Math.PI/3) - p/3));
//            return 2*Math.sqrt(-cua/3)*Math.cos(phi/3 + 2*Math.PI/3) - p/3;       // original code - KEEP
            return 2*Math.sqrt(-cua/3)*Math.cos(phi/3) - p/3;           // fix fix temporary code
        }
        else
        {
            System.out.println("1 cubic d > 0 : " + (Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3));
            return Math.cbrt(-cub/2 + Math.sqrt(cud)) + Math.cbrt(-cub/2 - Math.sqrt(cud)) - p/3;
        }
    }

    private static double getrotX(double m_x, double m_y, double m_theta)
    {
        return m_x*Math.cos(m_theta) + m_y*Math.sin(m_theta);
    }

    private static double getrotY(double m_x, double m_y, double m_theta)
    {
        return -m_x*Math.sin(m_theta) + m_y*Math.cos(m_theta);
    }
}
