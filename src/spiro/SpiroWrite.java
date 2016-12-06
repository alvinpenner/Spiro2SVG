
package spiro;

import java.awt.geom.*;
import java.awt.Shape;
import java.io.*;
import javax.swing.JOptionPane;

public final class SpiroWrite
{
    private static final double TOL = 0.00001;
    private static double numrot;                   // number of rotations
    private static int segs;                        // number of segments per rotation
    private static boolean isspiro;

    protected static void write_svg_file(String fname, boolean m_isspiro)
    {
        isspiro = m_isspiro;
        String strPath;
        if (fname.isEmpty())
        {
            JOptionPane.showMessageDialog(null, "No output file name was specified.\nPlease try again." , " No output location ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (!fname.endsWith(".svg"))
        {
            JOptionPane.showMessageDialog(null, "Output file type must be (*.svg).\nPlease try again." , " Wrong file type ", JOptionPane.WARNING_MESSAGE);
            return;
        }
        final File file = new File(fname);
        if (file.exists())
            if (JOptionPane.showConfirmDialog(null, "The file '" + file.getAbsolutePath() + "' already exists.\nDo you wish to over-write it?" , " SVG file exists ", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE) != JOptionPane.YES_OPTION)
                return;

        try
        {
            FileWriter out = new FileWriter(file);
            out.write(main.strHdr);
            strPath = "  <g\n"
                    + "     inkscape:label='Spiro2SVG'\n"
                    + "     inkscape:groupmode='layer'>\n";
            out.write(strPath);
            if (main.CanvasColor != 0)                          // background color
            {
                strPath = "    <rect\n"
                        + "        x='0.0'\n"
                        + "        y='0.0'\n"
                        + "        width='" + main.PAGE_WIDTH*main.mm2px + "'\n"
                        + "        height='" + main.PAGE_HEIGHT*main.mm2px + "'\n"
                        + "        style='"
                        + getspiroColor("fill", main.CanvasColor)
                        + "' />\n";
                out.write(strPath);
            }
            for (int iTemp = 0; iTemp < main.rowData[0].length; iTemp++)
            {
                float stroke_width = isspiro ? Float.parseFloat(main.rowData[10][iTemp])                        // CurvePenWidth
                                             : Float.parseFloat(main.rowData[main.rowData.length - 4][iTemp]);  // linewidth
                if (main.rowData[main.rowData.length - 1][iTemp].equals("Points"))                              // Edit Drawing Style
                    stroke_width = 1;
                strPath = "    <path\n"
                        + "        d='";
                out.write(strPath);
                if (isspiro)                                            // .spiro file
                    convertspiroParms(iTemp, out);
                else if (main.rowData.length == 13)                     // SpiroJ file
                    convertSpiroJParms(iTemp, out);
                else if (main.rowData.length == 14)                     // Farris Wheel file
                    convertFarrisParms(iTemp, out);
                strPath = "'\n"
                        + "        style='"
                        + (isspiro ? getspiroColor("fill", Integer.parseInt(main.rowData[13][iTemp]))           // FillArgb
                                   : getSpiroJColor("fill", main.rowData[main.rowData.length - 2][iTemp]))      // fillcolor
                        + (isspiro ? getspiroColor("stroke", Integer.parseInt(main.rowData[12][iTemp]))         // Argb
                                   : getSpiroJColor("stroke", main.rowData[main.rowData.length - 3][iTemp]))    // linecolor
                        + getStrokeWidth(stroke_width)                                                          // CurvePenWidth
                        + getFillRule(isspiro ? Integer.parseInt(main.rowData[14][iTemp]) : 0)                  // FillMode
                        + "stroke-miterlimit:" + main.MITER_LIMIT + "' />\n";                                   // arbitrary number
                out.write(strPath);
            }
//            SpiroCalc.write_test_quadratic(out);                    // fix fix for testing only
//            main.write_test_cubic(out);                             // fix fix for testing only
            out.write("  </g>\n");
            out.write(main.strFtr);
            out.close();
        }
        catch (IOException e)
            {System.out.println("saveSVG error = " + e);}
    }

    private static String getspiroColor(String type, long clr)
    {
        if (clr == 0) return type + ":none;";
        if (clr < 0) clr += 4294967296L;                            // 32 bit color
        if (main.FIT_POINTS_ONLY)               // show only the t_values for one lobe
            clr = 4294901760L;                  // red
        return type + ":#" + String.format("%06x", clr % 16777216) + ";"
             + type + "-opacity:" + (clr/16777216)/255.0F + ";";
    }

    private static String getSpiroJColor(String type, String clr)
    {
        if (clr.isEmpty()) return type + ":none;";
        return type + ":" + clr + ";"
             + type + "-opacity:1.0;";
    }

    private static String getStrokeWidth(float w)
    {
        return "stroke-width:" + w + "px;";
    }

    private static String getFillRule(int w)
    {
        return w == 0 ? "fill-rule:evenodd;" : "fill-rule:nonzero;";
    }

    private static void convertspiroParms(int index, FileWriter out)
    {
        // convert from Class 3 parameters - Spirograph 1.0.2.1
        // to Class 4 parameters - standard http://turnbull.mcs.st-and.ac.uk/

        double a, b, c;                          // standard spirograph parameters
        double newN = 0;                                            // new a/b
        double slide = Double.parseDouble(main.rowData[4][index])/100; // RotorSlide
        numrot = Integer.parseInt(main.rowData[2][index]);          // NumRotations
        segs = Integer.parseInt(main.rowData[3][index]);            // AnglesPerCycle

        if (main.rowData[9][index].equals("true"))                  // Lock
            main.rowData[8][index] = main.rowData[1][index];        // force hypocycloid behavior
        if (Double.parseDouble(main.rowData[8][index]) < 0)         // epiTrochoid
        {
            a = Double.parseDouble(main.rowData[0][index])          // StatorRadius
              - 2*Double.parseDouble(main.rowData[1][index]);
            b = Double.parseDouble(main.rowData[1][index]);         // RotorRadius
            c = -Double.parseDouble(main.rowData[8][index]);        // PenDistance
        }
        else                                                        // hypoTrochoid
        {
            a = Double.parseDouble(main.rowData[0][index]);         // StatorRadius
            b = -Double.parseDouble(main.rowData[1][index]);        // RotorRadius
            c = Double.parseDouble(main.rowData[8][index]);         // PenDistance
        }

        if (b == 0)
            System.out.println("Bad data at index " + index + " : " + a + ", " + b + ", " + c);
        else                                                        // compensate for slide
            newN = slide > 0 ? a/b + (a/b + 1)*slide : a/b + (a/b + 1)*slide/(1 - slide);
        b = (a + b)/(newN + 1);                                     // after sliding
        a = newN*b;                                                 // after sliding
        if (-b > a)                                                 // fix phase of hypo
        {                                                           // this is a custom fit, to satisfy .spiro format
            a = -a;
            b = -b;
        }
        if ((c == 0) || (a == 0))                                   // force a circle with 4 Bezier nodes
        {                                                           // this is a general fit, independent of .spiro
            b = -a - b - c;
            a = -2*b;
            c = 0;
        }
        if (a < 0)                                                  // transform to alternate representation
        {                                                           // this is a general fit, independent of .spiro
            double c_over_b = c/b;
            numrot = round2int(numrot*(1 + a/b));
            c = round2int(a + b);
            b = round2int(c*c_over_b);
            a = round2int(-a*c_over_b);
        }
        getPath(out, getspiroShape(index, a, b, c).getPathIterator(null));
    }

    private static void convertSpiroJParms(int index, FileWriter out)
    {
        // convert from Class 5 parameters - SpiroJ parameters
        // to Class 4 parameters - standard http://turnbull.mcs.st-and.ac.uk/
        // - alternatively -
        // render as generic roulette shapes

        numrot = 0;                                                         // number of rotations
        segs = Integer.parseInt(main.rowData[8][index]);                    // generator_steps
        if ((Math.abs(Float.parseFloat(main.rowData[0][index])) == Math.abs(Float.parseFloat(main.rowData[1][index])))
        &&  (Math.abs(Float.parseFloat(main.rowData[2][index])) == Math.abs(Float.parseFloat(main.rowData[3][index])))
        &&  (Math.abs(Float.parseFloat(main.rowData[4][index])) == Math.abs(Float.parseFloat(main.rowData[5][index])))
        &&  (Math.abs(Float.parseFloat(main.rowData[6][index])) == Math.abs(Float.parseFloat(main.rowData[7][index]))))
        {
            double ratio = Double.parseDouble(main.rowData[7][index])          // Frequency_y2
                         / Double.parseDouble(main.rowData[3][index]) - 1;     // Frequency_y1
            double b = Double.parseDouble(main.rowData[0][index])/(1 + ratio); // Radius_x1
            double a = b*ratio;
            double c = Double.parseDouble(main.rowData[4][index]);             // Radius_x2
            for (int i = 1; i < 100; i++)
                if (Math.round(i*ratio) == round2int(i*ratio))
                {
                    numrot = i;
                    break;
                }
            if (numrot > 0)
                getPath(out, getspiroShape(index, a, b, c).getPathIterator(null));
            else
                System.out.println("SpiroJ error: could not determine number of rotations");
        }
        else
            getPath(out, getSpiroJShape(Double.parseDouble(main.rowData[0][index]),
                                        Double.parseDouble(main.rowData[1][index]),
                                        Double.parseDouble(main.rowData[2][index]),
                                        Double.parseDouble(main.rowData[3][index]),
                                        Double.parseDouble(main.rowData[4][index]),
                                        Double.parseDouble(main.rowData[5][index]),
                                        Double.parseDouble(main.rowData[6][index]),
                                        Double.parseDouble(main.rowData[7][index])).getPathIterator(null));
    }

    private static void convertFarrisParms(int index, FileWriter out)
    {
        // render as three generic Farris Wheels

        numrot = 0;                                                         // not used
        segs = Integer.parseInt(main.rowData[9][index]);                    // generator_steps
        getPath(out, getFarrisShape(Double.parseDouble(main.rowData[0][index]),
                                    Double.parseDouble(main.rowData[1][index]),
                                    Double.parseDouble(main.rowData[2][index])*Math.PI/180,
                                    Double.parseDouble(main.rowData[3][index]),
                                    Double.parseDouble(main.rowData[4][index]),
                                    Double.parseDouble(main.rowData[5][index])*Math.PI/180,
                                    Double.parseDouble(main.rowData[6][index]),
                                    Double.parseDouble(main.rowData[7][index]),
                                    Double.parseDouble(main.rowData[8][index])*Math.PI/180).getPathIterator(null));
    }

    private static double round2int(double fData)
    {
        // round only if close to an int
        if (Math.abs(fData - Math.round(fData)) < TOL)
            fData = Math.round(fData);
        return fData;
    }

    private static void getPath(FileWriter out, PathIterator pit)
    {
//      very NB - use 'out' here because it is much faster than string concatenation

        int type;
        float[] coord = new float[6];
        float xend = 0, yend = 0;

        while (!pit.isDone())
        {
            type = pit.currentSegment(coord);
            switch (type)
            {
            case PathIterator.SEG_MOVETO:
                xend = coord[0];
                yend = coord[1];
                try {out.write(" M " + coord[0] + "," + coord[1]);}
                catch (IOException e) {System.out.println("saveSVG error = " + e);}
                break;
            case PathIterator.SEG_LINETO:
                if (Math.abs(xend-coord[0]) > 0.001 || Math.abs(yend-coord[1]) > 0.001)
                    try {out.write(" L " + coord[0] + "," + coord[1]);}
                    catch (IOException e) {System.out.println("saveSVG error = " + e);}
                break;
            case PathIterator.SEG_QUADTO:
                xend = coord[2];
                yend = coord[3];
                try {out.write(" Q " + coord[0] + "," + coord[1] + " " + coord[2] + "," + coord[3]);}
                catch (IOException e) {System.out.println("saveSVG error = " + e);}
                break;
            case PathIterator.SEG_CUBICTO:
                xend = coord[4];
                yend = coord[5];
                try {out.write(" C " + coord[0] + "," + coord[1] + " " + coord[2] + "," + coord[3] + " " + coord[4] + "," + coord[5]);}
                catch (IOException e) {System.out.println("saveSVG error = " + e);}
                break;
            case PathIterator.SEG_CLOSE:
                try {out.write(" Z");}
                catch (IOException e) {System.out.println("saveSVG error = " + e);}
                break;
            default:
                break;
            }
            pit.next();
        }
    }

    private static Shape getspiroShape(int index, double a, double b, double c)
    {
        Path2D path = new Path2D.Float();
        Rectangle2D.Float rect;
        Line2D.Float line;
        AffineTransform trans = AffineTransform.getTranslateInstance(
                                main.PAGE_WIDTH*main.mm2px/2,
                                main.PAGE_WIDTH*main.mm2px/2);
        if (isspiro)
        {
            trans = AffineTransform.getTranslateInstance(
                    main.PAGE_WIDTH*main.mm2px/2 + Float.parseFloat(main.rowData[5][index]),    // OriginX
                    main.PAGE_WIDTH*main.mm2px/2 + Float.parseFloat(main.rowData[6][index]));   // OriginY
            Float angle = Float.parseFloat(main.rowData[7][index])/500;                         // InitialAngle
            trans.scale(Float.parseFloat(main.rowData[11][index])/100,                          // Zoom
                        Float.parseFloat(main.rowData[11][index])/100);
            if (a != 0 || b < 0)
                trans.rotate(b < 0 ? angle : angle*(1 + 2*b/a));            // hypo/epi toggle
            else
                System.out.println("Cannot apply rotation transform to a circle");
        }

        if (main.IS_DEBUG)
            System.out.println("std spiro parms = " + numrot + ", " + segs + ", " + a + ", " + b + ", " + c);
        if (numrot*segs == 0)
            return path;
        int i;
        if (main.rowData[main.rowData.length - 1][index].equals("Lines"))   // Edit Drawing Style
        {
            boolean closed = false;
            for (i = 1; i < numrot*segs; i++)
            {
                line = new Line2D.Float((float)((a + b)*Math.cos(2*Math.PI*(i - 1)/segs) + c*Math.cos(2*Math.PI*(i - 1)*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.sin(2*Math.PI*(i - 1)/segs) + c*Math.sin(2*Math.PI*(i - 1)*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.cos(2*Math.PI*i/segs) + c*Math.cos(2*Math.PI*i*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.sin(2*Math.PI*i/segs) + c*Math.sin(2*Math.PI*i*(a/b + 1)/segs)));
                path.append(line, true);
                if (((i + 1) % segs == 0)
                &&  ((Math.round((i + 1)*a/b/segs) == round2int((i + 1)*a/b/segs)) || (c == 0)))
                {
                    if (main.IS_DEBUG)
                        System.out.println("breaking out at i = " + i);
                    closed = true;
                    break;
                }
            }
            if (closed)
                path.closePath();
            else
            {
                i = Math.round((float)numrot*segs);
                line = new Line2D.Float((float)((a + b)*Math.cos(2*Math.PI*(i - 1)/segs) + c*Math.cos(2*Math.PI*(i - 1)*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.sin(2*Math.PI*(i - 1)/segs) + c*Math.sin(2*Math.PI*(i - 1)*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.cos(2*Math.PI*i/segs) + c*Math.cos(2*Math.PI*i*(a/b + 1)/segs)),
                                        (float)((a + b)*Math.sin(2*Math.PI*i/segs) + c*Math.sin(2*Math.PI*i*(a/b + 1)/segs)));
                path.append(line, true);
            }
        }
        else if (main.rowData[main.rowData.length - 1][index].equals("Points"))
            for (i = 0; i < numrot*segs + 1; i++)
            {
                rect = new Rectangle2D.Float(-0.5F + (float)((a + b)*Math.cos(2*Math.PI*i/segs) + c*Math.cos(2*Math.PI*i*(a/b + 1)/segs)),
                                             -0.5F + (float)((a + b)*Math.sin(2*Math.PI*i/segs) + c*Math.sin(2*Math.PI*i*(a/b + 1)/segs)),
                                              1,
                                              1);
                path.append(rect, false);
                //double numer = SpiroCalc.getdX(2*Math.PI*i/segs)*SpiroCalc.getd2Y(2*Math.PI*i/segs) - SpiroCalc.getdY(2*Math.PI*i/segs)*SpiroCalc.getd2X(2*Math.PI*i/segs); // fix fix generate curvature
                //double denom = SpiroCalc.getdX(2*Math.PI*i/segs)*SpiroCalc.getdX(2*Math.PI*i/segs) + SpiroCalc.getdY(2*Math.PI*i/segs)*SpiroCalc.getdY(2*Math.PI*i/segs);   // fix fix generate curvature
                //denom = -Math.pow(denom, 1.5);
                //System.out.println(i + "/" + segs + " = ," + numer + ", " + denom + ", " + numer/denom);
                if (((i + 1) % segs == 0)
                &&  ((Math.round((i + 1)*a/b/segs) == round2int((i + 1)*a/b/segs)) || (c == 0)))
                {
                    if (main.IS_DEBUG)
                        System.out.println("breaking out at i = " + i);
                    break;
                }
            }
        else if (main.rowData[main.rowData.length - 1][index].equals("Bezier"))
        {
            double[] t_values = new double[12];
            double t_new = 0, t_old = 0;
            int N = SpiroCalc.get_t_values(t_values, a, b, c);

            if (main.IS_DEBUG)
            {
                System.out.println("Converting to Beziers : " + N);
                if (N > 0)
                    for (i = 0; i < N; i++)
                        System.out.print(" " + t_values[i]);
                System.out.println("\n\nt sequence :");
            }
            if (main.FIT_POINTS_ONLY)               // show only the t_values for one lobe
                for (i = 0; i < N; i++)
                {
                    rect = new Rectangle2D.Float(-1.5F + (float)((a + b)*Math.cos(t_values[i]) + c*Math.cos(t_values[i]*(a/b + 1))),
                                                 -1.5F + (float)((a + b)*Math.sin(t_values[i]) + c*Math.sin(t_values[i]*(a/b + 1))),
                                                  3,
                                                  3);
                    path.append(rect, false);
                }
            else
            {
                for (i = 0; i < round2int(N*numrot*a/Math.abs(b) + 1); i++)
                {
                    t_new = t_values[i % N] + 2*(int)(i/N)*Math.PI*Math.abs(b)/a;
                    if (t_new > 2*Math.PI*numrot)
                        t_new = 2*Math.PI*numrot;
                    if (t_new > t_old)
                    {
                        if (main.IS_DEBUG)
                            System.out.println(((int)(i/N)) + ", " + (i % N) + " :\t " + t_old + "\t " + t_new);
                        path.append(SpiroCalc.getBezier(t_old, t_new), true);
                        t_old = t_new;
                    }
                }
                if (Math.round(numrot*a/b) == round2int(numrot*a/b))
                    path.closePath();
            }
        }
        return trans.createTransformedShape(path);
    }

    private static Shape getSpiroJShape(double rx1, double ry1, double wx1, double wy1, double rx2, double ry2, double wx2, double wy2)
    {
        Path2D path = new Path2D.Float();
        Rectangle2D.Float rect;
        Line2D.Float line;
        AffineTransform trans = AffineTransform.getTranslateInstance(
                                main.PAGE_WIDTH*main.mm2px/2,
                                main.PAGE_WIDTH*main.mm2px/2);

        if (main.IS_DEBUG)
            System.out.println("SpiroJ parms = " + segs + ", " + rx1 + ", " + ry1 + ", " + wx1 + ", " + wy1 + ", " + rx2 + ", " + ry2 + ", " + wx2 + ", " + wy2);
        if (segs == 0)
            return path;

        int i;
        if (main.rowData[main.rowData.length - 1][0].equals("Lines"))
        {
            for (i = 1; i < segs + 1; i++)
                if ((segs*(2*i*(int)wy1/segs) == 2*i*wy1)
                &&  (segs*(2*i*(int)wy2/segs) == 2*i*wy2)
                &&  (segs*(i*(int)wx1/segs) == i*wx1)
                &&  (segs*(i*(int)wx2/segs) == i*wx2))
                    path.closePath();
                else
                {
                    line = new Line2D.Float((float)(rx1*Math.cos(2*Math.PI*(i - 1)*wx1/segs) + rx2*Math.cos(2*Math.PI*(i - 1)*wx2/segs)),
                                            (float)(ry1*Math.sin(2*Math.PI*(i - 1)*wy1/segs) + ry2*Math.sin(2*Math.PI*(i - 1)*wy2/segs)),
                                            (float)(rx1*Math.cos(2*Math.PI*i*wx1/segs) + rx2*Math.cos(2*Math.PI*i*wx2/segs)),
                                            (float)(ry1*Math.sin(2*Math.PI*i*wy1/segs) + ry2*Math.sin(2*Math.PI*i*wy2/segs)));
                    path.append(line, true);
                }
        }
        else if (main.rowData[main.rowData.length - 1][0].equals("Points"))
        {
            //double[] t_values = new double[100];            // fix fix just for testing purposes only
            //SpiroJCalc.get_t_values(t_values, rx1, ry1, wx1, wy1, rx2, ry2, wx2, wy2);  // fix fix testing only
            for (i = 0; i < segs; i++)
            {
                rect = new Rectangle2D.Float(-0.5F + (float)(rx1*Math.cos(2*Math.PI*i*wx1/segs) + rx2*Math.cos(2*Math.PI*i*wx2/segs)),
                                             -0.5F + (float)(ry1*Math.sin(2*Math.PI*i*wy1/segs) + ry2*Math.sin(2*Math.PI*i*wy2/segs)),
                                              1,
                                              1);
                path.append(rect, false);
                //double ptx = SpiroJCalc.getX(2*Math.PI*i/segs);     // fix fix testing only
                //double pty = SpiroJCalc.getY(2*Math.PI*i/segs);
                //double numer = SpiroJCalc.getdX(2*Math.PI*i/segs)*SpiroJCalc.getd2Y(2*Math.PI*i/segs) - SpiroJCalc.getdY(2*Math.PI*i/segs)*SpiroJCalc.getd2X(2*Math.PI*i/segs); // fix fix generate curvature
                //double denom = SpiroJCalc.getdX(2*Math.PI*i/segs)*SpiroJCalc.getdX(2*Math.PI*i/segs) + SpiroJCalc.getdY(2*Math.PI*i/segs)*SpiroJCalc.getdY(2*Math.PI*i/segs);   // fix fix generate curvature
                //denom = -Math.pow(denom, 1.5);
                //System.out.println(i + "/" + segs + " = ," + ptx + ", " + pty + ", " + numer + ", " + denom + ", " + numer/denom);
            }
        }
        else if (main.rowData[main.rowData.length - 1][0].equals("Bezier"))
        {
            double[] t_values = new double[500];
            int N = SpiroJCalc.get_t_values(t_values, rx1, ry1, wx1, wy1, rx2, ry2, wx2, wy2);
            t_values[N] = t_values[0] + 2*Math.PI;  // terminate the sequence, avoid rollover

            if (main.IS_DEBUG)
            {
                System.out.println("Converting to Beziers : " + N);
                if (N > 0)
                    for (i = 0; i <= N; i++)
                        System.out.print(" " + t_values[i]*180/Math.PI);
                System.out.println();
            }
//            for (i = 0; i < N; i++)
//            {
//                rect = new Rectangle2D.Float(-0.5F + rx1*(float)Math.cos(wx1*t_values[i]) + rx2*(float)Math.cos(wx2*t_values[i]),
//                                             -0.5F + ry1*(float)Math.sin(wy1*t_values[i]) + ry2*(float)Math.sin(wy2*t_values[i]),
//                                              1,
//                                              1);
//                path.append(rect, false);
//            }
            for (i = 0; i < N; i++)
                path.append(SpiroJCalc.getBezier(t_values[i], t_values[i + 1]), true);
            path.closePath();
        }
        return trans.createTransformedShape(path);
    }

    private static Shape getFarrisShape(double r1, double w1, double phi1, double r2, double w2, double phi2, double r3, double w3, double phi3)
    {
        Path2D path = new Path2D.Float();
        Rectangle2D.Float rect;
        Line2D.Float line;
        AffineTransform trans = AffineTransform.getTranslateInstance(
                                main.PAGE_WIDTH*main.mm2px/2,
                                main.PAGE_WIDTH*main.mm2px/2);

        if (main.IS_DEBUG)
            System.out.println("Farris parms = " + segs + ", " + r1 + ", " + w1 + ", " + phi1 + ", " + r2 + ", " + w2 + ", " + phi2 + ", " + r3 + ", " + w3 + ", " + phi3);
        if (segs == 0)
            return path;

        int i;
        if (main.rowData[main.rowData.length - 1][0].equals("Lines"))
        {
            for (i = 1; i < segs; i++)
            {
                line = new Line2D.Float((float)(r1*Math.cos(2*Math.PI*(i - 1)*w1/segs + phi1) + r2*Math.cos(2*Math.PI*(i - 1)*w2/segs + phi2) + r3*Math.cos(2*Math.PI*(i - 1)*w3/segs + phi3)),
                                        (float)(r1*Math.sin(2*Math.PI*(i - 1)*w1/segs + phi1) + r2*Math.sin(2*Math.PI*(i - 1)*w2/segs + phi2) + r3*Math.sin(2*Math.PI*(i - 1)*w3/segs + phi3)),
                                        (float)(r1*Math.cos(2*Math.PI*i*w1/segs + phi1) + r2*Math.cos(2*Math.PI*i*w2/segs + phi2) + r3*Math.cos(2*Math.PI*i*w3/segs + phi3)),
                                        (float)(r1*Math.sin(2*Math.PI*i*w1/segs + phi1) + r2*Math.sin(2*Math.PI*i*w2/segs + phi2) + r3*Math.sin(2*Math.PI*i*w3/segs + phi3)));
                path.append(line, true);
            }
            path.closePath();
        }
        else if (main.rowData[main.rowData.length - 1][0].equals("Points"))
        {
            for (i = 0; i < segs; i++)
            {
                rect = new Rectangle2D.Float(-0.5F + (float)(r1*Math.cos(2*Math.PI*i*w1/segs + phi1) + r2*Math.cos(2*Math.PI*i*w2/segs + phi2) + r3*Math.cos(2*Math.PI*i*w3/segs + phi3)),
                                             -0.5F + (float)(r1*Math.sin(2*Math.PI*i*w1/segs + phi1) + r2*Math.sin(2*Math.PI*i*w2/segs + phi2) + r3*Math.sin(2*Math.PI*i*w3/segs + phi3)),
                                              1,
                                              1);
                path.append(rect, false);
            }
        }
        else if (main.rowData[main.rowData.length - 1][0].equals("Bezier"))
        {
            double[] t_values = new double[1000];
            int N = FarrisCalc.get_t_values(t_values, r1, w1, phi1, r2, w2, phi2, r3, w3, phi3);
            t_values[N] = t_values[0] + 2*Math.PI;  // terminate the sequence, avoid rollover

            if (main.IS_DEBUG)
            {
                System.out.println("Converting to Beziers : " + N);
                if (N > 0)
                    for (i = 0; i <= N; i++)
                        System.out.print(" " + (float) (t_values[i]*180/Math.PI));
                System.out.println();
            }
//            for (i = 0; i < N; i++)
//            {
//                rect = new Rectangle2D.Float(-0.5F + (float)(r1*Math.cos(w1*t_values[i] + phi1) + r2*Math.cos(w2*t_values[i] + phi2) + r3*Math.cos(w3*t_values[i] + phi3)),
//                                             -0.5F + (float)(r1*Math.sin(w1*t_values[i] + phi1) + r2*Math.sin(w2*t_values[i] + phi2) + r3*Math.sin(w3*t_values[i] + phi3)),
//                                              1,
//                                              1);
//                path.append(rect, false);
//            }
            for (i = 0; i < N; i++)
                path.append(FarrisCalc.getBezier(t_values[i], t_values[i + 1]), true);
            path.closePath();
        }
        return trans.createTransformedShape(path);
    }
}
