
package spiro;

import java.awt.geom.*;
import java.awt.Shape;
import java.io.*;
import javax.swing.JOptionPane;

public final class SpiroWrite
{
    private static float a, b, c;                   // standard spirograph parameters
    private static float numrot;                    // number of rotations
    private static int segs;                        // number of segments per rotation

    public static void write_svg_file(String fname)
    {
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
                        + getColor("fill"  , main.CanvasColor)
                        + "' />\n";
                out.write(strPath);
            }
            for (int iTemp = 0; iTemp < main.rowData[0].length; iTemp++)
            {
                float stroke_width = Float.parseFloat(main.rowData[10][iTemp]);         // CurvePenWidth
                if (main.rowData[15][iTemp].equals("Points"))                           // Edit Drawing Style
                    stroke_width = 1;
                strPath = "    <path\n"
                        + "        d='";
                out.write(strPath);
                getStandardParms(iTemp);
                getPath(out, iTemp);
                strPath = "'\n"
                        + "        style='"
                        + getColor("fill"  , Integer.parseInt(main.rowData[13][iTemp])) // FillArgb
                        + getColor("stroke", Integer.parseInt(main.rowData[12][iTemp])) // Argb
                        + getStrokeWidth(stroke_width)                                  // CurvePenWidth
                        + getFillRule(Integer.parseInt(main.rowData[14][iTemp]))        // FillMode
                        + "stroke-miterlimit:" + main.MITER_LIMIT + "' />\n";           // arbitrary number
                out.write(strPath);
            }
            // SpiroCalc.write_test_quadratic(out, a, b, c);        // fix fix for testing only
            // SpiroCalc.write_test_cubic(out, a, b, c);            // fix fix for testing only
            out.write("  </g>\n");
            out.write(main.strFtr);
            out.close();
//            JOptionPane.showMessageDialog(null, "Saved as '" + file.getPath() + "' ", " Export SVG", JOptionPane.INFORMATION_MESSAGE);
        }
        catch (IOException e)
            {System.out.println("saveSVG error = " + e);}
    }

    private static String getColor(String type, long clr)
    {
        if (clr == 0) return type + ":none;";
        if (clr < 0) clr += 4294967296L;                            // 32 bit color
        return type + ":#" + String.format("%06x", clr % 16777216) + ";"
             + type + "-opacity:" + (clr/16777216)/255.0F + ";";
    }

    private static String getStrokeWidth(float w)
    {
        return "stroke-width:" + w + "px;";
    }

    private static String getFillRule(int w)
    {
        return w == 0 ? "fill-rule:evenodd;" : "fill-rule:nonzero;";
    }

    private static void getStandardParms(int index)
    {
        // convert from Class 3 parameters - Spirograph 1.0.2.1
        // to Class 4 parameters - standard http://turnbull.mcs.st-and.ac.uk/

        float newN = 0;                                             // new a/b
        float slide = Float.parseFloat(main.rowData[4][index])/100; // RotorSlide
        numrot = Integer.parseInt(main.rowData[2][index]);          // NumRotations
        segs = Integer.parseInt(main.rowData[3][index]);            // AnglesPerCycle

        if (main.rowData[9][index].equals("true"))                  // Lock
            main.rowData[8][index] = main.rowData[1][index];        // force hypocycloid behavior
        if (Float.parseFloat(main.rowData[8][index]) < 0)           // epiTrochoid
        {
            a = Float.parseFloat(main.rowData[0][index])            // StatorRadius
              - 2*Float.parseFloat(main.rowData[1][index]);
            b = Float.parseFloat(main.rowData[1][index]);           // RotorRadius
            c = Float.parseFloat(main.rowData[8][index]);           // PenDistance
        }
        else                                                        // hypoTrochoid
        {
            a = Float.parseFloat(main.rowData[0][index]);           // StatorRadius
            b = -Float.parseFloat(main.rowData[1][index]);          // RotorRadius
            c = -Float.parseFloat(main.rowData[8][index]);          // PenDistance
        }

        if (b == 0)
            System.out.println("Bad data at index " + index + " : " + a + ", " + b + ", " + c);
        else                                                        // compensate for slide
            newN = slide > 0 ? a/b + (a/b + 1)*slide : a/b + (a/b + 1)*slide/(1 - slide);
        b = (a + b)/(newN + 1);                                     // after sliding
        a = newN*b;                                                 // after sliding
        if ((c == 0) || (a == 0))                                   // force a circle with 4 Bezier nodes
        {                                                           // this is a general fit, independent of .spiro
            b = -a - b + c;
            a = -2*b;
            c = 0;
        }
        if (-b > a)                                                 // fix phase of hypo
        {                                                           // this is a custom fit, to satisfy .spiro format
            System.out.println("changing sign");
            a = -a;
            b = -b;
        }
        if (a < 0)                                                  // transform to alternate representation
        {                                                           // this is a general fit, independent of .spiro
            float c_over_b = c/b;
            System.out.println("alternate representation numrot = " + a + ", " + b + ", " + c + ", " + numrot);
            numrot = round2int(numrot*(1 + a/b));
            c = round2int(-a - b);
            b = round2int(c*c_over_b);
            a = round2int(a*c_over_b);
            System.out.println("new       representation numrot = " + a + ", " + b + ", " + c + ", " + numrot);
        }
    }

    private static float round2int(float fData)
    {
        // round only if close to an int
        if (Math.abs(fData - Math.round(fData)) < 0.00001)
            fData = Math.round(fData);
        return fData;
    }

    private static void getPath(FileWriter out, int index)
    {
//      very NB - use 'out' here because it is much faster than string concatenation

        int type;
        float[] coord = new float[6];
        float xend = 0, yend = 0;
        PathIterator pit = getShape(index).getPathIterator(null);

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

    private static Shape getShape(int index)
    {
        Path2D path = new Path2D.Float();
        Rectangle2D.Float rect;
        Line2D.Float line;
        AffineTransform trans = AffineTransform.getTranslateInstance(
                                main.PAGE_WIDTH*main.mm2px/2 + Float.parseFloat(main.rowData[5][index]),    // OriginX
                                main.PAGE_WIDTH*main.mm2px/2 + Float.parseFloat(main.rowData[6][index]));   // OriginY
        Float angle = Float.parseFloat(main.rowData[7][index])/500;         // InitialAngle

        trans.scale(Float.parseFloat(main.rowData[11][index])/100,          // Zoom
                    Float.parseFloat(main.rowData[11][index])/100);
        if (a != 0 || b < 0)
            trans.rotate(b < 0 ? angle : angle*(1 + 2*b/a));                // hypo/epi toggle
        else
            System.out.println("Cannot apply rotation transform to a circle");

        System.out.println("std parms = " + numrot + ", " + segs + ", " + a + ", " + b + ", " + c + ", " + angle);
        if (numrot*segs == 0)
            return null;
        int i;
        if (main.rowData[15][index].equals("Lines"))                        // Edit Drawing Style
        {
            boolean closed = false;
            for (i = 1; i < numrot*segs; i++)
            {
                line = new Line2D.Float((a + b)*(float)Math.cos(2*Math.PI*(i - 1)/segs) - c*(float)Math.cos(2*Math.PI*(i - 1)*(a/b + 1)/segs),
                                        (a + b)*(float)Math.sin(2*Math.PI*(i - 1)/segs) - c*(float)Math.sin(2*Math.PI*(i - 1)*(a/b + 1)/segs),
                                        (a + b)*(float)Math.cos(2*Math.PI*i/segs) - c*(float)Math.cos(2*Math.PI*i*(a/b + 1)/segs),
                                        (a + b)*(float)Math.sin(2*Math.PI*i/segs) - c*(float)Math.sin(2*Math.PI*i*(a/b + 1)/segs));
                path.append(line, true);
                if (((i + 1) % segs == 0)
                &&  ((Math.round((i + 1)*a/b/segs) == round2int((i + 1)*a/b/segs)) || (c == 0)))
                {
                    System.out.println("breaking out at i = " + i);
                    closed = true;
                    break;
                }
            }
            if (closed)
                path.closePath();
            else
            {
                i = Math.round(numrot*segs);
                line = new Line2D.Float((a + b)*(float)Math.cos(2*Math.PI*(i - 1)/segs) - c*(float)Math.cos(2*Math.PI*(i - 1)*(a/b + 1)/segs),
                                        (a + b)*(float)Math.sin(2*Math.PI*(i - 1)/segs) - c*(float)Math.sin(2*Math.PI*(i - 1)*(a/b + 1)/segs),
                                        (a + b)*(float)Math.cos(2*Math.PI*i/segs) - c*(float)Math.cos(2*Math.PI*i*(a/b + 1)/segs),
                                        (a + b)*(float)Math.sin(2*Math.PI*i/segs) - c*(float)Math.sin(2*Math.PI*i*(a/b + 1)/segs));
                path.append(line, true);
            }
        }
        else if (main.rowData[15][index].equals("Points"))
            for (i = 0; i < numrot*segs + 1; i++)
            {
                rect = new Rectangle2D.Float(-0.5F + (a + b)*(float)Math.cos(2*Math.PI*i/segs) - c*(float)Math.cos(2*Math.PI*i*(a/b + 1)/segs),
                                             -0.5F + (a + b)*(float)Math.sin(2*Math.PI*i/segs) - c*(float)Math.sin(2*Math.PI*i*(a/b + 1)/segs),
                                              1,
                                              1);
                path.append(rect, false);
                if (((i + 1) % segs == 0)
                &&  ((Math.round((i + 1)*a/b/segs) == round2int((i + 1)*a/b/segs)) || (c == 0)))
                {
                    System.out.println("breaking out at i = " + i);
                    break;
                }
            }
        else if (main.rowData[15][index].equals("Bezier"))
        {
            double[] t_values = new double[12];
            double t_new = 0, t_old = 0;
            int N = SpiroCalc.get_t_values(t_values, a, b, c);

            System.out.println("Converting to Beziers : " + N);
            if (N > 0)
                for (i = 0; i < N; i++)
                    System.out.print(" " + t_values[i]);
            System.out.println("\n\nt sequence :");
            for (i = 0; i < round2int(N*numrot*a/Math.abs(b) + 1); i++)
            {
                t_new = t_values[i % N] + 2*(int)(i/N)*Math.PI*Math.abs(b)/a;
                if (t_new > 2*Math.PI*numrot)
                    t_new = 2*Math.PI*numrot;
                if (t_new > t_old)
                {
                    System.out.println(((int)(i/N)) + ", " + (i % N) + " :\t " + t_old + "\t " + t_new);
                    path.append(SpiroCalc.getBezier(a, b, c, t_old, t_new), true);
                    t_old = t_new;
                }
            }
            if (Math.round(numrot*a/b) == round2int(numrot*a/b))
                path.closePath();
        }
        return trans.createTransformedShape(path);
    }
}
