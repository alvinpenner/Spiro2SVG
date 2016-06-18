
package spiro;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;

public final class SpiroConfig
{
    private static final Font fntBold = new Font(Font.DIALOG, Font.BOLD, 12);
    private static final Font fntPlain = new Font(Font.DIALOG, Font.PLAIN, 12);
    private static boolean ok = false;

    public static boolean showDialog(JDialog parent, String fname)
    {
        final JDialog cfg = new JDialog(parent, " Spiro2SVG v" + main.VERSION_NO + " : " + fname, true);

        final JTable jt = new JTable(new AbstractTableModel() {
            public int getRowCount()    { return main.rowData.length;  }
            public int getColumnCount() { return main.rowData[0].length; }

            public Object getValueAt(int r, int c) {
                if (r == 12 || r == 13)                 // fields Argb and FillArgb
                    return Integer.toHexString(Integer.parseInt(main.rowData[r][c])).toUpperCase();
                else
                    return main.rowData[r][c];
            }

            @Override public Class getColumnClass(int c) {
                return String.class;
            }

            @Override public boolean isCellEditable(int r, int c) {
                return r == main.rowData.length - 1;        // only last row is editable
            }

            @Override public void setValueAt(Object value, int r, int c) {
                main.rowData[r][c] = value.toString();
            }
        });

        jt.setFont(fntPlain);
        jt.getTableHeader().setFont(fntBold);
        jt.getTableHeader().setReorderingAllowed(false);
        jt.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        jt.setDefaultEditor(String.class, new DefaultCellEditor(new JComboBox(new String[] {main.STYLE_POINTS, main.STYLE_LINES, main.STYLE_BEZIER})));

        DefaultTableCellRenderer highlightRenderer = new DefaultTableCellRenderer() {
	    @Override public void setValue(Object value) {
                setText(value.toString());
	        if (value.equals("Points") || value.equals("Lines") || value.equals("Bezier"))
                    setBackground(new Color(128, 255, 128));        // light green
                else
                    setBackground(new Color(232, 232, 232));        // light gray
	    }
        };
        for (int i = 0; i < main.rowData[0].length; i++)
        {
            jt.getTableHeader().getColumnModel().getColumn(i).setHeaderValue(i);
            jt.getColumn(i).setCellRenderer(highlightRenderer);
        }

        jt.addKeyListener(new KeyAdapter()
        {
            @Override public void keyPressed(KeyEvent e)
            {
                if (e.getKeyCode() == KeyEvent.VK_F1)
                    {HelpAbout.showDialog(null);}
                else if (e.getKeyCode() == KeyEvent.VK_ENTER)
                {
                    ok = true;
                    cfg.dispose();
                }
            }
        });

        final JTable headerColumn = new JTable(main.rowNames, new String[] {""});
        headerColumn.setFont(fntBold);
        headerColumn.setBackground(jt.getTableHeader().getBackground());
        headerColumn.setEnabled(false);
        JViewport jtview = new JViewport();
        jtview.setView(headerColumn);
        jtview.setPreferredSize(new Dimension(120, 300));

        final JScrollPane jsp = new JScrollPane(jt);
        jsp.setRowHeader(jtview);
        jsp.setCorner(ScrollPaneConstants.UPPER_LEFT_CORNER, headerColumn.getTableHeader());

        JButton btnHelp = new JButton("Help");
        btnHelp.setPreferredSize(new Dimension(75, 25));
        btnHelp.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent event)
            {
                HelpAbout.showDialog(null);
            }
        });
        JButton btnOK = new JButton("OK");
        btnOK.setPreferredSize(new Dimension(75, 25));
        btnOK.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent event)
            {
                ok = true;
                cfg.dispose();
            }
        });
        JButton btnCancel = new JButton("Cancel");
        btnCancel.setPreferredSize(new Dimension(75, 25));
        btnCancel.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent event)
            {
                cfg.dispose();
            }
        });
        JLabel lblLeft = new JLabel("");
        lblLeft.setPreferredSize(new Dimension(10, 25));
        JLabel lblRight = new JLabel("");
        lblRight.setPreferredSize(new Dimension(10, 25));
        final JPanel pnl = new JPanel();
        pnl.setBorder(BorderFactory.createEmptyBorder(6,0,6,0));
        pnl.add(btnHelp);
        pnl.add(lblLeft);
        pnl.add(btnOK);
        pnl.add(lblRight);
        pnl.add(btnCancel);
        pnl.setPreferredSize(new Dimension(120, 50));
        cfg.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        cfg.getContentPane().setLayout(new BorderLayout());
        cfg.getContentPane().add(jsp, BorderLayout.CENTER);
        cfg.getContentPane().add(pnl, BorderLayout.SOUTH);
        cfg.getRootPane().setDefaultButton(btnOK);
        if (main.rowData[0].length < 3)
            cfg.setMinimumSize(new Dimension(140 + 75*3, 368));
        else
            cfg.setMinimumSize(new Dimension(140 + 75*main.rowData[0].length, 368));
        cfg.setLocationByPlatform(true);
        cfg.setVisible(true);
        return ok;
    }
}
