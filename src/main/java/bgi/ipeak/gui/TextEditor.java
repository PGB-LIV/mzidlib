package bgi.ipeak.gui;

import java.awt.Component;
import java.util.EventObject;

import javax.swing.AbstractCellEditor;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.TableCellEditor;

public class TextEditor extends AbstractCellEditor implements TableCellEditor {

    private static final long serialVersionUID = -1L;
    private JLabel label;

    public TextEditor() {
        label = new JLabel();
    }

    @Override
    public String getCellEditorValue() {
        return label.getText();
    }

    @Override
    public boolean isCellEditable(EventObject anEvent) {
        return false;
    }

    @Override
    public Component getTableCellEditorComponent(JTable table, Object value,
            boolean isSelected, int row, int column) {
        label.setText(value == null ? "" : (String) value);
        return label;
    }

}
