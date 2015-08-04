package bgi.ipeak.gui;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.UIManager;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputMethodListener;
import java.beans.EventHandler;
import java.io.File;
import java.util.Date;
import javax.swing.JTextField;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JTextArea;
import javax.swing.border.TitledBorder;
import javax.swing.ScrollPaneConstants;

public class IPeakGui extends JFrame{
	@SuppressWarnings("unused")
	private JButton msgfButton, tandemButton, omssaButton, fastaButton, outputButton, modButton, userModButton,
					 checkAdvancedOptions;
	private static JButton stopButton, startButton;
	private JTextField msgfText, tandemText, omssaText, fastaText, decoyRegexText, modfileText, userModText,
						 FileRegexText, outputText, maxFdrText;
	private JCheckBox debug;
	private JFileChooser fc;
	private JPanel inputLogOptions, advancedOptions;
	private static final long serialVersionUID = -1L;
	public static final String lookAndFeel =UIManager.getSystemLookAndFeelClassName();
	private JComboBox<String> accSplitRegexText;
	private static JTextArea log;
	private RunIPeak runIPeak;
	private String delim="/ /";
	private String fdr_level ="PeptideLevel";
	private boolean bc;
//	private JLabel PercolatorPath;
//	private JTextField percolatorText;
//	private JButton btn_percolatorOpen;
	private JLabel lblFdrLevel;
	private JComboBox<String> FDRlevelComboBox;
	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                new IPeakGui();
            }
        });
	}

	/**
	 * Create the application.
	 */
	public IPeakGui() {
		setDefaultLookAndFeelDecorated(true);
		setSize(680, 710);
		//setPreferredSize(new Dimension(600,600));
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		try {
			UIManager.setLookAndFeel(lookAndFeel);
		} catch (Exception e) {
			e.printStackTrace();
		} 
		setTitle("IPeak");
		fc=new JFileChooser();
		fc.setMultiSelectionEnabled(false);
		initialize();
		pack();
		setLocationRelativeTo(null);
		setVisible(true);
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[] {600};
                //set each region height
		gridBagLayout.rowHeights = new int[] {300, 70, 70};
		gridBagLayout.columnWeights = new double[]{1.0};
		gridBagLayout.rowWeights = new double[]{0.0, 0.0, 0.0};
		getContentPane().setLayout(gridBagLayout);
		
		inputLogOptions=new JPanel();
		inputLogOptions.setBorder(new TitledBorder(null, "Input", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagLayout gbl_inputlog=new GridBagLayout();
		gbl_inputlog.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		gbl_inputlog.columnWidths = new int[] {500, 80};
		gbl_inputlog.columnWeights=new double[]{1.0, 0.0};
		inputLogOptions.setLayout(gbl_inputlog);
		GridBagConstraints gbc_inputLogOptions = new GridBagConstraints();
		gbc_inputLogOptions.weighty = 1.0;
		gbc_inputLogOptions.fill = GridBagConstraints.HORIZONTAL;
		gbc_inputLogOptions.gridx = 0;
		gbc_inputLogOptions.gridy = 0;
		getContentPane().add(inputLogOptions, gbc_inputLogOptions);
		
		
		JLabel msgf=new JLabel("MSGF+ (search result,file format is .mzid)");
		GridBagConstraints gbc_msgf=new GridBagConstraints();
		gbc_msgf.anchor=GridBagConstraints.WEST;
		gbc_msgf.fill=GridBagConstraints.NONE;
		gbc_msgf.gridx=0;
		gbc_msgf.gridy=0;
		gbc_msgf.insets=new Insets(0, 5, 5, 5);
		inputLogOptions.add(msgf, gbc_msgf);
		
		msgfText=new JTextField();
		GridBagConstraints gbc_msgfText=new GridBagConstraints();
		gbc_msgfText.anchor=GridBagConstraints.CENTER;
		gbc_msgfText.fill=GridBagConstraints.HORIZONTAL;
		gbc_msgfText.gridx=0;
		gbc_msgfText.gridy=1;
		gbc_msgfText.insets=new Insets(0, 10, 5, 5);
		EventHandler.create(InputMethodListener.class, msgfText, "text", "source.text");
		inputLogOptions.add(msgfText, gbc_msgfText);
		
		msgfButton=new JButton("Open");
		GridBagConstraints gbc_msgfbtn=new GridBagConstraints();
		gbc_msgfbtn.anchor=GridBagConstraints.WEST;
		gbc_msgfbtn.fill=GridBagConstraints.BOTH;
		gbc_msgfbtn.gridx=1;
		gbc_msgfbtn.gridy=1;
		gbc_msgfbtn.insets=new Insets(0, 0, 5, 0);
		msgfButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
				fc.setFileFilter(new ExtFilter("mzid"));
				String f=chooseFile();
				if(f!=null)
					msgfText.setText(f);
			}
		});
		inputLogOptions.add(msgfButton, gbc_msgfbtn);
		
		JLabel tandem=new JLabel("X!Tandem (search result,file format is *.xml)");
		GridBagConstraints gbc_tandem=new GridBagConstraints();
		gbc_tandem.anchor=GridBagConstraints.WEST;
		gbc_tandem.fill=GridBagConstraints.NONE;
		gbc_tandem.gridx=0;
		gbc_tandem.gridy=2;
		gbc_tandem.insets=new Insets(0,5,5,5);
		inputLogOptions.add(tandem, gbc_tandem);
		
		tandemText=new JTextField();
		GridBagConstraints gbc_tandemText=new GridBagConstraints();
		gbc_tandemText.anchor=GridBagConstraints.CENTER;
		gbc_tandemText.fill=GridBagConstraints.HORIZONTAL;
		gbc_tandemText.gridx=0;
		gbc_tandemText.gridy=3;
		gbc_tandemText.insets=new Insets(0,10,5,5);
		EventHandler.create(InputMethodListener.class, tandemText, "text", "source.text");
		inputLogOptions.add(tandemText, gbc_tandemText);
		
		tandemButton=new JButton("Open");
		GridBagConstraints gbc_tandembtn=new GridBagConstraints();
		gbc_tandembtn.insets = new Insets(0, 0, 5, 0);
		gbc_tandembtn.anchor=GridBagConstraints.WEST;
		gbc_tandembtn.fill=GridBagConstraints.BOTH;
		gbc_tandembtn.gridx=1;
		gbc_tandembtn.gridy=3;
		gbc_msgfbtn.insets=new Insets(0, 0, 5, 5);
		tandemButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
				fc.setFileFilter(new ExtFilter("xml/pepxml"));
				String f=chooseFile();
				if(f!=null)
					tandemText.setText(f);
			}
		});
		inputLogOptions.add(tandemButton, gbc_tandembtn);
		
		JLabel omssa=new JLabel("OMSSA (search result,file format is *.omx)");
		GridBagConstraints gbc_omssa=new GridBagConstraints();
		gbc_omssa.anchor=GridBagConstraints.WEST;
		gbc_omssa.fill=GridBagConstraints.NONE;
		gbc_omssa.gridx=0;
		gbc_omssa.gridy=4;
		gbc_omssa.insets=new Insets(0,5,5,5);
		inputLogOptions.add(omssa, gbc_omssa);
		
		omssaText =new JTextField();
		GridBagConstraints gbc_omssaText=new GridBagConstraints();
		gbc_omssaText.anchor=GridBagConstraints.CENTER;
		gbc_omssaText.fill=GridBagConstraints.HORIZONTAL;
		gbc_omssaText.gridx=0;
		gbc_omssaText.gridy=5;
		gbc_omssaText.insets=new Insets(0,10,5,5);
		EventHandler.create(InputMethodListener.class, omssaText, "text", "source.text");
		inputLogOptions.add(omssaText, gbc_omssaText);
		
		omssaButton=new JButton("Open");
		GridBagConstraints gbc_omssaButton=new GridBagConstraints();
		gbc_omssaButton.anchor=GridBagConstraints.WEST;
		gbc_omssaButton.fill=GridBagConstraints.BOTH;
		gbc_omssaButton.gridx=1;
		gbc_omssaButton.gridy=5;
		gbc_omssaButton.insets=new Insets(0, 0, 5, 0);
		omssaButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
				fc.setFileFilter(new ExtFilter("omx"));
				String f=chooseFile();
				if(f!=null)
					omssaText.setText(f);
			}
		});
		inputLogOptions.add(omssaButton, gbc_omssaButton);
		
		JLabel fasta=new JLabel("FASTA (database,file format *.fasta)");
		GridBagConstraints gbc_fasta=new GridBagConstraints();
		gbc_fasta.anchor=GridBagConstraints.WEST;
		gbc_fasta.fill=GridBagConstraints.NONE;
		gbc_fasta.gridx=0;
		gbc_fasta.gridy=6;
		gbc_fasta.insets=new Insets(0,5,5,5);
		inputLogOptions.add(fasta, gbc_fasta);
		
		fastaText=new JTextField();
		GridBagConstraints gbc_fastaText=new GridBagConstraints();
		gbc_fastaText.anchor=GridBagConstraints.CENTER;
		gbc_fastaText.fill=GridBagConstraints.BOTH;
		gbc_fastaText.gridx=0;
		gbc_fastaText.gridy=7;
		gbc_fastaText.insets=new Insets(0,10,5,5);
		EventHandler.create(InputMethodListener.class, fastaText, "text", "source.text");
		inputLogOptions.add(fastaText, gbc_fastaText);
		
		fastaButton=new JButton("Open");
		GridBagConstraints gbc_fastabtn=new GridBagConstraints();
		gbc_fastabtn.anchor=GridBagConstraints.WEST;
		gbc_fastabtn.fill=GridBagConstraints.BOTH;
		gbc_fastabtn.gridx=1;
		gbc_fastabtn.gridy=7;
		gbc_fastabtn.insets=new Insets(0, 0, 5, 0);
		fastaButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
				fc.setFileFilter(new ExtFilter("fasta"));
				String f=chooseFile();
				if(f!=null)
					fastaText.setText(f);
			}
		});
		inputLogOptions.add(fastaButton, gbc_fastabtn);
		
		JLabel output=new JLabel("Output Directory");
		GridBagConstraints gbc_output=new GridBagConstraints();
		gbc_output.anchor=GridBagConstraints.WEST;
		gbc_output.fill=GridBagConstraints.NONE;
		gbc_output.gridx=0;
		gbc_output.gridy=8;
		gbc_output.insets=new Insets(0,5,5,5);
		inputLogOptions.add(output, gbc_output);
		
		outputText=new JTextField();
		GridBagConstraints gbc_outputText=new GridBagConstraints();
		gbc_outputText.anchor=GridBagConstraints.CENTER;
		gbc_outputText.fill=GridBagConstraints.BOTH;
		gbc_outputText.gridx=0;
		gbc_outputText.gridy=9;
		gbc_outputText.insets=new Insets(0,10,5,5);
		EventHandler.create(InputMethodListener.class, outputText, "text", "source.text");
		inputLogOptions.add(outputText, gbc_outputText);
		
		outputButton=new JButton("Open");
		GridBagConstraints gbc_outputbtn=new GridBagConstraints();
		gbc_outputbtn.anchor=GridBagConstraints.WEST;
		gbc_outputbtn.fill=GridBagConstraints.BOTH;
		gbc_outputbtn.gridx=1;
		gbc_outputbtn.gridy=9;
		gbc_outputbtn.insets=new Insets(0, 0, 5, 0);
		outputButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				fc.resetChoosableFileFilters();
				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				fc.setFileFilter(null);
				String f=chooseFile();
				if(f!=null)
					outputText.setText(f);
			}
		});
		inputLogOptions.add(outputButton, gbc_outputbtn);
		
                
		JLabel fileName=new JLabel("Result File Name");
		GridBagConstraints gbc_filename=new GridBagConstraints();
		gbc_filename.anchor=GridBagConstraints.WEST;
		gbc_filename.fill=GridBagConstraints.NONE;
		gbc_filename.gridx=0;
		gbc_filename.gridy=10;
		gbc_filename.insets=new Insets(0,5,5,5);
		inputLogOptions.add(fileName,gbc_filename);
		
		FileRegexText = new JTextField("iPeak_Result_"+String.format("%tF", new Date()));
		GridBagConstraints gbc_filenameText=new GridBagConstraints();
		gbc_filenameText.anchor=GridBagConstraints.CENTER;
		gbc_filenameText.fill=GridBagConstraints.BOTH;
		gbc_filenameText.gridx=0;
		gbc_filenameText.gridy=11;
		gbc_filenameText.insets=new Insets(0, 5, 0, 5);
		inputLogOptions.add(FileRegexText,gbc_filenameText);
		EventHandler.create(InputMethodListener.class, FileRegexText, "text", "source.text");
		
		//collapsiblePane=new JXCollapsiblePane();
		//collapsiblePane.getContentPane().setLayout(new BorderLayout());
		advancedOptions=new JPanel();
		//advancedOptions=new JXCollapsiblePane();
		advancedOptions.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Advanced Options"));
		GridBagLayout gbl_advanced=new GridBagLayout();
		gbl_advanced.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0};
		gbl_advanced.columnWidths = new int[] {80, 180, 80, 80, 50};
		gbl_advanced.columnWeights = new double[]{0.0, 1.0, 0.0};
		advancedOptions.setLayout(gbl_advanced);
		//collapsiblePane.getContentPane().add("Center",advancedOptions);
		//collapsiblePane.setVisible(false);
		//getContentPane().add("Center",collapsiblePane);
		GridBagConstraints gbc_advancedOptions = new GridBagConstraints();
		gbc_advancedOptions.fill = GridBagConstraints.BOTH;
		gbc_advancedOptions.gridx = 0;
		gbc_advancedOptions.gridy = 1;
		getContentPane().add(advancedOptions, gbc_advancedOptions);
		
		
		lblFdrLevel = new JLabel("FDR Level");
		GridBagConstraints gbc_lblFdrLevel = new GridBagConstraints();
		gbc_lblFdrLevel.anchor = GridBagConstraints.EAST;
		gbc_lblFdrLevel.insets = new Insets(0, 0, 5, 5);
		gbc_lblFdrLevel.gridx = 0;
		gbc_lblFdrLevel.gridy = 0;
		advancedOptions.add(lblFdrLevel, gbc_lblFdrLevel);
		
		FDRlevelComboBox = new JComboBox<String>();
		FDRlevelComboBox.setEditable(true);
		FDRlevelComboBox.setModel(new DefaultComboBoxModel<String>(new String[] {"PeptideLevel", "PSMLevel"}));
		GridBagConstraints gbc_FDRlevelComboBox = new GridBagConstraints();
		gbc_FDRlevelComboBox.insets = new Insets(0, 0, 5, 5);
		gbc_FDRlevelComboBox.fill = GridBagConstraints.HORIZONTAL;
		gbc_FDRlevelComboBox.gridx = 1;
		gbc_FDRlevelComboBox.gridy = 0;
		FDRlevelComboBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unchecked")
				JComboBox<String> cb=(JComboBox<String>)e.getSource();
				fdr_level = (String)cb.getSelectedItem();
			}
		});
		advancedOptions.add(FDRlevelComboBox, gbc_FDRlevelComboBox);
                
		JLabel decoyRegex = new JLabel("Decoy Regular Expression");
		GridBagConstraints gbc_decoyRegex = new GridBagConstraints();
		gbc_decoyRegex.anchor = GridBagConstraints.CENTER;
		gbc_decoyRegex.insets = new Insets(0, 5, 5, 5);
		gbc_decoyRegex.gridx = 0;
		gbc_decoyRegex.gridy = 1;
		advancedOptions.add(decoyRegex, gbc_decoyRegex);
		
		decoyRegexText = new JTextField();
		decoyRegexText.setText("###REV###");
		GridBagConstraints gbc_decoyRegexText = new GridBagConstraints();
		gbc_decoyRegexText.insets = new Insets(0, 0, 5, 5);
		gbc_decoyRegexText.fill = GridBagConstraints.HORIZONTAL;
		gbc_decoyRegexText.gridx = 1;
		gbc_decoyRegexText.gridy = 1;
		EventHandler.create(InputMethodListener.class, decoyRegexText, "text", "source.text");
		advancedOptions.add(decoyRegexText, gbc_decoyRegexText);
		decoyRegexText.setColumns(10);
		
		JLabel accDelim = new JLabel("Accession Delimiter");
		GridBagConstraints gbc_accDelim = new GridBagConstraints();
		gbc_accDelim.anchor = GridBagConstraints.CENTER;
		gbc_accDelim.insets = new Insets(0, 5, 5, 5);
		gbc_accDelim.gridx = 0;
		gbc_accDelim.gridy = 2;
		advancedOptions.add(accDelim, gbc_accDelim);
		
		accSplitRegexText = new JComboBox<String>();
		accSplitRegexText.setEditable(true);
		accSplitRegexText.setModel(new DefaultComboBoxModel<String>(new String[] {"/ /", "/|/"}));
		GridBagConstraints gbc_accSplitRegexText = new GridBagConstraints();
		gbc_accSplitRegexText.insets = new Insets(0, 0, 5, 5);
		gbc_accSplitRegexText.fill = GridBagConstraints.HORIZONTAL;
		gbc_accSplitRegexText.gridx = 1;
		gbc_accSplitRegexText.gridy = 2;
		accSplitRegexText.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				@SuppressWarnings("unchecked")
				JComboBox<String> cb=(JComboBox<String>)e.getSource();
				delim = (String)cb.getSelectedItem();
			}
		});
		advancedOptions.add(accSplitRegexText, gbc_accSplitRegexText);
		
		JLabel fdr=new JLabel("False Discovery Rate");
		GridBagConstraints gbc_fdr=new GridBagConstraints();
		gbc_fdr.anchor=GridBagConstraints.CENTER;
		gbc_fdr.fill=GridBagConstraints.NONE;
		gbc_fdr.gridx=0;
		gbc_fdr.gridy=3;
		gbc_fdr.insets=new Insets(0, 5, 0, 5);
		advancedOptions.add(fdr, gbc_fdr);
		
		maxFdrText=new JTextField("0.01");
		GridBagConstraints gbc_fdrText=new GridBagConstraints();
		gbc_fdrText.anchor=GridBagConstraints.CENTER;
		gbc_fdrText.fill=GridBagConstraints.HORIZONTAL;
		gbc_fdrText.gridx=1;
		gbc_fdrText.gridy=3;
		gbc_fdrText.insets=new Insets(0, 0, 0, 5);
		EventHandler.create(InputMethodListener.class, maxFdrText, "text", "source.text");
		advancedOptions.add(maxFdrText, gbc_fdrText);
		
		JPanel logButton=new JPanel();
		GridBagLayout gbl_logbtn=new GridBagLayout();
		gbl_logbtn.rowHeights = new int[] {70, 70, 70};
		gbl_logbtn.rowWeights = new double[]{1.0, 1.0, 1.0};
		gbl_logbtn.columnWeights = new double[]{0.0, 1.0};
		gbl_logbtn.columnWidths=new int[] {100, 500};
		logButton.setLayout(gbl_logbtn);
		
		startButton=new JButton("Start");
		GridBagConstraints gbc_start=new GridBagConstraints();
		gbc_start.anchor=GridBagConstraints.CENTER;
		gbc_start.fill=GridBagConstraints.BOTH;
		gbc_start.gridx=0;
		gbc_start.gridy=0;
		gbc_start.insets=new Insets(0, 5, 5, 5);
		startButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				runIPeak=checkParams();
				//System.out.println(runIPeak);
				if(runIPeak!=null){
					//System.out.println("Yes!");
					runIPeak.setLog(log);
					startButton.setEnabled(false);
					advancedOptions.setEnabled(false);
					runIPeak.execute();
				}else{
					//System.out.println("No!");
					startButton.setEnabled(true);
					advancedOptions.setEnabled(true);
				}
			}
		});
		logButton.add(startButton, gbc_start);
		
		JScrollPane jsp = new JScrollPane();
		jsp.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		jsp.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		GridBagConstraints gbc_jsp = new GridBagConstraints();
		gbc_jsp.weighty = 1.0;
		gbc_jsp.gridheight = 3;
		gbc_jsp.insets = new Insets(0, 0, 5, 5);
		gbc_jsp.fill = GridBagConstraints.BOTH;
		gbc_jsp.gridx = 1;
		gbc_jsp.gridy = 0;
		logButton.add(jsp, gbc_jsp);
		
		log = new JTextArea();
		log.setEditable(false);
		log.setTabSize(4);
		jsp.setViewportView(log);
		
		stopButton=new JButton("Stop");
		GridBagConstraints gbc_stop=new GridBagConstraints();
		gbc_stop.anchor=GridBagConstraints.CENTER;
		gbc_stop.fill=GridBagConstraints.BOTH;
		gbc_stop.gridx=0;
		gbc_stop.gridy=1;
		gbc_stop.insets=new Insets(0, 5, 5, 5);
		stopButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				startButton.setEnabled(true);
				advancedOptions.setEnabled(true);
				if(runIPeak!=null){
					runIPeak.cancel(true);
					runIPeak=null;
				}
			}
		});
		logButton.add(stopButton, gbc_stop);
		
		GridBagConstraints gbc_logButton = new GridBagConstraints();
		gbc_logButton.weighty = 1.0;
		gbc_logButton.anchor = GridBagConstraints.NORTH;
		gbc_logButton.fill = GridBagConstraints.HORIZONTAL;
		gbc_logButton.gridx = 0;
		gbc_logButton.gridy = 2;
		getContentPane().add(logButton, gbc_logButton);
	}
	
	

	private RunIPeak checkParams() {
		RunIPeak ipeak=new RunIPeak();
		if(!msgfText.getText().isEmpty() && msgfText.getText()!=null && new File(msgfText.getText()).exists())
			ipeak.setMsgf(msgfText.getText());
		else{
			ipeak=null;
			JOptionPane.showMessageDialog(null, "Missing MSGF+ Results!", "ERROR", JOptionPane.ERROR_MESSAGE);
			return ipeak;
		}
		if(!omssaText.getText().isEmpty() && omssaText.getText()!=null && new File(omssaText.getText()).exists())
			ipeak.setOmssa(omssaText.getText());
		else{
			ipeak=null;
			JOptionPane.showMessageDialog(null, "Missing OMSSA Results!", "ERROR", JOptionPane.ERROR_MESSAGE);
			return ipeak;
		}
		if(!tandemText.getText().isEmpty() && tandemText.getText()!=null && new File(tandemText.getText()).exists())
			ipeak.setXtandem(tandemText.getText());
		else{
			ipeak=null;
			JOptionPane.showMessageDialog(null, "Missing X!Tandem Results!", "ERROR", JOptionPane.ERROR_MESSAGE);
			return ipeak;
		}
		if(!fastaText.getText().isEmpty() && fastaText.getText()!=null 
				&& new File(fastaText.getText()).exists())
			ipeak.setDatabase(fastaText.getText());
		else{
			ipeak=null;
			JOptionPane.showMessageDialog(null, "Missing Database File!", "ERROR", JOptionPane.ERROR_MESSAGE);
			return ipeak;
		}
		if(!outputText.getText().isEmpty() && outputText.getText()!=null 
				&& new File(outputText.getText()).exists())
			ipeak.setOut_file_dir(outputText.getText());
		else{
			ipeak.setOut_file_dir(new File(".").getAbsolutePath());
			//System.out.println(ipeak.getOut_file_dir());
			JOptionPane.showMessageDialog(null, "Output Directory is Missing!\r\n" +
					"Use current direcotry instead!", "Warning", JOptionPane.WARNING_MESSAGE);
		}
		if(!FileRegexText.getText().isEmpty() && FileRegexText.getText()!=null)
			ipeak.setSpectrum_file_regex(FileRegexText.getText());
		else{
			String s="iPeak_Result_"+String.format("%tF", new Date());
			ipeak.setSpectrum_file_regex(s);
			JOptionPane.showMessageDialog(null, "Result File Name is Missing!\r\nUse "+s+" instead!"
							, "Warning", JOptionPane.WARNING_MESSAGE);
		}
		if(!decoyRegexText.getText().isEmpty() && decoyRegexText.getText()!=null)
			ipeak.setDecoyregrex(decoyRegexText.getText());
		else{
			ipeak=null;
			JOptionPane.showMessageDialog(null, "Missing Decoy Regular Expression!", 
					"ERROR", JOptionPane.ERROR_MESSAGE);
			return ipeak;
		}
                
		//set the path of percolator
		ipeak.setPercolator_path(bgi.ipeak.util.Properties.getPercolator_path());

		ipeak.setMod(bgi.ipeak.util.Properties.getModFile_path());
		ipeak.setUmod(bgi.ipeak.util.Properties.getUserModFile_path());
		ipeak.setDebug("false");
		ipeak.setAccessionSplitRegex(delim);
		if(fdr_level.equals("PeptideLevel")){
			ipeak.setUse_peptidefdr("true");
		}
		else{
			ipeak.setUse_peptidefdr("false");
		}
		if(maxFdrText.getText()!=null && !maxFdrText.getText().isEmpty()){
			Double fdr=Double.parseDouble(maxFdrText.getText());
			if(fdr>1){
				JOptionPane.showMessageDialog(null, "FDR is larger than 1", "ERROR", JOptionPane.ERROR_MESSAGE);
				ipeak=null;
				return ipeak;
			}
			ipeak.setMaxfdr(fdr);
		}else{
			ipeak.setMaxfdr(0.01);
			JOptionPane.showMessageDialog(null, "FDR is missing! Use 0.01 instead!",
					"Warning", JOptionPane.WARNING_MESSAGE);
		}
		ipeak.setBc(bc);
//		System.out.println(ipeak.getAccessionSplitRegex()+"\tssss"+delim+"\t"+ipeak.getOut_file_dir()+"\t"+outputText.getText()+"\t"+new File(outputText.getText()).exists());
		return ipeak;
	}

	private String chooseFile(){
		int returnVal = fc.showOpenDialog(IPeakGui.this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            fc.setCurrentDirectory(file);
            return file.getAbsolutePath();
		}
		fc.setSelectedFile(null);
		return null;
	}

	public static JButton getStopButton() {
		return stopButton;
	}

	public static JButton getStartButton() {
		return startButton;
	}
}
