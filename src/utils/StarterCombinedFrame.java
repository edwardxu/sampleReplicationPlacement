package utils;

import java.awt.Frame;
import java.awt.event.AdjustmentEvent;
import java.awt.event.WindowEvent;
import java.applet.Applet;
import java.awt.Event;
import java.awt.Label;
import java.awt.ScrollPane;
import java.awt.Scrollbar;
import java.awt.event.*;

import org.jgrapht.ListenableGraph;

public class StarterCombinedFrame extends Frame implements AdjustmentListener{
		
	public StarterCombinedFrame(String frameTitle, ListenableGraph graph) {
		super(frameTitle);
	    
	    ScrollPane pane = new ScrollPane();  
	    this.add(pane);
	    
		DisplayGraph applet = new DisplayGraph(graph);
		applet.init();
		applet.start();
		pane.add(applet, "Center");
		
		this.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent event){
				dispose();
				System.exit(0);
			}
		});
	}

	@Override
	public void adjustmentValueChanged(AdjustmentEvent e) {
	}
}