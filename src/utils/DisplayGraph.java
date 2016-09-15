package utils;

import system.DataCenter;
import graph.Node;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Map;
import javax.swing.JApplet;
import org.jgraph.JGraph;
import org.jgraph.graph.DefaultGraphCell;
import org.jgraph.graph.GraphConstants;
import org.jgrapht.ListenableGraph;
import org.jgrapht.ext.JGraphModelAdapter;
import org.jgrapht.graph.ListenableDirectedGraph;
import org.jgrapht.graph.ListenableUndirectedGraph;

import system.DataCenter;
import graph.InternetLink;

/**
 * A demo applet that shows how to use JGraph to visualize JGraphT graphs.
 *
 * @author Barak Naveh
 *
 * @since Aug 3, 2003
 */
public class DisplayGraph extends JApplet {
    /**
	 * 
	 */
	private static final long serialVersionUID = -2335680456809438001L;
	
	private static final Color     DEFAULT_BG_COLOR = Color.decode( "#FAFBFF" );
    private static final Dimension DEFAULT_SIZE = new Dimension( 2000, 2000 );

    // 
    private JGraphModelAdapter<Node, InternetLink> m_jgAdapter;
    
    private ListenableGraph<Node, InternetLink> graph;

    public DisplayGraph(ListenableGraph<Node, InternetLink> graph){
    	this.graph = graph;
    }
    /**
     * @see java.applet.Applet#init().
     */
    public void init() {
    	
    	//ListenableGraph g = new ListenableDirectedGraph<Node, MinCostFlowEdge>( MinCostFlowEdge.class );
    	ListenableGraph<Node, InternetLink> g = new ListenableUndirectedGraph<Node, InternetLink>( InternetLink.class );
    	
        // create a visualization using JGraph, via an adapter
        m_jgAdapter = new JGraphModelAdapter<Node, InternetLink>( g );

        JGraph jgraph = new JGraph( m_jgAdapter );

        adjustDisplaySettings( jgraph );
        getContentPane(  ).add( jgraph );
        resize( DEFAULT_SIZE );

        // add vertex
        for (Node node : this.graph.vertexSet()){
        	g.addVertex(node);
        }
        for (InternetLink edge : this.graph.edgeSet()){
        	InternetLink ed = g.addEdge((Node)edge.getEdgeSource(), (Node)edge.getEdgeTarget());
        }
        
        // position vertices nicely within JGraph component
        Toolkit kit = Toolkit.getDefaultToolkit();
        Dimension screenSize = kit.getScreenSize();
        int width = screenSize.width;
        int gridOneSideLength = 30;// 25 actually, for good look, we set to 30;
        int xStep = width / gridOneSideLength;
        
          
        for (Object nd : g.vertexSet()){
        	Node node = (Node)nd;
        	if (node instanceof DataCenter){
        		DataCenter sn = (DataCenter) node;
        		//positionVertexAt(sn, (sn.getGridX() + 1) * xStep, (sn.getGridY() + 1) * xStep);
        	}
        }
    }

    private void adjustDisplaySettings( JGraph jg ) {
        jg.setPreferredSize( DEFAULT_SIZE );

        Color c = DEFAULT_BG_COLOR;
        String colorStr = null;

        try {
            colorStr = getParameter( "bgcolor" );
        }
        catch( Exception e ) {}

        if( colorStr != null ) {
            c = Color.decode( colorStr );
        }
        
        
        jg.setBackground( c );
    }

//    private void positionVertexAt( Object vertex, int x, int y ) {
//        DefaultGraphCell cell = m_jgAdapter.getVertexCell( vertex );
//        Map attr = cell.getAttributes();
//        Rectangle2D b = GraphConstants.getBounds(attr);
//
//        GraphConstants.setBounds( attr, new Rectangle( x, y, (int)b.getWidth(), (int)b.getHeight() ) );
//        
//        Map cellAttr = new HashMap();
//        cellAttr.put(cell, attr);
//        m_jgAdapter.edit(cellAttr, null, null, null);
//    }
}