/* -------------------------
 * NetworkGenerator.java
 * -------------------------
 *
 * Original Author:  Zichuan Xu.
 *
 * The node class V need to have four parameters: id, x location, y location, rest energy and current energy.
 *
 * Changes
 * -------
 * 2-Dec-2011 : Initial revision (GB);
 *
 */

package graph;

//generate a graph by reading data from external files

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

import org.jgrapht.Graph;
import org.jgrapht.VertexFactory;
import org.jgrapht.generate.GraphGenerator;

import simulation.Parameters;
import system.DataCenter;
import utils.RanNum;
import graph.Node;
import graph.NodeFactory;
import graph.NodeInitialParameters;

public class NetworkGenerator<V, E> implements GraphGenerator<V, E, V> {
	// ~ Instance fields
	// --------------------------------------------------------

	private int size; // the number of data centers

	private int generateType = 0; // 0: generate metro-WLAN networks;
								  // 1: generate virtual networks;
								  // 2: copy 
	
	private boolean specialCase = false; 
		
	private String networkIndexPostFix;

	// ~ Constructors
	// -----------------------------------------------------------

	/**
	 * Construct a new NetworkGenerator.
	 * 
	 * @param size
	 *            number of vertices to be generated
	 * 
	 * @throws IllegalArgumentException
	 *             if the specified size is negative.
	 */
	public NetworkGenerator() {
	}

	// ~ Methods
	// ----------------------------------------------------------------

	/**
	 * {@inheritDoc}
	 */
	public void generateGraph(Graph<V, E> target,
			VertexFactory<V> vertexFactory, Map<String, V> resultMap) {
		
		
		if (0 == this.getGenerateType()) {
			
			ArrayList<V> nodeList = new ArrayList<V>();
			
			String fileName = ".//Data//" + this.size + "-25-25" + this.networkIndexPostFix + ".txt";
			try {
				File file = new File(fileName); 
				BufferedReader reader = new BufferedReader(new FileReader(file)); 
				String lineString = null; 
				int readStatus = -1; // 0: reading vertices data; 1: reading edges data
				int numOfAPRead = 0;
				while ((lineString = reader.readLine()) != null){
					if (lineString.startsWith("#"))
						continue;
					
					if (lineString.contains("VERTICES")){//start to parse vertices data
						readStatus = 0;
						continue;
					} else if (lineString.contains("EDGES")){
						readStatus = 1;
						continue;
					}
					if (0 == readStatus){
						lineString.trim();
						String [] attrs = lineString.split(" ");
						
						NodeInitialParameters ni = new NodeInitialParameters();
						
						ni.id = Integer.parseInt(attrs[0]);
						ni.name = "DataCenter" + ni.id;
						ni.gridX = Integer.parseInt(attrs[2]);
						ni.gridY = Integer.parseInt(attrs[3]);
						ni.processingCost = RanNum.getRandomDoubleRange(Parameters.processCostUnitDataMax, Parameters.processCostUnitDataMin);
						ni.storageCost = RanNum.getRandomDoubleRange(Parameters.storageCostUnitDataMax, Parameters.storageCostUnitDataMin);
						
						NodeFactory<Node> nf = new NodeFactory<Node>(DataCenter.class, ni);
						V newVertex = (V) nf.createVertex();
						nodeList.add(newVertex);
						target.addVertex(newVertex);
					}
					if (1 == readStatus){
						
						lineString.trim();
						String [] attrs = lineString.split(" ");
						
						double fromNodeId = Double.parseDouble(attrs[0]);
						double toNodeId = Double.parseDouble(attrs[1]);
						V fromNode = null;
						V toNode = null;
						
						//randomly generate capacity between Parameters.SNodeMax and Parameters.SNodeMin						
						//LinkInitialParameters li = new LinkInitialParameters();
						
						for (V node : nodeList){
							if (node instanceof Node){
								Node nn = (Node)node;
								
								if ( nn.getID() == fromNodeId){
									fromNode = node;
									if (toNode != null)
										break;
									else 
										continue;
								} 
								if ( nn.getID() == toNodeId){
									toNode = node;
									if (fromNode != null)
										break;
									else 
										continue;
								}
							}
						}
						target.addEdge(fromNode, toNode);
					}
				} 
				reader.close();
			} catch (IOException e) { 
				e.printStackTrace(); 
			}
			
//			// add front end servers. 
//			
//			ArrayList<Integer> feConnectDCsIndex = RanNum.getDistinctInts(nodeList.size(), 0, Parameters.numOfFrontEndServers);
//			
//			for (Integer index : feConnectDCsIndex){
//				
//				NodeInitialParameters ni = new NodeInitialParameters();
//				
//				ni.id = BigDataManageSimulator.idAllocator.nextId();
//				ni.name = "front end server" + ni.id;
//				
//				NodeFactory<Node> nf = new NodeFactory<Node>(FrontEnd.class, ni);
//				V newVertex = (V) nf.createVertex();
//				target.addVertex(newVertex);
//				target.addEdge(newVertex, nodeList.get(index));
//			}
			
		} else if (1 == this.getGenerateType()){
			//TODO if necessary
		}
	}

	public int getGenerateType() {
		return generateType;
	}

	public void setGenerateType(int generateType) {
		this.generateType = generateType;
	}
	
	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public boolean isSpecialCase() {
		return specialCase;
	}

	public void setSpecialCase(boolean specialCase) {
		this.specialCase = specialCase;
	}

	public String getNetworkIndexPostFix() {
		return networkIndexPostFix;
	}

	public void setNetworkIndexPostFix(String networkIndexPostFix) {
		this.networkIndexPostFix = networkIndexPostFix;
	}
}
