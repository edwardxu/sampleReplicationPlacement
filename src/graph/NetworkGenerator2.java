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
//
//generate a graph randoml
package graph;

//generate a graph randomly

import java.util.ArrayList;
import java.util.Map;
import java.util.Random;

import org.jgrapht.Graph;
import org.jgrapht.VertexFactory;
import org.jgrapht.generate.GraphGenerator;

import system.DataCenter;
import graph.Node;
import graph.NodeFactory;
import graph.NodeInitialParameters;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;

public class NetworkGenerator2<V, E> implements GraphGenerator<V, E, V> {

	private int generateType = 0; //0: generate data center networks;
									//1: generate virtual machine networks
									//2: copy
	
	private ArrayList<V> vmList = null;
	
	private int numOfDataCenters = Parameters.numOfDataCenters;
	
	public NetworkGenerator2() {
	}

	public void generateGraph(Graph<V, E> target, VertexFactory<V> vertexFactory, Map<String, V> resultMap) {
		Random ran = new Random();
		if (0 == this.getGenerateType()) {
			ArrayList<V> dataCenterList = new ArrayList<V>();
			// Add all the vertices to the set
			// Add data center nodes
			//int numOfDataCenters = ran.nextInt(Parameters.maxNumOfDataCenters - Parameters.minNumOfDataCenters) + Parameters.minNumOfDataCenters;
			
			for (int i = 0; i < numOfDataCenters; i++) {
				NodeInitialParameters ni = new NodeInitialParameters();				
				ni.id = SamplePlacementSimulator.idAllocator.nextId();
				ni.name = "DC-" + ni.id;
				
				NodeFactory<Node> nf = new NodeFactory<Node>(DataCenter.class, ni);
				V newVertex = (V) nf.createVertex();
				dataCenterList.add(newVertex);
				target.addVertex(newVertex);
			}
			
			for (int i = 0; i < dataCenterList.size(); i ++){
				V dataC1 = dataCenterList.get(i);
				for (int j = i; j < dataCenterList.size(); j ++){
					if (j == i) continue;
					if (RanNum.getZeroOneRandomNumbers(Parameters.connectivityProbality)) {
						V dataC2 = dataCenterList.get(j);
						target.addEdge(dataC1, dataC2);
					}
				}
			}
		}
	}

	public int getGenerateType() {
		return generateType;
	}

	public void setGenerateType(int generateType) {
		this.generateType = generateType;
	}

	public ArrayList<V> getVmList() {
		return vmList;
	}

	public void setVmList(ArrayList<V> vmList) {
		this.vmList = vmList;
	}

	public int getNumOfDataCenters() {
		return numOfDataCenters;
	}

	public void setNumOfDataCenters(int numOfDataCenters) {
		this.numOfDataCenters = numOfDataCenters;
	}

}

