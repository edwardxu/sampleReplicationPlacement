package graph;

import flow.MinCostFlowEdge;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.ListenableDirectedGraph;
import org.jgrapht.graph.SimpleWeightedGraph;

public class CheckConnectivity {
	
	private SimpleWeightedGraph<Node, InternetLink> network; // our network

	public CheckConnectivity(SimpleWeightedGraph<Node, InternetLink> simpleWeightedGraph) {
		if (simpleWeightedGraph == null) {
			throw new NullPointerException("network is null");
		}
		this.network = simpleWeightedGraph;
	
		
		checkResult(simpleWeightedGraph);
	}
	
	public void checkResult(SimpleWeightedGraph<Node, InternetLink> simpleWeightedGraph){
		ConnectivityInspector g  = new ConnectivityInspector(simpleWeightedGraph);
		if(g.isGraphConnected())
			System.out.println("this network is connectivity graph");
		else
			System.out.println("this network is not a connectivity graph");
	}
	


}
