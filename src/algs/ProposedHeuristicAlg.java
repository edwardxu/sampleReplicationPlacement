package algs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.SimpleWeightedGraph;

import flow.Commodity;
import flow.DemandNode;
import flow.MCMC;
import flow.MinCostFlowEdge;
import graph.InternetLink;
import graph.Node;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import system.DataCenter;
import system.Dataset;
import system.Query;
import system.User;

public class ProposedHeuristicAlg {
	/************* parameters *****************/
	private SamplePlacementSimulator simulator = null;

	private List<Integer> admittedNumOfQueriesPerTS = new ArrayList<Integer>();

	private List<Double> costPerTS = new ArrayList<Double>();
	private List<Double> storageCostPerTS = new ArrayList<Double>();
	private List<Double> updateCostPerTS = new ArrayList<Double>();
	private List<Double> accessCostPerTS = new ArrayList<Double>();
	private List<Double> processCostPerTS = new ArrayList<Double>();

	/************* construction function *****************/
	public ProposedHeuristicAlg(SamplePlacementSimulator simulator) {
		this.simulator = simulator;
		for (int i = 0; i < Parameters.numOfTSs; i++) {
			this.admittedNumOfQueriesPerTS.add(i, 0);
			this.costPerTS.add(i, 0.0);
			this.storageCostPerTS.add(i, 0.0);
			this.updateCostPerTS.add(i, 0.0);
			this.accessCostPerTS.add(i, 0.0);
			this.processCostPerTS.add(i, 0.0);
		}
	}

	/************* void function *****************/
	public void run() {
		SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = this.simulator.getDatacenterNetwork();
		setEdgeWeightDataCenterNetwork(datacenterNetwork);

		int numOfAdmittedQueriesPerTimeSlot = 0;

		String costPerTimeSlot = "";
		String storageCostPerTimeSlot = "";
		String updateCostPerTimeSlot = "";
		String accessCostPerTimeSlot = "";
		String processCostPerTimeSlot = "";

		for (int timeslot = 0; timeslot < Parameters.numOfTSs; timeslot++) {

			if (timeslot > 0)
				this.simulator.modifyCosts();

			List<Commodity> commodities = new ArrayList<Commodity>();
			SimpleWeightedGraph<Node, MinCostFlowEdge> flowNet = this.initializeFlowNetwork(commodities, timeslot);

			for (Commodity comm : commodities){// route commodities one by one, and update flow network weights and capacities
				
				
			}

			resetDataCenterNetwork(datacenterNetwork);
		}
	}

	// functions
	// return the constructed flow network (edges and vertice in datacenter
	// network + the auxiliaty virtual-sink node)
	private SimpleWeightedGraph<Node, MinCostFlowEdge> initializeFlowNetwork(List<Commodity> commodities,
			int timeslot) {

		SimpleWeightedGraph<Node, MinCostFlowEdge> flowNetwork = new SimpleWeightedGraph<Node, MinCostFlowEdge>(
				MinCostFlowEdge.class);
		// create a virtual sink node, and add it to the flownetwork
		Node virtualSink = new Node(SamplePlacementSimulator.idAllocator.nextId(), "Virtual-Sink");
		flowNetwork.addVertex(virtualSink);

		SimpleWeightedGraph<Node, InternetLink> dcNetwork = simulator.getDatacenterNetwork();// construct the flow network.		
		List<DataCenter> dcNodes = simulator.getDataCenters();//get all datacenter nodes of the datacenter network
		List<Query> queriesAccessData = simulator.getQueries().get(timeslot);//get all queries at time slot "timeslot"
		List<DemandNode> demandNodes = new ArrayList<DemandNode>();//create demand nodes and add all demand nodes into the flow network
		
		for(Query query : queriesAccessData) {
			List<Dataset> dsListOfAQuery = query.getDatasets();//get the list of datasets that query will access
			for (Dataset dsQuery : dsListOfAQuery) {
				DemandNode demandNode = new DemandNode(SamplePlacementSimulator.idAllocator.nextId(), "Demand Node", query, dsQuery);
				flowNetwork.addVertex(demandNode);
			}
		}
		//create virtual datacenter nodes and add them into the flow network
		List<DataCenter> virtualDCNodes = new ArrayList<DataCenter>();
		double totalCapacity = 0d; 
		for(DataCenter dc : dcNodes){
			flowNetwork.addVertex(dc);
			DataCenter vDC = new DataCenter(SamplePlacementSimulator.idAllocator.nextId(), "Virtual Data Center", dc);
			virtualDCNodes.add(vDC);
			flowNetwork.addVertex(vDC);
			MinCostFlowEdge edge = flowNetwork.addEdge(dc, vDC);
			
			edge.setCost(0d);
			edge.setCapacity(dc.getComputingCapacity());//TODO double check this. 
			totalCapacity += dc.getComputingCapacity();
			
			MinCostFlowEdge edge2 = flowNetwork.addEdge(vDC, virtualSink);			
			edge2.setCost(0d);
			edge2.setCapacity(dc.getComputingCapacity());//TODO double check this. 
		}
		
		for (Node deNode : demandNodes){
			DemandNode demandNode = (DemandNode) deNode;
			for (Node dcNode : dcNodes){
				//check the delay from the home datacenter of the query to dcNode, if yes, there is an edge
				DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(dcNetwork, demandNode.getDataset().getDatacenter(), dcNode);
				double delay = Double.MAX_VALUE;
				double cost = Double.MAX_VALUE;
				for (int i = 0; i < shortestPath.getPathEdgeList().size(); i ++){
					if (0 == i ) {
						delay = 0d;
						cost = 0d; 
					}
					delay += shortestPath.getPathEdgeList().get(i).getDelay();
					cost += dcNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(i));
				}
				
				if (delay <= demandNode.getQuery().getDelayRequirement()){
					MinCostFlowEdge edge = flowNetwork.addEdge(deNode, dcNode);
					edge.setCost(cost);
					edge.setCapacity(totalCapacity);//TODO double check this. 
				}
			}
		}
		
		return flowNetwork;
	}

	private void setEdgeWeightDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		for (InternetLink il : datacenterNetwork.edgeSet()) {
			datacenterNetwork.setEdgeWeight(il, il.getCost());
		}
	}



	private void resetDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> dcNetwork) {
		for (Node node : dcNetwork.vertexSet()) {
			if (node instanceof DataCenter) {
				((DataCenter) node).clear();
			}
		}
		for (InternetLink il : dcNetwork.edgeSet()) {
			il.clear();
		}
	}

	/********* setter and getter **************/
	public SamplePlacementSimulator getSimulator() {
		return simulator;
	}

	public void setSimulator(SamplePlacementSimulator simulator) {
		this.simulator = simulator;
	}

	public List<Integer> getAdmittedNumOfQueriesPerTS() {
		return admittedNumOfQueriesPerTS;
	}

	public void setAdmittedNumOfQueriesPerTS(List<Integer> admittedNumOfQueriesPerTS) {
		this.admittedNumOfQueriesPerTS = admittedNumOfQueriesPerTS;
	}

	public List<Double> getCostPerTS() {
		return costPerTS;
	}

	public void setCostPerTS(List<Double> costPerTS) {
		this.costPerTS = costPerTS;
	}

	public List<Double> getStorageCostPerTS() {
		return storageCostPerTS;
	}

	public void setStorageCostPerTS(List<Double> storageCostPerTS) {
		this.storageCostPerTS = storageCostPerTS;
	}

	public List<Double> getUpdateCostPerTS() {
		return updateCostPerTS;
	}

	public void setUpdateCostPerTS(List<Double> updateCostPerTS) {
		this.updateCostPerTS = updateCostPerTS;
	}

	public List<Double> getAccessCostPerTS() {
		return accessCostPerTS;
	}

	public void setAccessCostPerTS(List<Double> accessCostPerTS) {
		this.accessCostPerTS = accessCostPerTS;
	}

	public List<Double> getProcessCostPerTS() {
		return processCostPerTS;
	}

	public void setProcessCostPerTS(List<Double> processCostPerTS) {
		this.processCostPerTS = processCostPerTS;
	}

}
