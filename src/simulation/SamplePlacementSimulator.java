package simulation;

import java.util.ArrayList;
import java.util.List;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.SimpleWeightedGraph;

import graph.InternetLink;
import graph.NetworkGenerator;
import graph.Node;
import system.DataCenter;
import system.User;
import utils.IdAllocator;
import utils.RanNum;

public class SamplePlacementSimulator {
	private SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = null;
	public static IdAllocator idAllocator = new IdAllocator();
	private List<DataCenter> dataCenterList = new ArrayList<DataCenter>();
	private List<User> users = new ArrayList<User>();

	public SamplePlacementSimulator() {
		this.setDataCenterList(new ArrayList<DataCenter>());
	}

	public static void main(String[] s) {

	}

	/************************************
	 * Initialization functions
	 * 
	 * @return
	 ************************************/
	public SimpleWeightedGraph<Node, InternetLink> InitializeDataCenterNetwork(
			SimpleWeightedGraph<Node, InternetLink> dcNet, String networkIndexPostfix) {

		if (null == dcNet) {

			ConnectivityInspector<Node, InternetLink> connect = null;
			do {
				// Initialize the data center network
				datacenterNetwork = new SimpleWeightedGraph<Node, InternetLink>(InternetLink.class);
				NetworkGenerator<Node, InternetLink> networkGenerator = new NetworkGenerator<Node, InternetLink>();

				networkGenerator.setSize(Parameters.numOfDataCenters);
				networkGenerator.setGenerateType(0);// generate data center
													// networks
				networkGenerator.setNetworkIndexPostFix(networkIndexPostfix);
				networkGenerator.generateGraph(datacenterNetwork, null, null);
				// displayGraph(substrateNetwork);
				connect = new ConnectivityInspector<Node, InternetLink>(datacenterNetwork);
			} while (!connect.isGraphConnected());

			List<DataCenter> dcs = this.getDataCenters();

			for (DataCenter dc1 : dcs) {
				for (DataCenter dc2 : dcs) {
					InternetLink il = datacenterNetwork.getEdge(dc1, dc2);
					if (null != il) {
						il.setEdgeSource(dc1);
						il.setEdgeTarget(dc2);
					}
				}
			}
		} else {// clear some parameters of dcNet.
			// Initialize parameters of sNet.
			// modify.
			for (Node node : dcNet.vertexSet()) {
				if (node instanceof DataCenter) {
					DataCenter dc = (DataCenter) node;
					dc.clearAdmittedSourceData();
				} // TODO double check whether no need for front end servers.
			}
			for (InternetLink il : dcNet.edgeSet()) {
				il.clearUsers();
			}

			this.datacenterNetwork = dcNet;
		}
		return datacenterNetwork;
	}

	public void modifyCosts() {

		for (Node node : this.datacenterNetwork.vertexSet()) {

			if (node instanceof DataCenter) {

				double processingCost = RanNum.getRandomDoubleRange(Parameters.processCostUnitDataMax,
						Parameters.processCostUnitDataMin);
				double storageCost = RanNum.getRandomDoubleRange(Parameters.storageCostUnitDataMax,
						Parameters.storageCostUnitDataMin);

				DataCenter dc = (DataCenter) node;
				dc.setProcessingCost(processingCost);
				dc.setStorageCost(storageCost);
			}
		}

		for (InternetLink il : this.datacenterNetwork.edgeSet()) {
			il.setCost(RanNum.getRandomDoubleRange(Parameters.maxBandwidthCost, Parameters.minBandwidthCost));
		}
	}

	public List<DataCenter> getDataCenters() {
		if (!this.dataCenterList.isEmpty())
			return this.dataCenterList;

		for (Node dc : this.datacenterNetwork.vertexSet()) {
			if (dc instanceof DataCenter)
				this.dataCenterList.add((DataCenter) dc);
		}

		return this.dataCenterList;
	}

	public int getNumberOfDataCenters() {
		int numDC = 0;
		for (Node dc : this.datacenterNetwork.vertexSet()) {
			if (dc instanceof DataCenter)
				numDC++;
		}
		return numDC;
	}

	/******** getter and setter ************/
	public SimpleWeightedGraph<Node, InternetLink> getDatacenterNetwork() {
		return datacenterNetwork;
	}

	public void setDatacenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		this.datacenterNetwork = datacenterNetwork;
	}

	public List<DataCenter> getDataCenterList() {
		return dataCenterList;
	}

	public void setDataCenterList(List<DataCenter> dataCenterList) {
		this.dataCenterList = dataCenterList;
	}

	public List<User> getUsers() {
		return users;
	}

	public void setUsers(List<User> users) {
		this.users = users;
	}

}
