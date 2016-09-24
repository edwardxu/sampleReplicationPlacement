package simulation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.SimpleWeightedGraph;

import algs.ProposedApproximationAlg;
import algs.ProposedHeuristicAlg;
import graph.InternetLink;
import graph.NetworkGenerator;
import graph.Node;
import system.DataCenter;
import system.Dataset;
import system.Query;
import system.Sample;
import utils.IdAllocator;
import utils.RanNum;

public class SamplePlacementSimulator {
	private SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = null;
	public static IdAllocator idAllocator = new IdAllocator();
	
	private List<DataCenter> dataCenterList = new ArrayList<DataCenter>();
	private Map<Integer, List<Dataset>> datasets = new HashMap<Integer, List<Dataset>>();
	private Map<Integer, List<Sample>> samples = new HashMap<Integer, List<Sample>>();
	private Map<Integer, List<Query>> queries = new HashMap<Integer,List<Query>>();

	public SamplePlacementSimulator() {
		
	}

	public static void main(String[] s) {
		//performanceHeuristic();
		performanceApproximation();
	}
	
	public static void performanceHeuristic() {
		
		for(int round = 0; round < Parameters.roundNum; round ++) {
			String networkIndexPostFix = "";
			if (round > 0) 
				networkIndexPostFix = "-" + round;
			
			SamplePlacementSimulator simulator = new SamplePlacementSimulator();
			simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix);//get the data center network (cloud network)			
			simulator.InitializeDatasetsAndSamples();
			simulator.InitializeQueries();
			ProposedHeuristicAlg heuAlg = new ProposedHeuristicAlg(simulator);
			heuAlg.run();
		}
		
	}
	
	public static void performanceApproximation() {
		
		for(int round = 0; round < Parameters.roundNum; round ++) {
			String networkIndexPostFix = "";
			if (round > 0) 
				networkIndexPostFix = "-" + round;
			
			SamplePlacementSimulator simulator = new SamplePlacementSimulator();
			simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix);//get the data center network (cloud network)			
			simulator.InitializeDatasetsAndSamples();
			simulator.InitializeQueries();
			ProposedApproximationAlg approAlg = new ProposedApproximationAlg(simulator);
			approAlg.run();
		}
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
					dc.clearAdmittedSamples();
				} // TODO double check whether no need for front end servers.
			}
			for (InternetLink il : dcNet.edgeSet()) {
				il.clearUsers();
			}

			this.datacenterNetwork = dcNet;
		}
		return datacenterNetwork;
	}
	
	public Map<Integer, List<Dataset>> InitializeDatasetsAndSamples() {
		if (this.datasets.isEmpty()){
			for(int i = 0; i < Parameters.numOfTSs; i ++) {
				int numOfDatasetsPerTS = RanNum.getRandomIntRange(Parameters.maxNumOfDatasetsPerTS, Parameters.minNumOfDatasetsPerTS);
				List<Dataset> dss = new ArrayList<Dataset>();
				List<Sample> sams = new ArrayList<Sample>();
				for(int j = 0; j < numOfDatasetsPerTS; j ++) {
					Dataset ds = new Dataset(this.getDataCenters());
					dss.add(ds);
					for (int s = 0; s < Parameters.errorBounds.length; s ++){
						Sample sample = new Sample(ds, s);
						sams.add(sample);
						ds.getSamples().add(sample);
					}
				}
				this.datasets.put(i, dss);
				this.samples.put(i, sams);
			}
		}
		
		return this.datasets;
	}
	
	public Map<Integer, List<Query>> InitializeQueries() {
		
		//String numOfQueries = "";
		if (this.queries.isEmpty()){// generate queries
			for(int i = 0; i < Parameters.numOfTSs; i ++) {
				
				int numOfQueriesPerTS = RanNum.getRandomIntRange(Parameters.maxNumOfQueriesPerTS, Parameters.minNumOfQueriesPerTS);
				//numOfQueries += numOfQueriesPerTS + " ";
				List<Query> qus = new ArrayList<Query>();
				for(int j = 0; j < numOfQueriesPerTS; j ++){
					Query quer = new Query(this.getDataCenters(), this.getDatasets().get(i));
					qus.add(quer);
				}
				queries.put(i, qus);
			}
			//System.out.println(numOfQueries);
		}
		return this.queries;
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
	
	public List<Sample> samplesOfDataset(int timeslot, Dataset ds){
		
		List<Sample> foundSamples = new ArrayList<Sample>();
		for (Sample sample : this.getSamples().get(timeslot)){
			if (sample.getParentDataset().equals(ds))
				foundSamples.add(sample);
		}
		
		return foundSamples;
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

	public Map<Integer, List<Dataset>> getDatasets() {
		return datasets;
	}

	public void setDatasets(Map<Integer, List<Dataset>> datasets) {
		this.datasets = datasets;
	}

	public Map<Integer, List<Sample>> getSamples() {
		return samples;
	}

	public void setSamples(Map<Integer, List<Sample>> samples) {
		this.samples = samples;
	}

	public Map<Integer, List<Query>> getQueries() {
		return queries;
	}

	public void setQueries(Map<Integer, List<Query>> queries) {
		this.queries = queries;
	}

}
