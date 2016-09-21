package algs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.SimpleWeightedGraph;

import graph.InternetLink;
import graph.Node;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import system.DataCenter;
import system.Dataset;
import system.Query;
import system.Sample;

public class ProposedApproximationAlg {

	/************* parameters *****************/
	private SamplePlacementSimulator simulator = null;

	private List<Integer> admittedNumOfQueriesPerTS = new ArrayList<Integer>();

	private List<Double> costPerTS = new ArrayList<Double>();
	private List<Double> storageCostPerTS = new ArrayList<Double>();
	private List<Double> updateCostPerTS = new ArrayList<Double>();
	private List<Double> accessCostPerTS = new ArrayList<Double>();
	private List<Double> processCostPerTS = new ArrayList<Double>();
	
	private double optimalLowerBound = -1d; 

	/************* construction function *****************/
	/** 
	 * Class constructor specifying the simulator
	 */
	public ProposedApproximationAlg(SamplePlacementSimulator simulator) {
		this.setSimulator(simulator);
		for (int i = 0; i < Parameters.numOfTSs; i++) {
			this.admittedNumOfQueriesPerTS.add(i, 0);
			this.costPerTS.add(i, 0.0);
			this.storageCostPerTS.add(i, 0.0);
			this.updateCostPerTS.add(i, 0.0);
			this.accessCostPerTS.add(i, 0.0);
			this.processCostPerTS.add(i, 0.0);
		}
	}
	
	public void run(){
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
				this.simulator.modifyCosts();// double check this. 
			
			// the algorithm is divided into steps
			// step (1), solve the ILP
			
		}
	}
	
	/**
	 * Solve the relatexed version of the ILP. 
	 * 
	 * @param dcList			list of data centers
	 * @param sampleList		list of to be placed samples
	 * @param queryList 		list of to be assigned queries
	 * @return 				the query assignment and sample placement result <DataCenter, <Sample, Set of Queries>>
	 * 
	 * */
	public Map<Node, Map<Sample, Set<Query>>> solveLP(ArrayList<DataCenter> dcList, ArrayList<Sample> sampleList, ArrayList<Query> queryList, double error, 
			SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {

		Map<Node, Map<Sample, Set<Query>>> queryAssignment = null;
		
		List<Query> virtualQueries = new ArrayList<Query>();
		Map<Query, Integer> vQueryIndexMapInList = new HashMap<Query, Integer>();
		int tempIndex = 0; 
		for (Query query : queryList) {
			for (Dataset dataset : query.getDatasets()) {
				Query vQuery = new Query(query, dataset);
				virtualQueries.add(vQuery);
				vQueryIndexMapInList.put(vQuery, tempIndex);
				tempIndex ++; 
			}
		}
		
		Map<Integer, Integer> vQueryIndexToSampleIndex = new HashMap<Integer, Integer>();
		for (int j = 0; j < virtualQueries.size(); j ++) {
			Sample sample = virtualQueries.get(j).getDatasets().get(0).getSample(error);
			for (int o = 0; o < sampleList.size(); o ++){
				if (sample.equals(sampleList.get(o))){
					vQueryIndexToSampleIndex.put(j, o);
					break;
				}
			}
		}
				
		/**
		 * from now on, the problem deals with the set of virtual queries. 
		 * */
		double [][] X = new double[virtualQueries.size()][dcList.size()];// [queries][data centers]
		double [][] Y = new double[sampleList.size()][dcList.size()];// [samples][data centers]
		double[] objs = null; 
		try {
			/**
			 * constrains structure: query * datacenters + sample * datacenters
			 * example: location1 for query1, location 1 for query 2, location 2 for query1, location 2 for query 2, 			-------X
			 * 			location 1 for sample 1, location 1 for sample 2, location 2 for sample 1, location 2 for sample 2. 	-------Y
			 * 
			 * indexing rules: query j, datacenter i, sample o
			 * */
			int consSize = virtualQueries.size() * dcList.size() + sampleList.size() * dcList.size();
			LpSolve solver = LpSolve.makeLp(0, consSize);
			
			/**
			 * constraint (1) 
			 * */
			for (int j = 0; j < virtualQueries.size(); j ++) {
				double [] constraint = new double[consSize];
				
				for (int i = 0; i < dcList.size(); i ++) {
					constraint[virtualQueries.size() * i + j] = 1d; 
				}
				solver.addConstraint(constraint, LpSolve.GE, 1.0);
			}
			
			/**
			 * constraint (2) 
			 * */
			for (int j = 0; j < virtualQueries.size(); j ++) {
				for (int i = 0; i < dcList.size(); i ++){
					double [] constraint = new double[consSize];
					
					constraint[i * virtualQueries.size() + j] = 1d; 
					int oIndex = vQueryIndexToSampleIndex.get(j);
					constraint[virtualQueries.size() * dcList.size() + i * sampleList.size() + oIndex] = -1d; 
					solver.addConstraint(constraint, LpSolve.LE, 0.0);
				}
			}
			
			/**
			 * constraint (3) : capacity constraint of each data center
			 * */
			for (int i = 0; i < dcList.size(); i ++) {
				
				double [] constraint = new double[consSize];
				for (int o = 0; o < sampleList.size(); o ++) {
					constraint[virtualQueries.size() * dcList.size() + i * sampleList.size() + o] = sampleList.get(o).getVolume() * Parameters.computingAllocatedToUnitData; 
				}
				solver.addConstraint(constraint, LpSolve.LE, dcList.get(i).getComputingCapacity());
			}
			
			/**
			 * constraint (4) : positive variable constraint for X
			 * */
			for (int i = 0; i < dcList.size(); i ++) {
				for (int j = 0; j < virtualQueries.size(); j ++){
					double [] constraint = new double[consSize];
					constraint[virtualQueries.size() * i + j] = 1d; 
					solver.addConstraint(constraint, LpSolve.GE, 0d);
				}
			}
			
			/**
			 * constraint (5) : positive variable constraint for Y. 
			 * */
			for (int i = 0; i < dcList.size(); i ++) {
				for (int o = 0; o < sampleList.size(); o ++) {
					double [] constraint = new double[consSize];
					constraint[virtualQueries.size() * dcList.size() + i * sampleList.size() + o] = 1d; 
					solver.addConstraint(constraint, LpSolve.GE, 0d);
				}
			}

			// for (int i = 0; i < consSize; i ++)
			// solver.setInt(i, true);

			objs = new double[consSize];

			for (int i = 0; i < dcList.size(); i ++) {
				for (int j = 0; j < virtualQueries.size(); j ++) {
					// query j is assigned to its sample in data center i 
					DataCenter targetDC = dcList.get(i);
					DataCenter homeDC = virtualQueries.get(j).getHomeDataCenter(); 
					DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, homeDC, targetDC);
					double cost = Double.MAX_VALUE;
					for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++){
						if (0 == e ) 
							cost = 0d;
						cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(i));
					}
					objs[virtualQueries.size() * i + j] = cost * virtualQueries.get(j).getDatasets().get(0).getSample(error).getVolume();
				}
			}
			
			for (int i = 0; i < dcList.size(); i ++) {
				for (int o = 0; o < sampleList.size(); o ++) {
					
					// query j is assigned to its sample in data center i 
					DataCenter targetDC = dcList.get(i);
					DataCenter homeDC = (DataCenter) sampleList.get(o).getParentDataset().getDatacenter();
					
					DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, homeDC, targetDC);
					double cost = Double.MAX_VALUE;
					for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++){
						if (0 == e ) 
							cost = 0d;
						cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
					}
					
					objs[virtualQueries.size() * dcList.size() + i * sampleList.size() + o] = (cost * sampleList.get(o).getVolume() + targetDC.getStorageCost() + targetDC.getProcessingCost()); 
				}
			}

			//solver.setOutputfile(Parameters.LPOutputFile);

			solver.setObjFn(objs);
			// solver.setPresolve(1, 50000);
			//if (this.APs.size() > 500)
			//solver.setScaling(1);// ;set_scaling(solver, 1);
			solver.solve();

			this.optimalLowerBound = solver.getObjective();
			// print solution
			System.out.println("Lower bound of optimal cost : " + this.optimalLowerBound);
			double[] vars = solver.getPtrVariables();
			// String varString = "Value of var : ";
			// for (int i = 0; i < vars.length; i++) {
			// varString += vars[i] + " ";
			// }
			// //System.out.println("Number of Vars : " + vars.length);
			// System.out.println(varString);

			// delete the problem and free memory
			
			for (int i = 0; i < dcList.size(); i ++) {
				for (int j = 0; j < virtualQueries.size(); j ++) {
					X[j][i] = vars[i * virtualQueries.size() + j];
				}
			}
			
			for (int i = 0; i < dcList.size(); i ++) {
				for (int o = 0; o < sampleList.size(); o ++) {
					Y[o][i] = vars[virtualQueries.size() * dcList.size() + i * sampleList.size() + o];
				}
			}
			
			solver.deleteLp();

		} catch (LpSolveException e) {
			e.printStackTrace();
		}
		
		for (int j = 0; j < virtualQueries.size(); j ++){
			Query vQuery = virtualQueries.get(j); 
			double ilpCost = 0d; 
			for (int i = 0; i < dcList.size(); i ++){
				ilpCost += X[j][i] * objs[virtualQueries.size() * i + j];
			}
			vQuery.setILPcost(ilpCost);
		}
			
		
		// step (2) dealing with arrays X and Y. 
		ArrayList<Map<Query, ArrayList<Query>>> clusters = new ArrayList<Map<Query, ArrayList<Query>>>(sampleList.size());
		for (int o = 0; o < sampleList.size(); o ++) {
			Sample sample = sampleList.get(o);
			ArrayList<Query> G_o = new ArrayList<Query>();
			
			for (int j = 0; j <  virtualQueries.size(); j ++) {
				Query vQuery = virtualQueries.get(j);
				if (vQuery.getDatasets().get(0).getSample(error).equals(sample)) {
					G_o.add(vQuery);
					vQuery.setDemand_(0d);
					vQuery.setDemand(1d);
					break; 
				}
			}
			// <cluster center, cluster members>
			Map<Query, ArrayList<Query>> clustersG_o = new HashMap<Query, ArrayList<Query>>();
			for (Query vQuery : G_o) {
				clustersG_o.put(vQuery, new ArrayList<Query>());
				clustersG_o.get(vQuery).add(vQuery);
			}
			Collections.sort(G_o, Query.QueryILPCostComparator);
			
			for (int q1I = 0; q1I < G_o.size(); q1I ++ ) {
				Query vQuery1 = G_o.get(q1I);
				
				boolean foundK = false; 
				for (int q2I = q1I; q2I < G_o.size(); q2I ++ ) {
					
					Query vQuery2 = G_o.get(q2I);
					
					DataCenter sourceDC = vQuery1.getHomeDataCenter(); 
					DataCenter targetDC = vQuery2.getHomeDataCenter(); 
					
					if (vQuery2.getDemand_() > 0) {
						DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, sourceDC, targetDC);
						double cost = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++){
							if (0 == e ) 
								cost = 0d;
							cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
						}
					
						double maxILPCost = (vQuery1.getILPcost() > vQuery2.getILPcost())? vQuery1.getILPcost(): vQuery2.getILPcost();
						if (cost < 4 * maxILPCost) {
							foundK = true; 
							vQuery2.setDemand_(vQuery2.getDemand_() + vQuery1.getDemand());
							if (null == clustersG_o.get(vQuery2))
								clustersG_o.put(vQuery2, new ArrayList<Query>());
							
							if (null != clustersG_o.get(vQuery1)) {
								for (Query toBeMovedQuery : clustersG_o.get(vQuery1)){
									clustersG_o.get(vQuery2).add(toBeMovedQuery);
								}
								clustersG_o.remove(vQuery1);
							} else {
								clustersG_o.get(vQuery2).add(vQuery1);
							}
						}
					}
				}
				
				if (!foundK)
					vQuery1.setDemand_(vQuery1.getDemand());
			}
			clusters.add(clustersG_o);
		}
		
		double [][] new_X = new double[virtualQueries.size()][dcList.size()];// [queries][data centers]
		double [][] new_Y = new double[sampleList.size()][dcList.size()];// [samples][data centers]
		
		for (int o = 0; o < sampleList.size(); o ++) {
			Sample sample = sampleList.get(o);
			Map<Query, ArrayList<Query>> clustersG_o = clusters.get(o);
			
			for (Entry<Query, ArrayList<Query>> entry : clustersG_o.entrySet()) {
				
				ArrayList<DataCenter> fractionallyAssignedDCs = new ArrayList<DataCenter>();
				ArrayList<DataCenter> fractionallyAssignedDCs_ = new ArrayList<DataCenter>();
				double totalAssigned = 0d; 
				int j = vQueryIndexMapInList.get(entry.getKey());
				
				Map<DataCenter, Integer> dcToIndex = new HashMap<DataCenter, Integer>();
				for (int i = 0; i < dcList.size(); i ++) {
					DataCenter dc = dcList.get(i);
					
					if (X[j][i] > 0) {
						totalAssigned += X[j][i];
						fractionallyAssignedDCs.add(dc);
						dcToIndex.put(dc, i);
						DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, entry.getKey().getHomeDataCenter(), dc);
						double cost = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e++) {
							if (0 == e ) 
								cost = 0d;
							cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
						}
						cost *= sample.getVolume();
						
						if (cost <= 2 * entry.getKey().getILPcost()) {
							fractionallyAssignedDCs_.add(dc);
						}
					}
				}
				
				int totalToBeAssigned = 0;
				if (totalAssigned > 1d)
					totalToBeAssigned = (int) totalAssigned;
				else
					totalToBeAssigned = 1; 
				
				//double totalDemandGroup = entry.getKey().getDemand_();
				ArrayList<DataCenter> fullyAssignedDCs = new ArrayList<DataCenter>();
				if (totalToBeAssigned < fractionallyAssignedDCs_.size()) {
					for (int mm = 0; mm < totalToBeAssigned; mm ++){
						new_Y[o][dcToIndex.get(fractionallyAssignedDCs_.get(mm))] = 1d;  
						//new_X[j][dcToIndex.get(fractionallyAssignedDCs_.get(mm))] = 1d; 
						fullyAssignedDCs.add(fractionallyAssignedDCs_.get(mm));
					}
				} else {
					for (int mm = 0; mm < fractionallyAssignedDCs_.size(); mm ++) {
						new_Y[o][dcToIndex.get(fractionallyAssignedDCs_.get(mm))] = 1d;
						//new_X[j][dcToIndex.get(fractionallyAssignedDCs_.get(mm))] = 1d; 
						fullyAssignedDCs.add(fractionallyAssignedDCs_.get(mm));
					}
					for (int mm = 0; mm < totalToBeAssigned - fractionallyAssignedDCs_.size(); mm ++){
						new_Y[o][dcToIndex.get(fractionallyAssignedDCs.get(mm))] = 1d;
						//new_X[j][dcToIndex.get(fractionallyAssignedDCs.get(mm))] = 1d; 
						fullyAssignedDCs.add(fractionallyAssignedDCs.get(mm));
					}
				}
				
				ArrayList<Query> queriesToBeAssigned = entry.getValue();
				for (Query toBeAssignedQ : queriesToBeAssigned) {
					
					//if (toBeAssignedQ.equals(entry.getKey())) continue;
					
					double minCost = Double.MAX_VALUE;
					DataCenter minCostDC = null; 
					for (DataCenter dc : fullyAssignedDCs){
						DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, toBeAssignedQ.getHomeDataCenter(), dc);
						double cost = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e++) {
							if (0 == e ) 
								cost = 0d;
							cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
						}
						cost *= sample.getVolume();
						
						if (minCost > cost){
							minCost = cost; 
							minCostDC = dc; 
						}
					}
					new_X[vQueryIndexMapInList.get(toBeAssignedQ)][dcToIndex.get(minCostDC)] = 1d;
				}
			}
		}
		
		// now process with new_X and new_Y. 
		

		return queryAssignment;
	}	
	
	
	private void setEdgeWeightDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		for (InternetLink il : datacenterNetwork.edgeSet()) {
			datacenterNetwork.setEdgeWeight(il, il.getCost());
		}
	}

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
