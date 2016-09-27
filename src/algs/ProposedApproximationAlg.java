package algs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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

	// the number of admitted for different trials, where each trial has a different set of queries. 
	private List<Integer> admittedNumOfQueriesTrials = new ArrayList<Integer>();

	private List<Double> costTrials = new ArrayList<Double>();
	private List<Double> storageCostTrials = new ArrayList<Double>();
	private List<Double> updateCostTrials = new ArrayList<Double>();
	private List<Double> accessCostTrials = new ArrayList<Double>();
	private List<Double> processCostTrials = new ArrayList<Double>();
	
	private double optimalLowerBound = -1d;

	/************* construction function *****************/
	/** 
	 * Class constructor specifying the simulator
	 */
	public ProposedApproximationAlg(SamplePlacementSimulator simulator) {
		this.setSimulator(simulator);
		for (int i = 0; i < Parameters.numOfTrials; i++) {
			this.admittedNumOfQueriesTrials.add(i, 0);
			this.costTrials.add(i, 0.0);
			this.storageCostTrials.add(i, 0.0);
			this.updateCostTrials.add(i, 0.0);
			this.accessCostTrials.add(i, 0.0);
			this.processCostTrials.add(i, 0.0);
		}
	}
	
	public void run() {
		SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = this.simulator.getDatacenterNetwork();
		setEdgeWeightDataCenterNetwork(datacenterNetwork);

		for (int trial = 0; trial < Parameters.numOfTrials; trial++) {
			if (trial > 0) {
				this.simulator.modifyCosts();// double check this. 
				this.resetDataCenterNetwork(datacenterNetwork);
			}
			
			// the algorithm is divided into steps
			// step (1), solve the ILP
			int currErrorIndex = -1; 
			boolean success = false;
			if ((!success)&&(currErrorIndex < Parameters.errorBounds.length )){
				currErrorIndex ++; 
				double error = Parameters.errorBounds[currErrorIndex];
				
				ArrayList<Sample> sampleList = new ArrayList<Sample>();
				ArrayList<Query> queryList = new ArrayList<Query>();
				
				for (Query query : this.simulator.getQueries().get(trial)) {
					queryList.add(query);
					for (Dataset ds : query.getDatasets()) {
						sampleList.add(ds.getSample(error));
					}
				}
				
				// run the approximation algorithm. 
				success = this.approximation(this.simulator.getDataCenters(), sampleList, queryList, error, datacenterNetwork);
			}
			
			double totalStorageCostTrial = 0d; 
			double totalUpdateCostTrial = 0d; 
			double totalAccessCostTrial = 0d;
			double totalProcessCostTrial = 0d; 
			
			if (success) {
				for (DataCenter dc : this.simulator.getDataCenters()) {
					// storage cost for all placed samples
					for (Sample admittedSample : dc.getAdmittedSamples()) {
						totalStorageCostTrial += admittedSample.getVolume() * dc.getStorageCost();
						totalProcessCostTrial += admittedSample.getVolume() * dc.getProcessingCost();
						
						double updateCost = 0d; 
						if (!admittedSample.getParentDataset().getDatacenter().equals(dc)){
							DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, admittedSample.getParentDataset().getDatacenter(), dc);
							updateCost = Double.MAX_VALUE;
							for (int i = 0; i < shortestPath.getPathEdgeList().size(); i ++){
								if (0 == i ) 
									updateCost = 0d;
								updateCost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(i));
							}
						}
						
						if (Double.MAX_VALUE != updateCost)
							totalUpdateCostTrial += updateCost * admittedSample.getVolume(); 
						else 
							System.out.println("ERROR: path should exist!!");
					}
				
					for (Entry<Sample, Set<Query>> entry : dc.getAdmittedQueriesSamples().entrySet()) {
						Sample admittedSample = entry.getKey();
					
						for (Query accessQuery : entry.getValue()) {
							
							double accessCost = 0d;
							if (!accessQuery.getHomeDataCenter().equals(dc)){
								DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, accessQuery.getHomeDataCenter(), dc);
								accessCost = Double.MAX_VALUE;
								for (int i = 0; i < shortestPath.getPathEdgeList().size(); i ++){
									if (0 == i ) 
										accessCost = 0d;
									accessCost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(i));
								}
							}
							
							if (Double.MAX_VALUE != accessCost)
								totalAccessCostTrial += accessCost * admittedSample.getVolume(); 
							else 
								System.out.println("ERROR: path should exist!!");
						}
					}
				}
			}
			
			this.getCostTrials().set(trial, totalAccessCostTrial + totalStorageCostTrial + totalUpdateCostTrial + totalProcessCostTrial);
			this.getAccessCostTrials().set(trial, totalAccessCostTrial);
			this.getStorageCostTrials().set(trial, totalStorageCostTrial);
			this.getUpdateCostTrials().set(trial, totalUpdateCostTrial);
			this.getProcessCostTrials().set(trial, totalProcessCostTrial);	
		}
	}
	
	/**
	 * Solve the relatexed version of the ILP. 
	 * 
	 * @param dcList			list of data centers
	 * @param sampleList		list of to be placed samples
	 * @param queryList 		list of to be assigned queries
	 * @return 					whether the fractional version of the ILP is successful.
	 * 
	 * */
	public boolean approximation(List<DataCenter> dcList, ArrayList<Sample> sampleList, ArrayList<Query> queryList, double error, 
			SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {

		
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
		
		Map<Sample, Integer> sampleIndexMapInList = new HashMap<Sample, Integer>();
		for (int o = 0; o < sampleList.size(); o ++) {
			sampleIndexMapInList.put(sampleList.get(o), o);
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
			 * example: location 1 for query1, location 1 for query 2, location 2 for query1, location 2 for query 2, 			-------X
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
				solver.addConstraint(constraint, LpSolve.GE, 1d);
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
				solver.addConstraint(constraint, LpSolve.LE, dcList.get(i).getAvailableComputing());
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
					double cost = 0d; 
					if (!homeDC.equals(targetDC)) {
						cost = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++){
							if (0 == e ) 
								cost = 0d;
							cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
						}
					}
					
					objs[virtualQueries.size() * i + j] = cost * virtualQueries.get(j).getDatasets().get(0).getSample(error).getVolume();
				}
			}
			
			for (int i = 0; i < dcList.size(); i ++) {
				for (int o = 0; o < sampleList.size(); o ++) {
					
					// query j is assigned to its sample in data center i 
					DataCenter targetDC = dcList.get(i);
					DataCenter homeDC = (DataCenter) sampleList.get(o).getParentDataset().getDatacenter();
					
					double cost = 0; 
					if (!homeDC.equals(targetDC)){
						DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, homeDC, targetDC);
						cost = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++){
							if (0 == e ) 
								cost = 0d;
							cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
						}
					}
					objs[virtualQueries.size() * dcList.size() + i * sampleList.size() + o] = sampleList.get(o).getVolume() * (cost + targetDC.getStorageCost() + targetDC.getProcessingCost()); 
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
			
			if (this.optimalLowerBound == 0d || this.optimalLowerBound == 1e+30) {// ILP failure. 
				return false;
			}
			
			double [] vars = solver.getPtrVariables();
			// String varString = "Value of var : ";
			// for (int i = 0; i < vars.length; i++) {
			// varString += vars[i] + " ";
			// }
			// //System.out.println("Number of Vars : " + vars.length);
			// System.out.println(varString);

			// delete the problem and free memory
			
			for (int j = 0; j < virtualQueries.size(); j ++) {
				for (int i = 0; i < dcList.size(); i ++) {
					X[j][i] = vars[i * virtualQueries.size() + j];
				}
			}
			
			for (int o = 0; o < sampleList.size(); o ++) {
				for (int i = 0; i < dcList.size(); i ++) {
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
		ArrayList<Map<Query, Set<Query>>> clusters = new ArrayList<Map<Query, Set<Query>>>(sampleList.size());
		for (int o = 0; o < sampleList.size(); o ++) {
			Sample sample = sampleList.get(o);
			ArrayList<Query> G_o = new ArrayList<Query>();
			
			for (int j = 0; j <  virtualQueries.size(); j ++) {
				Query vQuery = virtualQueries.get(j);
				if (vQuery.getDatasets().get(0).getSample(error).equals(sample)) {
					G_o.add(vQuery);
					vQuery.setDemand_(0d);
					vQuery.setDemand(1d);
				}
			}
			// <cluster center, cluster members>
			Map<Query, Set<Query>> clustersG_o = new HashMap<Query, Set<Query>>();
			for (Query vQuery : G_o) {
				clustersG_o.put(vQuery, new HashSet<Query>());
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
						double cost = 0d;
						if (!sourceDC.equals(targetDC)) {
							DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, sourceDC, targetDC);
							cost = Double.MAX_VALUE;
							for (int e = 0; e < shortestPath.getPathEdgeList().size(); e ++) {
								if (0 == e ) 
									cost = 0d;
								cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
							}
						}
					
						double maxILPCost = (vQuery1.getILPcost() > vQuery2.getILPcost())? vQuery1.getILPcost(): vQuery2.getILPcost();
						if (cost < 4 * maxILPCost) {
							foundK = true;
							vQuery2.setDemand_(vQuery2.getDemand_() + vQuery1.getDemand());
							if (null == clustersG_o.get(vQuery2))
								clustersG_o.put(vQuery2, new HashSet<Query>());
							
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
				
				if (!foundK) {
					vQuery1.setDemand_(vQuery1.getDemand());
					if (null == clustersG_o.get(vQuery1))
						clustersG_o.put(vQuery1, new HashSet<Query>());
					clustersG_o.get(vQuery1).add(vQuery1);
				}
			}
			clusters.add(clustersG_o);
		}
		
		double [][] new_X = new double[virtualQueries.size()][dcList.size()];// [queries][data centers]
		double [][] new_Y = new double[sampleList.size()][dcList.size()];// [samples][data centers]
		
		for (int o = 0; o < sampleList.size(); o ++) {
			
			Sample sample = sampleList.get(o);
			Map<Query, Set<Query>> clustersG_o = clusters.get(o);
			
			for (Entry<Query, Set<Query>> entry : clustersG_o.entrySet()) {
				
				ArrayList<DataCenter> fractionallyAssignedDCs = new ArrayList<DataCenter>();
				ArrayList<DataCenter> fractionallyAssignedDCs_ = new ArrayList<DataCenter>();
				ArrayList<DataCenter> fractionallyAssignedDCs__ = new ArrayList<DataCenter>();

				double totalAssigned = 0d; 
				int j = vQueryIndexMapInList.get(entry.getKey());
				
				Map<DataCenter, Integer> dcToIndex = new HashMap<DataCenter, Integer>();
				for (int i = 0; i < dcList.size(); i ++) {
					DataCenter dc = dcList.get(i);
					dcToIndex.put(dc, i);
					
					if (X[j][i] > 0) {
						
						if (X[j][i] < 1d) 
							totalAssigned += X[j][i];
						else 
							totalAssigned += 1d;
						
						fractionallyAssignedDCs.add(dc);
						
						double cost = 0d; 
						if(!entry.getKey().getHomeDataCenter().equals(dc)) {
							DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, entry.getKey().getHomeDataCenter(), dc);
							cost = Double.MAX_VALUE;
							for (int e = 0; e < shortestPath.getPathEdgeList().size(); e++) {
								if (0 == e ) 
									cost = 0d;
								cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
							}
							cost *= sample.getVolume();
						}
						
						if (cost <= 2 * entry.getKey().getILPcost()) {
							fractionallyAssignedDCs_.add(dc);
						} else {
							fractionallyAssignedDCs__.add(dc);
						}
					}
				}
				
				if (totalAssigned == 0d)
					return false; 
				
				int totalToBeAssigned = 0;
				if (totalAssigned > 1d)
					totalToBeAssigned = (int) totalAssigned;
				else if (totalAssigned > fractionallyAssignedDCs.size() ){
					System.out.println(totalAssigned + " v.s. " + dcList.size());
					totalToBeAssigned = fractionallyAssignedDCs.size();
				} else 
					totalToBeAssigned = 1; 
				
				//double totalDemandGroup = entry.getKey().getDemand_();
				ArrayList<DataCenter> fullyAssignedDCs = new ArrayList<DataCenter>();
				if (totalToBeAssigned <= fractionallyAssignedDCs_.size()) {
					for (int mm = 0; mm < totalToBeAssigned; mm ++) {
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
					for (int mm = 0; mm < totalToBeAssigned - fractionallyAssignedDCs_.size(); mm ++) {
						new_Y[o][dcToIndex.get(fractionallyAssignedDCs__.get(mm))] = 1d;
						//new_X[j][dcToIndex.get(fractionallyAssignedDCs.get(mm))] = 1d; 
						fullyAssignedDCs.add(fractionallyAssignedDCs__.get(mm));
					}
				}
				
				Set<Query> queriesToBeAssigned = entry.getValue();
				for (Query toBeAssignedQ : queriesToBeAssigned) {
					
					//if (toBeAssignedQ.equals(entry.getKey())) continue;
					
					double minCost = Double.MAX_VALUE;
					DataCenter minCostDC = null; 
					for (DataCenter dc : fullyAssignedDCs) {
						double cost = 0d; 
						if (!toBeAssignedQ.getHomeDataCenter().equals(dc)) {
							DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, toBeAssignedQ.getHomeDataCenter(), dc);
							cost = Double.MAX_VALUE;
							for (int e = 0; e < shortestPath.getPathEdgeList().size(); e++) {
								if (0 == e ) 
									cost = 0d;
								cost += datacenterNetwork.getEdgeWeight(shortestPath.getPathEdgeList().get(e));
							}
							cost *= sample.getVolume();
						}
						
						if (minCost > cost) {
							minCost = cost; 
							minCostDC = dc; 
						}
					}
					new_X[vQueryIndexMapInList.get(toBeAssignedQ)][dcToIndex.get(minCostDC)] = 1d;
				}
			}
		}
		
		// now process with new_X and new_Y. 
		//ArrayList<Query> vQueriesViolatedDelay = new ArrayList<Query>();
		for (int j = 0; j < virtualQueries.size(); j ++) {
			Query vQ = virtualQueries.get(j);
			
			for (int i = 0; i < dcList.size(); i ++) {
				DataCenter targetDC = dcList.get(i);
				if (new_X[j][i] == 1d) {
					// check the delay. 
					double unitDelay = 0d; 
					if (!vQ.getHomeDataCenter().equals(targetDC)) {
						DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, vQ.getHomeDataCenter(), targetDC);
						unitDelay = Double.MAX_VALUE;
						for (int e = 0; e < shortestPath.getPathEdgeList().size(); e++) {
							if (0 == e ) 
								unitDelay = 0d;
							unitDelay += shortestPath.getPathEdgeList().get(e).getDelay();
						}
					}
					
					Sample sample = vQ.getDatasets().get(0).getSample(error);
					double delay = unitDelay * sample.getVolume();
					
					double newError = error; 
					while (delay > vQ.getDelayRequirement() && newError < Parameters.errorBounds[Parameters.errorBounds.length - 1] ){
						for (int eI = 0; eI < Parameters.errorBounds.length - 1; eI ++) {
							if (newError == Parameters.errorBounds[eI]) {
								newError = Parameters.errorBounds[eI + 1];
								break; 
							}
						}
						delay = unitDelay * vQ.getDatasets().get(0).getSample(newError).getVolume();
					}
					
					// replace the current sample with a sample with a higher error bound. 
					if (newError != error) {
						Sample newSample = vQ.getDatasets().get(0).getSample(newError);
						int sampleIndex = sampleIndexMapInList.get(sample);
						//sampleIndexMapInList.remove(sample);
						sampleIndexMapInList.put(newSample, sampleIndex);
						sampleList.set(sampleIndex, newSample);
					}
				}
			}
		}
		
		// now organize the results for sample placement and query assignments.
		for (int i = 0; i < dcList.size(); i ++) {
			DataCenter dc = dcList.get(i);
			for (int j = 0; j < virtualQueries.size(); j ++) {
				if (new_X[j][i] == 1d) {
					Query vQ = virtualQueries.get(j);
					int sampleIndex = vQueryIndexToSampleIndex.get(j);
				
					if (new_Y[sampleIndex][i] != 1d){
						System.out.println("A sample for this query should have been placed!!!");
					} else {
						dc.admitSample(sampleList.get(sampleIndex), vQ.getParent());
					}
				}
			}
		}
		
		return true; 
	}	
	
	private void resetDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> dcNetwork) {
		for (Node node : dcNetwork.vertexSet()) {
			if (node instanceof DataCenter) {
				((DataCenter) node).reset();
			}
		}
		for (InternetLink il : dcNetwork.edgeSet()) {
			il.clear();
		}
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
	

	public List<Integer> getAdmittedNumOfQueriesTrials() {
		return admittedNumOfQueriesTrials;
	}

	public void setAdmittedNumOfQueriesTrials(List<Integer> admittedNumOfQueriesPerTS) {
		this.admittedNumOfQueriesTrials = admittedNumOfQueriesPerTS;
	}

	public List<Double> getCostTrials() {
		return costTrials;
	}

	public void setCostTrials(List<Double> costPerTS) {
		this.costTrials = costPerTS;
	}

	public List<Double> getStorageCostTrials() {
		return storageCostTrials;
	}

	public void setStorageCostTrials(List<Double> storageCostPerTS) {
		this.storageCostTrials = storageCostPerTS;
	}

	public List<Double> getUpdateCostTrials() {
		return updateCostTrials;
	}

	public void setUpdateCostTrials(List<Double> updateCostPerTS) {
		this.updateCostTrials = updateCostPerTS;
	}

	public List<Double> getAccessCostTrials() {
		return accessCostTrials;
	}

	public void setAccessCostTrials(List<Double> accessCostPerTS) {
		this.accessCostTrials = accessCostPerTS;
	}

	public List<Double> getProcessCostTrials() {
		return processCostTrials;
	}

	public void setProcessCostTrials(List<Double> processCostPerTS) {
		this.processCostTrials = processCostPerTS;
	}
}
