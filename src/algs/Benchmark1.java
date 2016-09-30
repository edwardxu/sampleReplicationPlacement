package algs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.SimpleWeightedGraph;

import graph.InternetLink;
import graph.Node;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import system.DataCenter;
import system.Dataset;
import system.Query;
import system.Sample;

public class Benchmark1 {
	
	/************* parameters *****************/
	private SamplePlacementSimulator simulator = null;

	// the number of admitted for different trials, where each trial has a different set of queries. 
	private List<Integer> admittedNumOfQueriesTrials = new ArrayList<Integer>();

	private List<Double> costTrials = new ArrayList<Double>();
	private List<Double> storageCostTrials = new ArrayList<Double>();
	private List<Double> updateCostTrials = new ArrayList<Double>();
	private List<Double> accessCostTrials = new ArrayList<Double>();
	private List<Double> processCostTrials = new ArrayList<Double>();

	private List<Double> averageErrorTrials = new ArrayList<Double>();
	
	private double lowestError = Parameters.errorBounds[0];
	
	/************* construction function *****************/
	public Benchmark1(SamplePlacementSimulator simulator) {
		this.simulator = simulator;
		for (int i = 0; i < Parameters.numOfTrials; i++) {
			this.admittedNumOfQueriesTrials.add(i, 0);
			this.costTrials.add(i, 0.0);
			this.storageCostTrials.add(i, 0.0);
			this.updateCostTrials.add(i, 0.0);
			this.accessCostTrials.add(i, 0.0);
			this.processCostTrials.add(i, 0.0);
			this.averageErrorTrials.add(i, 0.0);
		}
	}

	/************* void function *****************/
	public void run() {
		SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = this.simulator.getDatacenterNetwork();
		setEdgeWeightDataCenterNetwork(datacenterNetwork);

		for (int trial = 0; trial < Parameters.numOfTrials; trial++) {
			
			if (trial > 0) {
				this.simulator.modifyCosts();// double check this.
				this.resetDataCenterNetwork(datacenterNetwork);
			}

			List<Query> toBeAssignedQueries = new ArrayList<Query>();
			for (Query query : simulator.getQueries().get(trial)) {				
				toBeAssignedQueries.add(query);
			}
			
			List<Query> rejectedQueries = new ArrayList<Query>();
			int errorIndexForUnAdmitted = 0;
			for (int i = 0; i < Parameters.errorBounds.length; i ++) {
				if (this.getLowestError() == Parameters.errorBounds[i]){
					errorIndexForUnAdmitted = i;
					break;
				}
			}
			
			boolean increaseAdmittedErrors = false;			
			while(!toBeAssignedQueries.isEmpty()) {
				
				int numOfDCsHaveSampleErrorIncreased = 0; 
				if (increaseAdmittedErrors) {
					for (DataCenter dc : this.simulator.getDataCenters()) {
						Set<Sample> adjustedSamples = new HashSet<Sample>();
						Map<Sample, Set<Query>> sampleQueries = new HashMap<Sample, Set<Query>>();
						
						boolean haveSampleErrorsIncreased = false; 
						for (Sample admittedSample : dc.getAdmittedSamples()) {
							
							int newErrorIndex = -1;
							
							for (int i = 0; i < Parameters.errorBounds.length - 1; i ++) {
								if (admittedSample.getError() == Parameters.errorBounds[i]){
									newErrorIndex = i + 1; 
									break; 
								}
							}
							
							if (-1 != newErrorIndex) {
								Sample newSample = admittedSample.getParentDataset().getSample(Parameters.errorBounds[newErrorIndex]);	
								adjustedSamples.add(newSample);
								sampleQueries.put(newSample, dc.getAdmittedQueriesSamples().get(admittedSample));
								haveSampleErrorsIncreased = true; 
							} else {
								adjustedSamples.add(admittedSample);
								sampleQueries.put(admittedSample, dc.getAdmittedQueriesSamples().get(admittedSample));
							}
						}
						
						if (haveSampleErrorsIncreased)
							numOfDCsHaveSampleErrorIncreased ++; 
						
						dc.setAdmittedSamples(adjustedSamples);
						dc.setAdmittedQueriesSamples(sampleQueries);
					}
				}
				
				double error = Parameters.errorBounds[errorIndexForUnAdmitted];
				
				for (Iterator<Query> iter = toBeAssignedQueries.iterator(); iter.hasNext(); ) {
					
					Query toBeAdmittedQuery = iter.next();
					Map<Sample, DataCenter> sampleAdmissions = new HashMap<Sample, DataCenter>();
					for (Dataset ds : toBeAdmittedQuery.getDatasets()) {
						
						Sample toBeAdmittedSample = ds.getSample(error);
						
						DataCenter maxAvailableDC = null; 
						double maxAvail = Double.MIN_VALUE;
						for (DataCenter dc : simulator.getDataCenters()) {
							
							double delay = 0d; 
							if (!toBeAdmittedSample.getParentDataset().getDatacenter().equals(dc)) {
								DijkstraShortestPath<Node, InternetLink> shortestPath = new DijkstraShortestPath<Node, InternetLink>(datacenterNetwork, toBeAdmittedQuery.getHomeDataCenter(), dc);
								delay = Double.MAX_VALUE;
								for (int i = 0; i < shortestPath.getPathEdgeList().size(); i ++){
									if (0 == i ) 
										delay = 0d;
									delay += shortestPath.getPathEdgeList().get(i).getDelay();
								}
							}
							
							delay *= toBeAdmittedSample.getVolume();
							
							double availR = dc.getAvailableComputing();
							
							if (!dc.isSampleAdmitted(toBeAdmittedSample)) {
								if (delay < toBeAdmittedQuery.getDelayRequirement() && availR > maxAvail){
									maxAvail= availR; 
									maxAvailableDC = dc; 
								}
							} else {
								if (delay < toBeAdmittedQuery.getDelayRequirement() && availR > maxAvail && (toBeAdmittedSample.getVolume() * Parameters.computingAllocatedToUnitData < availR)){
									maxAvail= availR; 
									maxAvailableDC = dc; 
								}
							}
						}
						
						if (null != maxAvailableDC) {
							sampleAdmissions.put(toBeAdmittedSample, maxAvailableDC);
							maxAvailableDC.admitSample(toBeAdmittedSample, toBeAdmittedQuery);
						}
					}
					
					if (sampleAdmissions.size() == toBeAdmittedQuery.getDatasets().size()) {
						iter.remove();//  admit this query
					} else {
						for (Entry<Sample, DataCenter> entry : sampleAdmissions.entrySet()) {
							entry.getValue().removeSample(entry.getKey(), toBeAdmittedQuery);
						}
						
						if (increaseAdmittedErrors && (numOfDCsHaveSampleErrorIncreased == 0)){
							rejectedQueries.add(toBeAdmittedQuery);
							iter.remove();
							continue; 
						}
					}
				}
				
				if (!toBeAssignedQueries.isEmpty()) {
					if (errorIndexForUnAdmitted < Parameters.errorBounds.length - 1) {
						errorIndexForUnAdmitted ++;
						increaseAdmittedErrors = false;
					} else 
						increaseAdmittedErrors = true;
				}
			}// end while
			
			double totalStorageCostTrial = 0d;
			double totalUpdateCostTrial = 0d;
			double totalAccessCostTrial = 0d;
			double totalProcessCostTrial = 0d;
			
			for (DataCenter dc : this.simulator.getDataCenters()) {
				// storage cost for all placed samples
				if (dc.getAdmittedSamples().isEmpty())
					continue; 
				
				for (Sample admittedSample : dc.getAdmittedSamples()) {
					totalStorageCostTrial += admittedSample.getVolume() * dc.getStorageCost();
					totalProcessCostTrial += admittedSample.getVolume() * dc.getProcessingCost();

					double updateCost = 0d; 
					if (!admittedSample.getParentDataset().getDatacenter().equals(dc)) {
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
						if (!accessQuery.getHomeDataCenter().equals(dc)) {
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
			
			this.getCostTrials().set(trial, totalAccessCostTrial + totalStorageCostTrial + totalUpdateCostTrial + totalProcessCostTrial);
			this.getAccessCostTrials().set(trial, totalAccessCostTrial);
			this.getStorageCostTrials().set(trial, totalStorageCostTrial);
			this.getUpdateCostTrials().set(trial, totalUpdateCostTrial);
			this.getProcessCostTrials().set(trial, totalProcessCostTrial);
			
			Map<Query, Set<Sample>> queryPlacedSamples = new HashMap<Query, Set<Sample>>();
			for (DataCenter dc : this.simulator.getDataCenters()) {
				for (Entry<Sample, Set<Query>> entry : dc.getAdmittedQueriesSamples().entrySet()){
					for (Query query : entry.getValue()) {
						if (null == queryPlacedSamples.get(query))
							queryPlacedSamples.put(query, new HashSet<Sample>());
						queryPlacedSamples.get(query).add(entry.getKey());					
					}
				}
			}
			
			double averageError = 0d; 
			for (Entry<Query, Set<Sample>> entry : queryPlacedSamples.entrySet()) {
				
				double totalSampleVolume = 0d; 
				for (Sample sam : entry.getValue())
					totalSampleVolume += sam.getVolume();
				
				double errorQ = 0d; 
				for (Sample sam : entry.getValue())
					errorQ += (sam.getError() * sam.getVolume() / totalSampleVolume);
				
				averageError += errorQ;
			}
			averageError /= queryPlacedSamples.size(); 
			
			this.getAverageErrorTrials().set(trial, averageError);
		}// end for trials;
	}
	
	private void setEdgeWeightDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		for (InternetLink il : datacenterNetwork.edgeSet()) {
			datacenterNetwork.setEdgeWeight(il, il.getCost());
		}
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

	/********* setter and getter **************/
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

	public List<Double> getAverageErrorTrials() {
		return averageErrorTrials;
	}

	public void setAverageErrorTrials(List<Double> averageErrorTrials) {
		this.averageErrorTrials = averageErrorTrials;
	}

	public double getLowestError() {
		return lowestError;
	}

	public void setLowestError(double lowestError) {
		this.lowestError = lowestError;
	}

}