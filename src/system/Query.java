package system;

import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import graph.Node;

public class Query extends Node {
	
	private List<Dataset> datasets = null;
	private User user = null; 
	private int rate = 0;
	private double delayRequirement = Double.MAX_VALUE;
	private DataCenter homeDataCenter = null; 
	/*
	 * used to create a virtual query in the approximation algorithm. 
	 * */
	private Query parent = null; 
	private double ILPcost = 0d;
	private double demand_ = 0d;
	private double demand = 1d; 
	
//	private int startTime;
//	private int occupyPeriod;
//	private double sourceDataVolume = -1d;// read from parameters. the total size of source data of evaluating query q

	/** 
	 * Class constructor to construct a query by specifying its depedant data center list and dataset
	 * 
	 * @param dcList		the list of data centers where the query's home datacenter is selected.
	 * @param allDataSets	the list of datasets where the query's required datasets are chosen. 
	 */
	public Query(List<DataCenter> dcList, List<Dataset> allDataSets) {
		
		super(SamplePlacementSimulator.idAllocator.nextId(), "Query");
		
		int numOfDatasets = RanNum.getRandomIntRange(Parameters.numOfDatasetPerQueryMax, Parameters.numOfDatasetPerQueryMin);
		List<Integer> indexOfDatasets = RanNum.getDistinctInts(allDataSets.size() - 1, 0, numOfDatasets);
		this.datasets = new ArrayList<Dataset>(numOfDatasets);		
		
		for (Integer index : indexOfDatasets){
			this.datasets.add(allDataSets.get(index));
		}
		
		this.setRate(RanNum.getRandomIntRange(Parameters.queryRateMax, Parameters.queryRateMin));
		this.setDelayRequirement((RanNum.getRandomDoubleRange(Parameters.queryDelayRequirementMax, Parameters.queryDelayRequirementMin)));
		
		int indexOfDCDatasetLocated = RanNum.getRandomIntRange(dcList.size(), 0);
		this.setHomeDataCenter(dcList.get(indexOfDCDatasetLocated));
	}
	
	/** 
	 * Class constructor to construct a virtual query by specifying its parent query, where each virtual query has a single dataset of its parent's datasets. 
	 * 
	 * @param query		its parent query. 
	 */
	public Query(Query query, Dataset dataset) {
		super(SamplePlacementSimulator.idAllocator.nextId(), "Virtual Query");
		this.setParent(query);
		this.datasets = new ArrayList<Dataset>(1);
		this.getDatasets().add(dataset);
		this.setDelayRequirement(query.getDelayRequirement());
		this.setHomeDataCenter(query.getHomeDataCenter());
		this.setRate(query.getRate());
	}
	
	public static Comparator<Query> QueryILPCostComparator = new Comparator<Query>() {
		public int compare(Query query1, Query query2) {
			Double ilpCost1 = query1.getILPcost();
			Double ilpCost2 = query2.getILPcost();
			//ascending order
			return ilpCost1.compareTo(ilpCost2);
			//descending order
			//return ilpCost2.compareTo(ilpCost1);
		}
	};
	
	
	/***********************setter and getter*******************************/
	public List<Dataset> getDatasets() {
		return datasets;
	}

	public void setDatasets(List<Dataset> datasets) {
		this.datasets = datasets;
	}

	public User getUser() {
		return user;
	}

	public void setUser(User user) {
		this.user = user;
	}

	public int getRate() {
		return rate;
	}

	public void setRate(int rate) {
		this.rate = rate;
	}

	public double getDelayRequirement() {
		return delayRequirement;
	}

	public void setDelayRequirement(double delayRequirement) {
		this.delayRequirement = delayRequirement;
	}

	public DataCenter getHomeDataCenter() {
		return homeDataCenter;
	}

	public void setHomeDataCenter(DataCenter homeDataCenter) {
		this.homeDataCenter = homeDataCenter;
	}

	public Query getParent() {
		return parent;
	}

	public void setParent(Query parent) {
		this.parent = parent;
	}

	public double getILPcost() {
		return ILPcost;
	}

	public void setILPcost(double iLPcost) {
		ILPcost = iLPcost;
	}

	public double getDemand() {
		return demand;
	}

	public void setDemand(double demand) {
		this.demand = demand;
	}

	public double getDemand_() {
		return demand_;
	}

	public void setDemand_(double demand_) {
		this.demand_ = demand_;
	}

}
