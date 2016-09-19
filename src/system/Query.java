package system;

import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;
import java.util.ArrayList;
import java.util.List;

import graph.Node;

public class Query extends Node {
	
	private List<Dataset> datasets = null;
	private User user = null; 
	private int rate = 0;
	private double delayRequirement = Double.MAX_VALUE;
	private DataCenter homeDataCenter = null; 
	
//	private int startTime;
//	private int occupyPeriod;
//	private double sourceDataVolume = -1d;// read from parameters. the total size of source data of evaluating query q

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

}
