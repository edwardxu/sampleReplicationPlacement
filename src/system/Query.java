package system;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import graph.Node;
public class Query extends Node {
	
	private List<Dataset> dataset = null;
	private User user = null;
//	private int startTime;
//	private int occupyPeriod;
//	private double sourceDataVolume = -1d;// read from parameters. the total size of source data of evaluating query q

	public Query( double id, String name, List<DataCenter> dcList) {
		super(SamplePlacementSimulator.idAllocator.nextId(), "Query");
		
		List<DataCenter> dataLocations = new ArrayList<DataCenter>();
		int numOfDatasetLocations = RanNum.getRandomIntRange(Parameters.numOfDatasetPerQueryMax, Parameters.numOfDatasetPerQueryMin);
		this.dataset = new ArrayList<Dataset>(numOfDatasetLocations);
		
		List<Integer> indexOfDCsHaveRequestedDataset = RanNum.getDistinctInts(dcList.size(), 0, numOfDatasetLocations);
		for (Integer index : indexOfDCsHaveRequestedDataset){
			dataLocations.add(dcList.get(index));
		}
	}
	
	/***********************setter and getter*******************************/
	public List<Dataset> getDataset() {
		return dataset;
	}

	public void setDataset(List<Dataset> dataset) {
		this.dataset = dataset;
	}

	public User getUser() {
		return user;
	}

	public void setUser(User user) {
		this.user = user;
	}

}
