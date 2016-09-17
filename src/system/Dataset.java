package system;

import java.util.List;

import graph.Node;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;

public class Dataset {
	private double ID;
	private Node datacenter;
	private double volume;
	
	public Dataset (List<DataCenter> dcList) {
		this.ID = SamplePlacementSimulator.idAllocator.nextId();
		
		int indexOfDCDatasetLocated = RanNum.getRandomIntRange(dcList.size(), 0);
		this.datacenter = dcList.get(indexOfDCDatasetLocated);
		// set the volume of this dataset
		this.volume = RanNum.getRandomDoubleRange(Parameters.sizePerDatasetMax, Parameters.sizePerDatasetMin);
	}
	
//	public boolean equals (Object another){
//		if (this == another)
//			return true;
//		
//		if (!(another instanceof Dataset))
//			return false;
		//?????????????????????????whether all the conditions are satisfied
//		if (this.user.equals(((Dataset) another).getUser()) )
//			return true;
//		else 
//			return false;
//	}
	
	public double getID() {
		return ID;
	}
	public void setID(double iD) {
		ID = iD;
	}
	public Node getDatacenter() {
		return datacenter;
	}
	public void setDatacenter(Node datacenter) {
		this.datacenter = datacenter;
	}
	public double getVolume() {
		return volume;
	}
	public void setVolume(double volume) {
		this.volume = volume;
	} 
}
