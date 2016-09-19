package system;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import graph.Node;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import utils.RanNum;

public class Dataset {
	private double ID;
	private Node datacenter;
	private double volume;
	private Set<Sample> samples = null; 
	private Set<Sample> placedSamples = null;
	
	public Dataset (List<DataCenter> dcList) {
		this.ID = SamplePlacementSimulator.idAllocator.nextId();
		
		int indexOfDCDatasetLocated = RanNum.getRandomIntRange(dcList.size(), 0);
		this.datacenter = dcList.get(indexOfDCDatasetLocated);
		// set the volume of this dataset
		this.volume = RanNum.getRandomDoubleRange(Parameters.sizePerDatasetMax, Parameters.sizePerDatasetMin);
		this.setSamples(new HashSet<Sample>());
		this.setPlacedSamples(new HashSet<Sample>());
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
	
	public Sample getSample(double error) {
		for (Sample sample : this.getSamples()){
			if (error == sample.getError())
				return sample; 
		}
		return null; 
	}
	
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

	public Set<Sample> getSamples() {
		return samples;
	}

	public void setSamples(Set<Sample> samples) {
		this.samples = samples;
	}

	public Set<Sample> getPlacedSamples() {
		return placedSamples;
	}

	public void setPlacedSamples(Set<Sample> placedSamples) {
		this.placedSamples = placedSamples;
	} 
}
