package system;

import graph.Node;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;

public class Sample extends Node {
	
	private double error;
	private double volume;
	private Dataset parentDataset;
	//private double lifeCycle;//how many time slots that this sample can live in the system
	private DataCenter toBePlaced = null;//the destination where this sample will finally be placed

	public Sample(Dataset parent, int errorIndex) {
		super(SamplePlacementSimulator.idAllocator.nextId(), "Sample");
		this.parentDataset = parent;
		//int choice = RanNum.getRandomIntRange(Parameters.errorBounds.length - 1, 0);		
		this.error = Parameters.errorBounds[errorIndex];
		this.volume = parentDataset.getVolume() * (1 - this.error);
		//this.lifeCycle = RanNum.getRandomIntRange(Parameters.lifeCycleMax, Parameters.lifeCycleMin);
	}
	
	//get the computing demand of processing this sample
	public double getComputingDemands(Sample sample){
		double dems = 0d;
		dems = sample.volume * Parameters.computingAllocatedToUnitData;
		return dems; 
	}
	
	public void reset(){
		this.toBePlaced = null;
	}
	
	//getter and setter
	public double getError() {
		return error;
	}

	public void setError(double error) {
		this.error = error;
	}

	public double getVolume() {
		return volume;
	}

	public void setVolume(double volume) {
		this.volume = volume;
	}

	public Dataset getParentDataset() {
		return parentDataset;
	}

	public void setParentDataset(Dataset parentDataset) {
		this.parentDataset = parentDataset;
	}
	
	public DataCenter getToBePlaced() {
		return toBePlaced;
	}

	public void setToBePlaced(DataCenter toBePlaced) {
		this.toBePlaced = toBePlaced;
	}

}
