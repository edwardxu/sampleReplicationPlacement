package system;

import simulation.Parameters;
import utils.RanNum;

public class Sample {
	
	private double error;
	private double volume;
	private Dataset parentDataset;
	private double lifeCycle;//how many time slots that this sample can live in the system
	private DataCenter toBePlaced = null;//the destination where this sample will finally be placed

	public Sample(Dataset parent) {
		this.parentDataset = parent;
		int choice = RanNum.getRandomIntRange(Parameters.errorBounds.length - 1, 0);		
		this.error = Parameters.errorBounds[choice];
		this.volume = parentDataset.getVolume() * (1 - this.error);
		this.lifeCycle = RanNum.getRandomIntRange(Parameters.lifeCycleMax, Parameters.lifeCycleMin);
	}
	
	//get the computing demand of processing this sample
	public double getComputingDemands(Sample sample){
		double dems = 0d;
		dems = sample.volume * Parameters.computingAllocatedToUnitData;
		return dems; 
	}
	
	//getter and setter
	public double getError() {
		return error;
	}

	public void setError(double error) {
		this.error = error;
	}


	public double getLifeCycle() {
		return lifeCycle;
	}


	public void setLifeCycle(double lifeCycle) {
		this.lifeCycle = lifeCycle;
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
