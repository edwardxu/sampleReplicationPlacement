package system;

import java.util.List;
import simulation.Parameters;
import utils.RanNum;

public class Sample extends Dataset {
	
	private double error;
	private double volume;
	private Dataset parentDataset;
	private double lifeCycle;//how many time slots that this sample can live in the system
	private DataCenter toBePlaced = null;//the destination where this sample will finally be placed

	public Sample(double ID, String name, List<DataCenter> dcList) {
		super(ID, name, dcList);
		// TODO Auto-generated constructor stub
		int choice = RanNum.getRandomIntRange(3,1);
		if(choice==1){
			this.error = Parameters.errorLow;
			this.volume = parentDataset.getVolume() * (1 - this.error);
		}
		if(choice==2){
			this.error = Parameters.errorMedium;
			this.volume = parentDataset.getVolume() * (1 - this.error);
		}
		if(choice==3){
			this.error = Parameters.errorHigh;
			this.volume = parentDataset.getVolume() * (1 - this.error);
		}
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
