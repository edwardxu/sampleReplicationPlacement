//the link between data centers
package graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import system.DataTransmission;
import system.User;
import graph.Node;

import org.jgrapht.graph.DefaultWeightedEdge;

import simulation.Parameters;
import utils.RanNum;

public class InternetLink extends DefaultWeightedEdge {
	
	private static final long serialVersionUID = 1L;
	private double capacity;
	private double cost;// for transmitting 1M data. 
	private double delay; 
	private Node edgeSource = null; 
	private Node edgeTarget = null;
	
	//<query, volumes of multiple intermediate result transmissions>
	private ArrayList<DataTransmission> intermediateResults = new ArrayList<DataTransmission>();
	
	//<query, copies of source data that are transferred through this link. 
	private Map<User, ArrayList<DataTransmission>> sourceData = new HashMap<User, ArrayList<DataTransmission>>();
	
	public InternetLink() {
		super();
		this.capacity = RanNum.getRandomDoubleRange(Parameters.maxBandwidthPerLink,  Parameters.minBandwidthPerLink); 
		this.setCost(RanNum.getRandomDoubleRange(Parameters.maxBandwidthCost, Parameters.minBandwidthCost));
		this.setDelay(RanNum.getRandomDoubleRange(Parameters.maxLinkDelay, Parameters.minLinkDelay));
	}
	
	public void clear(){
		//this.capacity = RanNum.getRandomDoubleRange(Parameters.maxBandwidthPerLink,  Parameters.minBandwidthPerLink); 
		intermediateResults = new ArrayList<DataTransmission>();
		sourceData = new HashMap<User, ArrayList<DataTransmission>>();
	}
	
	public void admitIntermediateDataTransmission(
			DataTransmission dt) {
		this.intermediateResults.add(dt);
	}
	
	public void admitSourceDataTransmission(User u, DataTransmission dt){
		if (null == this.sourceData.get(u))
			this.sourceData.put(u, new ArrayList<DataTransmission>());
		this.sourceData.get(u).add(dt);
	}

	public double getAvailableBandwidth() {
		double occupy = 0.0;
		
		for (Entry<User, ArrayList<DataTransmission>> entry : this.sourceData.entrySet()) {
			for (DataTransmission dt : entry.getValue()) {
				occupy += dt.getBandwidth();
			}
		}
		
		for (DataTransmission dt : this.intermediateResults) {
			occupy += dt.getBandwidth();
		}
		
		double availableNetworkRes = this.capacity - occupy;
		return availableNetworkRes;
	}
	
	public double getAvailableBandwidthSourceData() {
		double occupy = 0.0;
		
		for (Entry<User, ArrayList<DataTransmission>> entry : this.sourceData.entrySet()) {
			for (DataTransmission dt : entry.getValue()) {
				occupy += dt.getBandwidth();
			}
		}
		
		double availableNetworkRes = this.capacity - occupy;
		return availableNetworkRes;
	}
	
	public double getAvailableBandwidthIntermediateData() {
		double occupy = 0.0;
		
		for (DataTransmission dt : this.intermediateResults) {
			occupy += dt.getBandwidth();
		}
		
		double availableNetworkRes = this.capacity - occupy;
		return availableNetworkRes;
	}
	
	public double getBandwidthOccupiedBySourceData() {
		double occupy = 0.0;
		
		for (Entry<User, ArrayList<DataTransmission>> entry : this.sourceData.entrySet()) {
			for (DataTransmission dt : entry.getValue()) {
				occupy += dt.getBandwidth();
			}
		}
		
		return occupy;
	}
	
	public double getBandwidthOccupiedByIntermediateData() {
		double occupy = 0.0;
		
		for (DataTransmission dt : this.intermediateResults) {
			occupy += dt.getBandwidth();
		}
		
		return occupy;
	}
	
	
	public void releaseUser(User u){
		this.intermediateResults.remove(u);
		this.sourceData.remove(u);
	}
	

	public void clearUsers(){
		this.intermediateResults = new ArrayList<DataTransmission>();
		this.sourceData = new HashMap<User, ArrayList<DataTransmission>>();
	}
	
	@Override
	public String toString(){
		return "Capacity: " + this.capacity;
	}
	
	@Override 
	public boolean equals(Object another){
		//check for self-comparison
		if (this == another)
			return true;
		//use instanceof instead of getClass here for two reasons
		//1. if need be, it can match any supertype, an dnot just one class;
		//2. it renders an explict check for "that == null" redundant, since it does the check for null already--"null instanceof [type]" always returns false
		if(!(another instanceof InternetLink))
			return false;
		
		Object thisS = this.getEdgeSource();
		Object anotherS = ((InternetLink) another).getEdgeSource();
		
		Object thisT = this.getEdgeTarget();
		Object anotherT = ((InternetLink) another).getEdgeTarget();
		
		//The algorithm only accepts comparison between identical V types
		if((!(thisS instanceof Node)) && (!(anotherS instanceof Node)) && (!(thisT instanceof Node)) && (!(anotherT instanceof Node)))
			return false;
		
		if((((Node )thisS).getID() == ((Node)anotherS).getID()) && ( ((Node)thisT).getID() == ((Node)anotherT).getID() ) )
			return true;
		
		return false;
	} 
	
	
	/**********************setter and getter*****************************/
	
	public void setCapacity(double capacity) {
		this.capacity = capacity;
	}

	public double getCapacity() {
		return capacity;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public double getCost() {
		return cost;
	}

	public void setSourceData(Map<User, ArrayList<DataTransmission>> sourceData) {
		this.sourceData = sourceData;
	}

	public Map<User, ArrayList<DataTransmission>> getSourceData() {
		return sourceData;
	}

	public ArrayList<DataTransmission> getIntermediateResults() {
		return intermediateResults;
	}
	
	public void setIntermediateResults(ArrayList<DataTransmission> intermediateResults) {
		this.intermediateResults = intermediateResults;
	}

	public Node getEdgeSource() {
		return edgeSource;
	}

	public void setEdgeSource(Node edgeSource) {
		this.edgeSource = edgeSource;
	}

	public Node getEdgeTarget() {
		return edgeTarget;
	}

	public void setEdgeTarget(Node edgeTarget) {
		this.edgeTarget = edgeTarget;
	}

	public double getDelay() {
		return delay;
	}

	public void setDelay(double delay) {
		this.delay = delay;
	}
}
