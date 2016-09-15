package system;

import graph.NodeInitialParameters;
import java.util.Random;
import graph.Node;
import simulation.Parameters;

public class DataCenter extends Node{
	
	/****************static properties*****************/
	private double computingCapacity = -1d;
	private double availComputing = -1d;
	private double admittedSamples = 0d;
	private double processingCost = 0d;
	private double storageCost = 0d;
//	private Set<DataCenter> virtualDataCenters = new HashSet<DataCenter>();// set of virtual data centers.
//	private Map<Group, Set<InternetLink>> multicastTrees = new HashMap<Group, Set<InternetLink>>();
//	
	/****************properties for a virtual data center*****************/
	// used when this data center is a virtual data center. 
	private DataCenter parent = null;
//	private Set<InternetLink> tree = null; 
	
	/***********Initialization functions***********/
	public DataCenter(NodeInitialParameters ni){
		this(ni.id, ni.name);
		this.processingCost = ni.processingCost;
		this.storageCost = ni.storageCost;
	}
	
	// used to contruct a virtual data center
	public DataCenter(double id, String name, DataCenter parent){
		super(id, name);
		this.parent = parent; 
		//this.group = group;
	}
	
	private DataCenter(double id, String name){
		super(id, name);
		Random ran = new Random();
		this.computingCapacity = ran.nextDouble()* (Parameters.maxComputingPerDC - Parameters.minComputingPerDC) + Parameters.minComputingPerDC;
		this.availComputing = computingCapacity;
	}

	/*************functions*************/
	
	public void clear(){
		this.admittedSamples = 0d;
		this.parent = null;
	}
	public void admitSourceData(double dataVolume){
		this.admittedSamples += dataVolume;
	}
	
	public void removeSourceData (double dataVolume){
		this.admittedSamples -= dataVolume;
	}
	
	public double getAvailableComputing(){
		double occupiedComputing = this.admittedSamples * Parameters.computingAllocatedToUnitData;
		this.availComputing = this.computingCapacity - occupiedComputing;
		return this.availComputing;
	}
	public void clearAdmittedSourceData(){
		this.admittedSamples = 0d;
	}
	
	/*************getter and setter*************/
	public double getComputingCapacity() {
		return computingCapacity;
	}
	public void setComputingCapacity(double computingCapacity) {
		this.computingCapacity = computingCapacity;
	}
	public double getAdmittedSourceData() {
		return admittedSamples;
	}
	public void setAdmittedSourceData(double admittedSourceData) {
		this.admittedSamples = admittedSourceData;
	}
	public DataCenter getParent() {
		return parent;
	}
	public void setParent(DataCenter parent) {
		this.parent = parent;
	}
	
	public double getProcessingCost() {
		return processingCost;
	}
	public void setProcessingCost(double processingCost) {
		this.processingCost = processingCost;
	}
	public double getStorageCost() {
		return storageCost;
	}
	public void setStorageCost(double storageCost) {
		this.storageCost = storageCost;
	}
}
