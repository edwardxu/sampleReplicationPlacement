package system;

import graph.NodeInitialParameters;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import graph.Node;
import simulation.Parameters;

public class DataCenter extends Node{
	
	/****************static properties*****************/
	private double computingCapacity = -1d;
	private double availComputing = -1d;
	private Set<Sample> admittedSamples = null;
	private double volumeOfAdmittedSamples = 0d;
	private double processingCost = 0d;
	private double storageCost = 0d;
	
//	private Map<Group, Set<InternetLink>> multicastTrees = new HashMap<Group, Set<InternetLink>>();
	
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
	
	// used to construct a virtual data center, only parent is not null, and other properties are all null. 
	public DataCenter(double id, String name, DataCenter parent){
		super(id, name);
		this.parent = parent; 
	}
	
	private DataCenter(double id, String name){
		super(id, name);
		this.setAdmittedSamples(new HashSet<Sample>());
		Random ran = new Random();
		this.computingCapacity = ran.nextDouble()* (Parameters.maxComputingPerDC - Parameters.minComputingPerDC) + Parameters.minComputingPerDC;
		this.availComputing = computingCapacity;
	}

	/*************functions*************/
	
	public void clear(){
		this.getAdmittedSamples().clear();
		this.volumeOfAdmittedSamples = 0d;
		this.availComputing = this.getComputingCapacity();
		this.parent = null;
	}
	public void admitSample(Sample sample){
		this.getAdmittedSamples().add(sample);
		this.volumeOfAdmittedSamples += sample.getVolume();
	}
	
	public void removeSample(Sample sample){
		if (!this.getAdmittedSamples().contains(sample))
			System.out.println("Sample not exist! Removal failure!");
		else {
			this.getAdmittedSamples().remove(sample);
			this.volumeOfAdmittedSamples -= sample.getVolume();
		}
	}
	
	public double getAvailableComputing(){
		double occupiedComputing = this.volumeOfAdmittedSamples * Parameters.computingAllocatedToUnitData;
		this.availComputing = this.computingCapacity - occupiedComputing;
		return this.availComputing;
	}
	public void clearAdmittedSamples(){
		this.getAdmittedSamples().clear();
		this.volumeOfAdmittedSamples = 0d;
	}
	
	/*************getter and setter*************/
	public double getComputingCapacity() {
		return computingCapacity;
	}
	
	public void setComputingCapacity(double computingCapacity) {
		this.computingCapacity = computingCapacity;
	}
	
	public double getVolumeOfAdmittedSamples() {
		return volumeOfAdmittedSamples;
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

	public Set<Sample> getAdmittedSamples() {
		return admittedSamples;
	}

	public void setAdmittedSamples(Set<Sample> admittedSamples) {
		this.admittedSamples = admittedSamples;
	}
}
