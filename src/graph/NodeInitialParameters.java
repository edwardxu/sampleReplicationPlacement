/*
 * Double Auction for Relay Assignment
 * Copyright (C) 2011 Zichuan Xu
 *
 */

package graph;

import java.util.ArrayList;

public class NodeInitialParameters{
	
	// Common uses
	public double id;
	public String name;
	public int gridX;
	public int gridY;
	
	// used in big data management
	public double processingCost; 
	public double storageCost; 
	
	// Used in Cloudlet placement problems
	public int numOfReqs;
	
	public boolean potentialLocation;
	
	// Used in VNE algorithms
	public double occupied;
	public double capacity;
	
	public ArrayList<Double> perDems; // = new ArrayList<Double>(Parameters.numTSsPerPeriod);
	
	// Used in flow-based request distribution algorithms
	public double serviceRate;
	public int numberOfMachines;
	public int availableMachines;
	public double bandwidthCost;
	public ArrayList<Double> elecPrices = new ArrayList<Double>();
	public double energyPerReq ;
	public double energyWarmUpPerReq;
	public double delay;
	
	public double requestRate;
	public double bandwidthRate;
	public int indexInPara; 
	
	public int timezone;
	
}
