package simulation;

public class Parameters {
	
	/**************overall settings****************/
	final public static int printScale = 10;
	
	final public static int numOfTrials = 5;
	final public static int roundNum = 5;
	
	/**************data center and links****************/
	final public static double connectivityProbality = 0.2;// the connection probability between edges
	public static int numOfDataCenters = 20;//the number of data centers in the distributed cloud 
	
	//Data centre resource capacity ranges
	final public static double maxComputingPerDC = 2000;//GHz  //the maximum computing capacity of a data center
	final public static double minComputingPerDC = 100;//1000 GHz  //the minimum computing capacity of a data center

	final public static double minLinkDelay = 50; //ms
	final public static double maxLinkDelay = 100; //ms 50
	
	final public static double maxBandwidthPerLink = 1000;//1;//Gbps //the maximum network capacity of an inter-datacentre link
	final public static double minBandwidthPerLink = 100;//0.1;//Gbps   //the minimum network capacity of an inter-datacenter link
	
	/**************data sets and their samples****************/
	//the size of a dataset
	final public static double sizePerDatasetMax = 10;//GB
	final public static double sizePerDatasetMin = 5;//GB
	
	final public static double [] errorBounds = {0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};
	
	public static int maxNumOfDatasetsPerTS = 100;// to be reset
	public static int minNumOfDatasetsPerTS = 50;// to be reset
	
	//public static int numOfSamplesEachDataset = 3;
	/**************queries****************/
	public static int maxNumOfQueriesPerTS = 200; //100 to be reset
	public static int minNumOfQueriesPerTS = 50; // to be reset
	
	final public static int queryRateMax = 10;// to be reset
	final public static int queryRateMin = 1;// to be reset
	
	final public static int numOfDatasetPerQueryMax = 3;
	final public static int numOfDatasetPerQueryMin = 1;
	
	final public static double queryDelayRequirementMax = 1000; // ms 2000
	final public static double queryDelayRequirementMin = 500;// ms 300
	
	//final public static int lifeCycleMax = (int) ((int) numOfTrials * 0.2); 
	//final public static int lifeCycleMin = (int) ((int) numOfTrials * 0.02);
	
	/****************general cost settings ***********/
	// TODO please double check these cost settings. 
	final public static double storageCostUnitDataMax = 0.0035;// $ per GB data //c_s(v_i) 
	final public static double storageCostUnitDataMin = 0.0010;// $ per GB data //c_s(v_i) 
	
	final public static double processCostUnitDataMax = 0.1;//0.22 $ per GB data // c_p(v_i) 
	final public static double processCostUnitDataMin = 0.05;//0.15 $ per GB data // c_p(v_i)
	
	final public static double computingAllocatedToUnitData = 50;//10 GHz $ per GB //r_c
	
	final public static double bandwidthAllocatedToUnitData = 0.065;//Gbps $ per GB //r_t
	
	final public static double maxBandwidthCost = 0.12;// the maximum cost of transferring 1 GB data 
	final public static double minBandwidthCost = 0.05;// the minimum cost of transferring 1 GB data

	public static final String LPOutputFile = "./out/LPOutput.txt";
	
}
