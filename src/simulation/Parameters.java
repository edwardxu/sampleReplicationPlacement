package simulation;

public class Parameters {
	final public static int printScale = 10;
	
	final public static int numOfTSs = 101;

	final public static int roundNum = 5;
	/**************data center****************/
	final public static double connectivityProbality = 0.2;// the connection probability between edges
	public static int numOfDataCenters = 20;//the number of data centers in the distributed cloud 
	
	//Data centre resource capacity ranges
	final public static double maxComputingPerDC = 2000;//GHz  //the maximum computing capacity of a data center
	final public static double minComputingPerDC = 1000;//GHz  //the minimum computing capacity of a data center

	final public static double maxBandwidthPerLink = 10;//1;//Gbps //the maximum network capacity of an inter-datacentre link
	final public static double minBandwidthPerLink = 1;//0.1;//Gbps   //the minimum network capacity of an inter-datacenter link
	
	//the size of a dataset
	final public static double sizePerDatasetMax = 8;//GB
	final public static double sizePerDatasetMin = 4;//GB
	
	final public static int numOfDatasetPerQueryMax = 3;
	final public static int numOfDatasetPerQueryMin = 1;

	
	final public static double errorLow = 0.1;
	final public static double errorMedium = 0.2;
	final public static double errorHigh = 0.3;
	
	final public static int lifeCycleMax = (int) ((int) numOfTSs * 0.2); 
	final public static int lifeCycleMin = (int) ((int) numOfTSs * 0.02);
	
	/****************generaral cost settings ***********/
	// TODO please double check these cost settings. 
	final public static double storageCostUnitDataMax = 0.0035;// $ per GB data //c_s(v_i) 
	final public static double storageCostUnitDataMin = 0.0010;// $ per GB data //c_s(v_i) 
	
	final public static double processCostUnitDataMax = 0.22;// $ per GB data // c_p(v_i) 
	final public static double processCostUnitDataMin = 0.15;// $ per GB data // c_p(v_i)
	
	final public static double computingAllocatedToUnitData = 10;//GHz $ per GB //r_c
	
	final public static double bandwidthAllocatedToUnitData = 0.065;//Gbps $ per GB //r_t
	
	final public static double maxBandwidthCost = 0.12;// the maximum cost of transferring 1 GB data 
	final public static double minBandwidthCost = 0.05;// the minimum cost of transferring 1 GB data
	
	/*****************settings of users*************/
	final public static int numOfQueriesPerUserMax = 5;
	final public static int numOfQueriesPerUserMin = 1;
	
	
	
//	final public static int numOfGroupsMax = 10; // the number of groups in the cloud 
//	final public static int numOfGroupsMin = 5;
//	
//	public static int numOfUsersInAGroupMax = 20; //the number of users in a group 
//	final public static int numOfUsersInAGroupMin = 5;
//	
//	public static double percentUsersInAGroupMax = 0.3;// maximum percent of users in a group that has data to place into the cloud 
//	final public static double percentUsersInAGroupMin = 0.1;
//	
//	final public static int numOfCandidateDCsAUserMax = 5;
//	final public static int numOfCandidateDCsAUserMin = 1;

}
