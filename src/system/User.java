package system;

import graph.Node;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class User extends Node{
	/***********parameters************/
	private DataCenter homeDataCenter = null;
		
	//<timeslot, sourcedata>
	private HashMap<Integer, SourceData> sourceData = null;
	//private double sourceDataVolume = -1d;
	
	/*********Construction function**************/
	public User(FrontEnd fe, double ID) {
		super(ID, "front end server" + ID);
		this.setFrontEndServer(fe);
		setSourceData(new HashMap<Integer, SourceData>());
		setCandidateDataCenters(new ArrayList<DataCenter>());
	}
	
	public void clearSourceData(){
		this.setSourceData(new HashMap<Integer, SourceData>());
	}

	/**********getter and setter*************/
	public FrontEnd getFrontEndServer() {
		return frontEndServer;
	}

	public void setFrontEndServer(FrontEnd frontEndServer) {
		this.frontEndServer = frontEndServer;
	}

	public HashMap<Integer, SourceData> getSourceData() {
		return sourceData;
	}

	public void setSourceData(HashMap<Integer, SourceData> sourceData) {
		this.sourceData = sourceData;
	}

	public List<DataCenter> getCandidateDataCenters() {
		return candidateDataCenters;
	}

	public void setCandidateDataCenters(List<DataCenter> candidateDataCenters) {
		this.candidateDataCenters = candidateDataCenters;
	}
	
	public void addSourceData(int timeslot, SourceData sd){
		this.sourceData.put(timeslot, sd);
	}

	public DataCenter getHomeDataCenter() {
		return homeDataCenter;
	}

	public void setHomeDataCenter(DataCenter homeDataCenter) {
		this.homeDataCenter = homeDataCenter;
	}

}
