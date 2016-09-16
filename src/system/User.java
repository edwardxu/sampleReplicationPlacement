package system;

import graph.Node;
import utils.RanNum;
import java.util.HashMap;
import java.util.List;

public class User extends Node{
	/***********parameters************/
	private DataCenter homeDataCenter = null;
	// < timeslot, number of queries during this timeslot >
	private HashMap <Integer, Query> query = null;
	
	/*********Construction function**************/
	public User( double ID, List<DataCenter> dcList) {
		super(ID, "front end server" + ID);
		int indexOfDCDatasetLocated = RanNum.getRandomIntRange(dcList.size(), 0);
		this.homeDataCenter = dcList.get(indexOfDCDatasetLocated);
		setQuery(new HashMap<Integer, Query>());
	}
	
	public void clearQuery(){
		this.setQuery(new HashMap<Integer, Query>());
	}

	/**********getter and setter*************/
	public DataCenter getHomeDataCenter() {
		return homeDataCenter;
	}
	public void setHomeDataCenter(DataCenter homeDataCenter) {
		this.homeDataCenter = homeDataCenter;
	}
	public HashMap<Integer, Query> getQuery() {
		return query;
	}
	public void setQuery(HashMap<Integer, Query> query) {
		this.query = query;
	}
}
