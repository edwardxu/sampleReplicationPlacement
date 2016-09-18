package flow;

import graph.Node;
import system.Dataset;
import system.Query;

public class DemandNode extends Node {
	
	private Query query = null; 
	private Dataset dataset = null; 
	
	public DemandNode(double id, String name, Query query, Dataset dataset){
		super(id, name);
		this.setQuery(query); 
		this.setDataset(dataset); 
	}

	public Query getQuery() {
		return query;
	}

	public void setQuery(Query query) {
		this.query = query;
	}

	public Dataset getDataset() {
		return dataset;
	}

	public void setDataset(Dataset dataset) {
		this.dataset = dataset;
	}
}
