package flow;

import graph.Node;
import system.Dataset;
import system.Query;
import system.Sample;

public class DemandNode extends Node {
	
	private Query query = null; 
	private Dataset dataset = null;
	private Sample sample = null; 
	
	public DemandNode(double id, String name, Query query, Dataset dataset, Sample sample){
		super(id, name);
		this.setQuery(query); 
		this.setDataset(dataset); 
		this.setSample(sample);
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

	public Sample getSample() {
		return sample;
	}

	public void setSample(Sample sample) {
		this.sample = sample;
	}
}
