package flow;

import java.util.List;

public class FlowPath {
	
	private List<MinCostFlowEdge> edges;
	
	private double flow = 0d;
	
	private double ratio = 0d;
	
	FlowPath(List<MinCostFlowEdge> edges){
		this.setEdges(edges);
	}
	
	FlowPath(List<MinCostFlowEdge> edges, double flow){
		this.setEdges(edges);
		this.flow = flow;
	}

	public double getFlow() {
		return flow;
	}

	public void setFlow(double flow) {
		this.flow = flow;
	}
	
	public void scaleFlow(double scaleFactor){
		this.flow = this.flow / scaleFactor;
	}

	public double getRatio() {
		return ratio;
	}

	public void setRatio(double ratio) {
		this.ratio = ratio;
	}

	public List<MinCostFlowEdge> getEdges() {
		return edges;
	}

	public void setEdges(List<MinCostFlowEdge> edges) {
		this.edges = edges;
	}

}
