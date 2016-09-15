package flow;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import system.DataCenter;
import graph.InternetLink;
import graph.Node;


public class MinCostFlowEdge extends DefaultWeightedEdge{
		
	/**
	 * 
	 */
	private static final long serialVersionUID = -1L;
	
	private double cost;
	private double capacity;
	
	// potentials 
	private double length;
	private double costLength;	
	
	//private double flows = 0;
	private double flows = 0;
		
	private InternetLink il = null;
	
	private Node edgeSource = null;
	private Node edgeTarget = null;
	
	/**
	 * Default constructor 
	 */
	public MinCostFlowEdge(){
		super();
	}
	
	/**
	 * Retrieves the source of this edge. 
	 * 
	 * @return source of this edge
	 */
	public Object getSource() {
		return super.getSource();
	}

	/**
	 * Retrieves the target of this edge.
	 * 
	 * @return target of this edge
	 */
	public Object getTarget() {
		return super.getTarget();
	}
	
	@Override
	public String toString() {
		//return super.toString();
		return "Capacity: " + this.capacity + "; Length:"+ this.length + "; Cost:" + this.cost + "; Flow:" + this.flows;
	}
	
	@Override
	public boolean equals(Object another) {

		// Check for self-comparison
		if (this == another)
			return true;

		// Use instanceof instead of getClass here for two reasons
		// 1. if need be, it can match any supertype, and not just one class;
		// 2. it renders an explict check for "that == null" redundant, since
		// it does the check for null already - "null instanceof [type]" always
		// returns false.
		if (!(another instanceof MinCostFlowEdge))
			return false;
		
		
		Object thisS = this.getEdgeSource();
		Object anotherS = this.getEdgeSource();//((MinCostFlowEdge) another).getSource();

		Object thisT = this.getEdgeTarget();//this.getTarget();
		Object anotherT = this.getEdgeTarget();//((MinCostFlowEdge) another).getTarget();

		// The algorithm only accepts comparison between identical V types.
		if ((!(thisS instanceof Node)) && (!(anotherS instanceof Node))
				&& (!(thisT instanceof Node)) && (!(anotherT instanceof Node)))
			return false;
		
		
		System.out.println("thisS is : " + (Node)thisS);
		System.out.println("anotherS is : " + (Node)anotherS);
		System.out.println("thisT is : " + (Node)thisT);
		System.out.println("anotherT is : " + (Node)anotherT);
		
		if((((Node)thisS).getID() == ((Node)anotherS).getID())&&(((Node)thisT).getID() == ((Node)anotherT).getID()))
//		if((thisS.equals(anotherS))&&(thisT.equals(anotherT)))
			return true;
		
		return false;
	}

	public double getCost() {
		return cost;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

	public double getCapacity() {
		return capacity;
	}

	public void setCapacity(double capacity) {
		this.capacity = capacity;
	}

	public double getLength() {
		return length;
	}

	public void setLength(double length) {
		this.length = length;
	}

	public double getCostLength() {
		return costLength;
	}

	public void setCostLength(double costLength) {
		this.costLength = costLength;
	}
	
	public void addFlow(double amount){
		this.flows += amount;
	}
	
	public double getFlows(){
		return this.flows;
	}
	
	public void clearFlows(){
		this.flows = 0;
	}
	
	public void scaleFlow(double scaleFactor){
		this.flows = this.flows / scaleFactor;
	}

	public void setInternetLink(InternetLink il) {
		this.il = il;
	}

	public InternetLink getInternetLink() {
		return il;
	}

	public Node getEdgeSource() {
		return edgeSource;
	}

	public void setEdgeSource(Node edgeSource) {
		this.edgeSource = edgeSource;
	}

	public Node getEdgeTarget() {
		return edgeTarget;
	}

	public void setEdgeTarget(Node edgeTarget) {
		this.edgeTarget = edgeTarget;
	}
}
