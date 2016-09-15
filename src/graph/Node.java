package graph;

public class Node {
	private double ID;
	private String name = " ";
	private NodeType type;
	
	public Node (double id, String name){
		this.setID(id);
		this.setName(name);	
	}
	
	@Override
	public boolean equals (Object another){
		if (this == another)
			return true;
		
		if (!(another instanceof Node))
			return false;
		
		if (this.ID == ((Node) another).getID())
			return true;
		
		else 
			return false;
	}

	public void setID(double iD) {
		ID = iD;
	}

	public double getID() {
		return ID;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}

	public void setType(NodeType type) {
		this.type = type;
	}

	public NodeType getType() {
		return type;
	}
	
}
