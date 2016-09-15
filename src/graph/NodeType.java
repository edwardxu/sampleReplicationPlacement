package graph;

public enum NodeType {
	DCNODE(1111), DCNODEWITHSOURCEDATA(1112), VIRTUALNODE(1113); 
	private int value;
	
	private NodeType(int value){
		this.value = value;		
	}
}
