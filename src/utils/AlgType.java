package utils;

public enum AlgType {
	FINEGRAINED(9999), 
	FIXEDMAX(8888),
	FINEGRAINEDWITHOUTFUTURE(7777),
	FIXEDMAXPAGERANK(65975),
	FIXEDMAXPAGERANKBFS(65976),
	ONLINE(11111),
	ONLINEWITHOUTFUTURE(22222),
	ACTIONINCREASING(1111),
	ACTIONDECREASING(2222),
	ACTIONHOLDING(3333),
	STABLENETWORK(4444),
	STARTUPNETWORK(5555),
	FADINGNETWORK(6666);
	
	private int value;
	
	private AlgType (int value){
		this.setValue(value);
	}

	public void setValue(int value) {
		this.value = value;
	}

	public int getValue() {
		return value;
	}
	
}
