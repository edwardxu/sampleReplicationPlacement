package utils;

public class ReturnPair<T1, T2> {

	private T1 value1;
	private T2 value2;
	
	public ReturnPair (T1 v1, T2 v2){
		this.setValue1(v1);
		this.setValue2(v2);
	}

	public T1 getValue1() {
		return value1;
	}

	public void setValue1(T1 value1) {
		this.value1 = value1;
	}

	public T2 getValue2() {
		return value2;
	}

	public void setValue2(T2 value2) {
		this.value2 = value2;
	}
}
