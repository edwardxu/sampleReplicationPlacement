package utils;

import java.util.ArrayList;

public class ArrayOperations {
public static void scaleToSpecifiedTotalValue (ArrayList<Integer> array, int totalValue){
		
		int currTotalValue = ArrayOperations.calculateTotalValueInteger(array);
		double ratio = totalValue / (double) currTotalValue;
		for (int i = 0; i < array.size(); i ++){
			array.set(i, (int) Math.floor(array.get(i) * ratio));
		}
		// minor adjustment.
		currTotalValue = ArrayOperations.calculateTotalValueInteger(array);
		int diff = totalValue - currTotalValue;
		int ranIndex = RanNum.getRandomIntRange(array.size(), 0);
		array.set(ranIndex, array.get(ranIndex) + diff);
	}
	
	public static void scaleToSpecifiedTotalValue (ArrayList<Double> array, double totalValue){
		double currTotalValue = ArrayOperations.calculateTotalValue(array);
		double ratio = totalValue / currTotalValue;
		for (int i = 0; i < array.size(); i ++){
			array.set(i, array.get(i) * ratio);
		}
	}
	
	public static Integer calculateTotalValueInteger (ArrayList<Integer> array){
		int total = 0; 
		for (int i = 0; i < array.size(); i ++){
			total += array.get(i);
		}
		return total;
	}
	
	public static Double calculateTotalValue (ArrayList<Double> array){
		double total = 0; 
		for (int i = 0; i < array.size(); i ++){
			total += array.get(i);
		}
		return total;
	}

}
