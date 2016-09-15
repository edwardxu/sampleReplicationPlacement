/*
 * Double Auction for Relay Assignment
 * Copyright (C) 2011 Zichuan Xu
 *
 */

package utils;

import java.util.Random;

public class RanNumSpecProb {

	/**
	 * This function returns true with a specified probability
	 */
	public static boolean getZeroOneRandomNumbers(double probability){
		Random ran = new Random();
		double rN = ran.nextDouble()*100;
		
		if(rN < probability*100){
			return true;
		}else
			return false;
	}
}
