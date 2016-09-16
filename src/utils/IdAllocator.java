package utils;

/**
 * utility class for handing out unique ids
 * A using class wishing to assign unique ids to each of its instances should declare a static member variable:
 * <code>
 *   static IdAllocator idAllocator = new IdAllocator();
 * <code>
 *   
 *   In its constructor it should use something like"
 * <code>
 *   id = idAllocator.nextId();
 * <code>
 * @author u5141193
 *
 */

public class IdAllocator {
	protected long nextId = 1000; //starts with 1000
	
	public IdAllocator(){
		
	}

	public synchronized long nextId(){
		return nextId ++;
	}
}
