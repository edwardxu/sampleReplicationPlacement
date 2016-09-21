package algs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.SimpleWeightedGraph;

import graph.InternetLink;
import graph.Node;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import simulation.Parameters;
import simulation.SamplePlacementSimulator;
import system.DataCenter;
import system.Query;
import system.Sample;

public class ProposedApproximationAlg {

	/************* parameters *****************/
	private SamplePlacementSimulator simulator = null;

	private List<Integer> admittedNumOfQueriesPerTS = new ArrayList<Integer>();

	private List<Double> costPerTS = new ArrayList<Double>();
	private List<Double> storageCostPerTS = new ArrayList<Double>();
	private List<Double> updateCostPerTS = new ArrayList<Double>();
	private List<Double> accessCostPerTS = new ArrayList<Double>();
	private List<Double> processCostPerTS = new ArrayList<Double>();

	/************* construction function *****************/
	/** 
	 * Class constructor specifying the simulator
	 */
	public ProposedApproximationAlg(SamplePlacementSimulator simulator) {
		this.setSimulator(simulator);
		for (int i = 0; i < Parameters.numOfTSs; i++) {
			this.admittedNumOfQueriesPerTS.add(i, 0);
			this.costPerTS.add(i, 0.0);
			this.storageCostPerTS.add(i, 0.0);
			this.updateCostPerTS.add(i, 0.0);
			this.accessCostPerTS.add(i, 0.0);
			this.processCostPerTS.add(i, 0.0);
		}
	}
	
	public void run(){
		SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = this.simulator.getDatacenterNetwork();
		setEdgeWeightDataCenterNetwork(datacenterNetwork);

		int numOfAdmittedQueriesPerTimeSlot = 0;

		String costPerTimeSlot = "";
		String storageCostPerTimeSlot = "";
		String updateCostPerTimeSlot = "";
		String accessCostPerTimeSlot = "";
		String processCostPerTimeSlot = "";

		for (int timeslot = 0; timeslot < Parameters.numOfTSs; timeslot++) {
			if (timeslot > 0)
				this.simulator.modifyCosts();// double check this. 
			
			// the algorithm is divided into steps
			// step (1), solve the ILP
			
		}
	}
	
	/*
	 * Solve the relatexed version of the ILP. 
	 * 
	 * @para dcList			list of data centers
	 * @para sampleList		list of to be placed samples
	 * @para queryList 		list of to be assigned queries
	 * @return 				the query assignment and sample placement result <DataCenter, <Sample, Set of Queries>>
	 * 
	 * */
	public Map<Node, Map<Sample, Set<Query>>> solveLP(ArrayList<DataCenter> dcList, ArrayList<Sample> sampleList, ArrayList<Query> queryList) {

		Map<Node, Map<Sample, Set<Query>>> queryAssignment = null;
		
		try {

			double [][] X = new double[dcList.size()][queryList.size()];// [data centers][queries]
			double [][] Y = new double[sampleList.size()][dcList.size()];// [samples][data centers]
			
			int n = this.APs.size();

			for (int i = 0; i < n; i++)
				this.totalNumReqs += ((AccessPoint) this.APs.get(i)).getExpectedNumberOfUsers();

			int L = this.potentialLocations.size();// the number of potential
													// locations.
			int K = this.cloudlets.size(); // the number of cloudlets.

			int cloudletCap = this.cloudlets.get(0).getCapacity();

			int consSize = L + n * L;

			LpSolve solver = LpSolve.makeLp(0, consSize);

			// group A, constraints (2), (3), and (4).

			// (2) K * L p_{i, l}
			for (int j = 0; j < n; j++) {
				double[] cons_2 = new double[consSize];
				for (int l = 0; l < L; l++) {
					cons_2[L + j * L + l] = 1.0;
				}
				solver.addConstraint(cons_2, LpSolve.EQ, 1.0);
			}

			// (3) 
			for (int l = 0; l < L; l++) {
				for (int j = 0; j < n; j++) {
					double[] cons_3 = new double[consSize];
					cons_3[L + j * L + l] = 1.0;
					cons_3[l] = -1.0;
					solver.addConstraint(cons_3, LpSolve.LE, 0.0);
				}
			}

			// (4)
			for (int l = 0; l < L; l++) {

				double[] cons_4 = new double[consSize];
				cons_4[l] = -cloudletCap;

				for (int j = 0; j < n; j++) {
					int numReqs = this.APs.get(j).getExpectedNumberOfUsers();
					cons_4[L + j * L + l] = -numReqs;
				}
				solver.addConstraint(cons_4, LpSolve.LE, 0.0);
			}

			// (5)

			double[] cons_5 = new double[consSize];
			for (int l = 0; l < L; l++) {
				cons_5[l] = 1.0;
			}
			solver.addConstraint(cons_5, LpSolve.EQ, K);

			// (6, 7) :
			for (int i = 0; i < consSize; i++) {
				double[] cons_6 = new double[consSize];
				cons_6[i] = 1.0;
				solver.addConstraint(cons_6, LpSolve.GE, 0d);
				// solver.addConstraint(cons_6, LpSolve.LE, 1d);
			}

			for (int i = 0; i < consSize; i++) {
				double[] cons_7 = new double[consSize];
				cons_7[i] = 1.0;
				// solver.addConstraint(cons_7, LpSolve.GE, 0d);
				solver.addConstraint(cons_7, LpSolve.LE, 1d);
			}

			// for (int i = 0; i < consSize; i ++)
			// solver.setInt(i, true);

			double[] objs = new double[consSize];

			for (int j = 0; j < n; j++) {
				AccessPoint ap = this.APs.get(j);
				int numReqs = ap.getExpectedNumberOfUsers();
				for (int l = 0; l < L; l++) {

					Node location = this.potentialLocations.get(l);
					double accWeight = 0;
					if (!ap.equals(location))
						accWeight = this.metricGraph.getEdgeWeight(this.metricGraph.getEdge(ap, location));
					objs[L + j * L + l] = (numReqs * accWeight / totalNumReqs);
				}
			}

			solver.setOutputfile(Parameters.LPOutputFile);

			solver.setObjFn(objs);
			// solver.setPresolve(1, 50000);
			//if (this.APs.size() > 500)
			//solver.setScaling(1);// ;set_scaling(solver, 1);
			solver.solve();

			this.optimalLowerBound = solver.getObjective();
			// print solution
			System.out.println("Lower bound of optimal cost : " + this.optimalLowerBound);
			double[] vars = solver.getPtrVariables();
			// String varString = "Value of var : ";
			// for (int i = 0; i < vars.length; i++) {
			// varString += vars[i] + " ";
			// }
			// //System.out.println("Number of Vars : " + vars.length);
			// System.out.println(varString);

			double totalPlacedCloudlets = 0;
			double numIntegralPlacements = 0;
			for (int l = 0; l < L; l++) {
				AccessPoint location = (AccessPoint) this.potentialLocations.get(l);
				if (vars[l] + this.getEpsilon() >= 1d) {
					location.setPlacedCloudletLP(1d);
					totalPlacedCloudlets += location.getPlacedCloudletLP();
					numIntegralPlacements++;
				} else {
					location.setPlacedCloudletLP(vars[l]);
					totalPlacedCloudlets += location.getPlacedCloudletLP();
				}
			}

			if (numIntegralPlacements == K)
				this.setIntegralSolution(true);

			// System.out.println("Total number of cloudlets that are placed by
			// LP :" + totalPlacedCloudlets);

			for (int j = 0; j < n; j++) {
				AccessPoint ap = this.APs.get(j);
				double delay = 0;
				for (int l = 0; l < L; l++) {
					Node location = this.potentialLocations.get(l);
					double unitDelay = 0;
					if (!ap.equals(location))
						unitDelay = this.metricGraph.getEdgeWeight(this.metricGraph.getEdge(ap, location));
					double assignPortion = vars[L + j * L + l];
					if (assignPortion > 0)
						delay += (assignPortion * unitDelay);
				}
				ap.setUnitDelayLP(delay);
			}

			// <Location, <AccessPoint, Number of Reqs>>
			queryAssignment = new HashMap<Node, Map<AccessPoint, Double>>();

			for (int l = 0; l < L; l++) {
				if (0d == vars[l]) // no cloudlet placed in this location l.
					continue;

				Node location = this.potentialLocations.get(l);
				if (null == queryAssignment.get(location))
					queryAssignment.put(location, new HashMap<AccessPoint, Double>());

				for (int j = 0; j < n; j++) {
					if (0d == vars[L + j * L + l])
						continue;

					AccessPoint currAP = this.APs.get(j);
					double assigned = vars[L + j * L + l];
					queryAssignment.get(location).put(currAP, assigned);
				}
			}

			// delete the problem and free memory
			solver.deleteLp();

		} catch (LpSolveException e) {
			e.printStackTrace();
		}

		return queryAssignment;
	}	
	
	
	private void setEdgeWeightDataCenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		for (InternetLink il : datacenterNetwork.edgeSet()) {
			datacenterNetwork.setEdgeWeight(il, il.getCost());
		}
	}

	public SamplePlacementSimulator getSimulator() {
		return simulator;
	}

	public void setSimulator(SamplePlacementSimulator simulator) {
		this.simulator = simulator;
	}
	

	public List<Integer> getAdmittedNumOfQueriesPerTS() {
		return admittedNumOfQueriesPerTS;
	}

	public void setAdmittedNumOfQueriesPerTS(List<Integer> admittedNumOfQueriesPerTS) {
		this.admittedNumOfQueriesPerTS = admittedNumOfQueriesPerTS;
	}

	public List<Double> getCostPerTS() {
		return costPerTS;
	}

	public void setCostPerTS(List<Double> costPerTS) {
		this.costPerTS = costPerTS;
	}

	public List<Double> getStorageCostPerTS() {
		return storageCostPerTS;
	}

	public void setStorageCostPerTS(List<Double> storageCostPerTS) {
		this.storageCostPerTS = storageCostPerTS;
	}

	public List<Double> getUpdateCostPerTS() {
		return updateCostPerTS;
	}

	public void setUpdateCostPerTS(List<Double> updateCostPerTS) {
		this.updateCostPerTS = updateCostPerTS;
	}

	public List<Double> getAccessCostPerTS() {
		return accessCostPerTS;
	}

	public void setAccessCostPerTS(List<Double> accessCostPerTS) {
		this.accessCostPerTS = accessCostPerTS;
	}

	public List<Double> getProcessCostPerTS() {
		return processCostPerTS;
	}

	public void setProcessCostPerTS(List<Double> processCostPerTS) {
		this.processCostPerTS = processCostPerTS;
	}
}
