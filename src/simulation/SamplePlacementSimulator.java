package simulation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.SimpleWeightedGraph;

import algs.Benchmark1;
import algs.ProposedApproximationAlg;
import algs.ProposedHeuristicAlg;
import graph.InternetLink;
import graph.NetworkGenerator;
import graph.Node;
import system.DataCenter;
import system.Dataset;
import system.Query;
import system.Sample;
import utils.IdAllocator;
import utils.RanNum;

public class SamplePlacementSimulator {
	private SimpleWeightedGraph<Node, InternetLink> datacenterNetwork = null;
	public static IdAllocator idAllocator = new IdAllocator();
	
	private List<DataCenter> dataCenterList = new ArrayList<DataCenter>();
	// <trial, list of datasets>
	private Map<Integer, List<Dataset>> datasets = new HashMap<Integer, List<Dataset>>();
	private Map<Integer, List<Sample>> samples = new HashMap<Integer, List<Sample>>();
	private Map<Integer, List<Query>> queries = new HashMap<Integer,List<Query>>();

	public SamplePlacementSimulator() {
		
	}

	public static void main(String[] s) {
		//performanceHeuristic();
		//performanceApproximation();
		//performance();
		performanceHeuAppro();
		//impactOfErrorOnCost();
	}
	
	public static void performance() {
		
		int numOfAlgs = 2; 
		//int [] network_sizes = {20, 30, 40, 50, 100, 150, 200}; 
		int [] network_sizes = {50};
		double [][] aveCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveStorageCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveUpdateCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveAccessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveProcessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveError = new double [numOfAlgs][network_sizes.length];

		
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			
			for(int round = 0; round < Parameters.roundNum; round ++) {// different network toplolgies.
				String networkIndexPostFix = "";
				if (round > 0) 
					networkIndexPostFix = "-" + round;
				
				SamplePlacementSimulator simulator = new SamplePlacementSimulator();
				simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix, network_size);//get the data center network (cloud network)			
				simulator.InitializeDatasetsAndSamples();
				simulator.InitializeQueries();
				ProposedHeuristicAlg heuAlg = new ProposedHeuristicAlg(simulator);
				heuAlg.run();
				
				double averageCostT = 0d;
				double averageStorageCostT = 0d;
				double averageUpdateCostT = 0d;
				double averageAccessCostT = 0d;
				double averageProcessCostT = 0d;
				double averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (heuAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (heuAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (heuAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (heuAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (heuAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
					averageErrorT += (heuAlg.getAverageErrorTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[0][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[0][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[0][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[0][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[0][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[0][i] += (averageErrorT / Parameters.roundNum);
				
				for (Node node : simulator.getDatacenterNetwork().vertexSet()) {
					if (node instanceof DataCenter) {
						((DataCenter) node).reset();
					}
				}
				for (InternetLink il : simulator.getDatacenterNetwork().edgeSet()) {
					il.clear();
				}
				
				for (int t = 0; t < Parameters.numOfTrials; t ++){
					for (Dataset ds : simulator.getDatasets().get(t)){
						ds.reset();
					}
				}
				
				Benchmark1 benchmarkAlg = new Benchmark1(simulator);
				benchmarkAlg.run();
				
				averageCostT = 0d;
				averageStorageCostT = 0d;
				averageUpdateCostT = 0d;
				averageAccessCostT = 0d;
				averageProcessCostT = 0d;
				averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (benchmarkAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (benchmarkAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (benchmarkAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (benchmarkAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (benchmarkAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
					averageErrorT += (benchmarkAlg.getAverageErrorTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[1][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[1][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[1][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[1][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[1][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[1][i] += (averageErrorT / Parameters.roundNum);
				//aveCost[2][i] += (averageCostLowerboundT / Parameters.roundNum);
			}
		}
		
		System.out.println("total costs and lower bound");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("average errors");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveError[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("storage costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveStorageCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("update costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveUpdateCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("access costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveAccessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("process costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs ; j ++){
				out += aveProcessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
	}
	
	public static void performanceHeuAppro() {
		
		int numOfAlgs = 2; 
		//int [] network_sizes = {20, 30, 40, 50, 100, 150, 200}; 
		int [] network_sizes = {50};
		double [][] aveCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveStorageCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveUpdateCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveAccessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveProcessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveError = new double [numOfAlgs][network_sizes.length];

		
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			
			for(int round = 0; round < Parameters.roundNum; round ++) {// different network toplolgies.
				String networkIndexPostFix = "";
				if (round > 0) 
					networkIndexPostFix = "-" + round;
				
				SamplePlacementSimulator simulator = new SamplePlacementSimulator();
				simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix, network_size);//get the data center network (cloud network)			
				simulator.InitializeDatasetsAndSamples();
				simulator.InitializeQueries();
				ProposedHeuristicAlg heuAlg = new ProposedHeuristicAlg(simulator);
				heuAlg.run();
				
				double averageCostT = 0d;
				double averageStorageCostT = 0d;
				double averageUpdateCostT = 0d;
				double averageAccessCostT = 0d;
				double averageProcessCostT = 0d;
				double averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (heuAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (heuAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (heuAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (heuAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (heuAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
					averageErrorT += (heuAlg.getAverageErrorTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[0][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[0][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[0][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[0][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[0][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[0][i] += (averageErrorT / Parameters.roundNum);
				
				for (Node node : simulator.getDatacenterNetwork().vertexSet()) {
					if (node instanceof DataCenter) {
						((DataCenter) node).reset();
					}
				}
				for (InternetLink il : simulator.getDatacenterNetwork().edgeSet()) {
					il.clear();
				}
				
				for (int t = 0; t < Parameters.numOfTrials; t ++){
					for (Dataset ds : simulator.getDatasets().get(t)){
						ds.reset();
					}
				}
				
				ProposedApproximationAlg approAlg = new ProposedApproximationAlg(simulator);
				approAlg.run();
				
				averageCostT = 0d;
				averageStorageCostT = 0d;
				averageUpdateCostT = 0d;
				averageAccessCostT = 0d;
				averageProcessCostT = 0d;
				int numOfInvalidTrials = 0;	
				averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					if (approAlg.getCostTrials().get(t) == 0d) numOfInvalidTrials ++; 
					averageCostT += (approAlg.getCostTrials().get(t));
					averageStorageCostT += (approAlg.getStorageCostTrials().get(t));
					averageUpdateCostT += (approAlg.getUpdateCostTrials().get(t));
					averageAccessCostT += (approAlg.getAccessCostTrials().get(t));
					averageProcessCostT += (approAlg.getProcessCostTrials().get(t));
					averageErrorT += (approAlg.getAverageErrorTrials().get(t));
				}
				
				averageCostT /= (Parameters.numOfTrials - numOfInvalidTrials);
				averageStorageCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageUpdateCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageAccessCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageProcessCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageErrorT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				
				aveCost[1][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[1][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[1][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[1][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[1][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[1][i] += (averageErrorT / Parameters.roundNum);
				//aveCost[2][i] += (averageCostLowerboundT / Parameters.roundNum);
			}
		}
		
		System.out.println("total costs and lower bound");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("average errors");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveError[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("storage costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveStorageCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("update costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveUpdateCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("access costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveAccessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("process costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs ; j ++){
				out += aveProcessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
	}

	
	public static void impactOfErrorOnCost() {
		int numOfAlgs = 2; 
		//int [] network_sizes = {20, 30, 40, 50, 100, 150, 200}; 
		double [] lowestErrors = {0.05, 0.075, 0.1, 0.125};
		
		double [][] aveCost = new double [numOfAlgs][lowestErrors.length];
		double [][] aveStorageCost = new double [numOfAlgs][lowestErrors.length];
		double [][] aveUpdateCost = new double [numOfAlgs][lowestErrors.length];
		double [][] aveAccessCost = new double [numOfAlgs][lowestErrors.length];
		double [][] aveProcessCost = new double [numOfAlgs][lowestErrors.length];
		double [][] aveError = new double [numOfAlgs][lowestErrors.length];

		
		for (int i = 0; i < lowestErrors.length; i ++) {
			int network_size = 50;
			
			double lowestError = lowestErrors[i];
			
			for(int round = 0; round < Parameters.roundNum; round ++) {// different network toplolgies.
				String networkIndexPostFix = "";
				if (round > 0) 
					networkIndexPostFix = "-" + round;
				
				SamplePlacementSimulator simulator = new SamplePlacementSimulator();
				simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix, network_size);//get the data center network (cloud network)			
				simulator.InitializeDatasetsAndSamples();
				simulator.InitializeQueries();
				ProposedHeuristicAlg heuAlg = new ProposedHeuristicAlg(simulator);
				heuAlg.setLowestError(lowestError);
				heuAlg.run();
				
				double averageCostT = 0d;
				double averageStorageCostT = 0d;
				double averageUpdateCostT = 0d;
				double averageAccessCostT = 0d;
				double averageProcessCostT = 0d;
				double averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (heuAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (heuAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (heuAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (heuAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (heuAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
					averageErrorT += (heuAlg.getAverageErrorTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[0][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[0][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[0][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[0][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[0][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[0][i] += (averageErrorT / Parameters.roundNum);
				
				for (Node node : simulator.getDatacenterNetwork().vertexSet()) {
					if (node instanceof DataCenter) {
						((DataCenter) node).reset();
					}
				}
				for (InternetLink il : simulator.getDatacenterNetwork().edgeSet()) {
					il.clear();
				}
				
				for (int t = 0; t < Parameters.numOfTrials; t ++){
					for (Dataset ds : simulator.getDatasets().get(t)){
						ds.reset();
					}
				}
				
				Benchmark1 benchmarkAlg = new Benchmark1(simulator);
				benchmarkAlg.setLowestError(lowestError);
				benchmarkAlg.run();
				
				averageCostT = 0d;
				averageStorageCostT = 0d;
				averageUpdateCostT = 0d;
				averageAccessCostT = 0d;
				averageProcessCostT = 0d;
				averageErrorT = 0d; 
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (benchmarkAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (benchmarkAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (benchmarkAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (benchmarkAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (benchmarkAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
					averageErrorT += (benchmarkAlg.getAverageErrorTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[1][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[1][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[1][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[1][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[1][i] += (averageProcessCostT / Parameters.roundNum);
				aveError[1][i] += (averageErrorT / Parameters.roundNum);
				//aveCost[2][i] += (averageCostLowerboundT / Parameters.roundNum);
			}
		}
		
		System.out.println("total costs and lower bound");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double lowestError = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveCost[j][i] + " ";
			}
			System.out.println("" + lowestError + " " + out);
		}
		
		System.out.println("average errors");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double lowestError = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveError[j][i] + " ";
			}
			System.out.println("" + lowestError + " " + out);
		}
		
		System.out.println("storage costs");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double lowestError = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveStorageCost[j][i] + " ";
			}
			System.out.println("" + lowestError + " " + out);
		}
		
		System.out.println("update costs");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double lowestError = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveUpdateCost[j][i] + " ";
			}
			System.out.println("" + lowestError + " " + out);
		}
		
		System.out.println("access costs");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double lowestError = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveAccessCost[j][i] + " ";
			}
			System.out.println("" + lowestError + " " + out);
		}
		
		System.out.println("process costs");
		for (int i = 0; i < lowestErrors.length; i ++) {
			double network_size = lowestErrors[i];
			String out = "";
			for (int j = 0; j < numOfAlgs ; j ++){
				out += aveProcessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
	}
	
	public static void performanceHeuristic() {
		
		int numOfAlgs = 1; 
		//int [] network_sizes = {20, 30, 40, 50, 100, 150, 200}; 
		int [] network_sizes = {20, 30};
		double [][] aveCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveStorageCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveUpdateCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveAccessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveProcessCost = new double [numOfAlgs][network_sizes.length];
		
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			
			for(int round = 0; round < Parameters.roundNum; round ++) {// different network toplolgies.
				String networkIndexPostFix = "";
				if (round > 0) 
					networkIndexPostFix = "-" + round;
				
				SamplePlacementSimulator simulator = new SamplePlacementSimulator();
				simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix, network_size);//get the data center network (cloud network)			
				simulator.InitializeDatasetsAndSamples();
				simulator.InitializeQueries();
				ProposedHeuristicAlg heuAlg = new ProposedHeuristicAlg(simulator);
				heuAlg.run();
				
				double averageCostT = 0d;
				double averageStorageCostT = 0d; 
				double averageUpdateCostT = 0d; 
				double averageAccessCostT = 0d;
				double averageProcessCostT = 0d;
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					averageCostT += (heuAlg.getCostTrials().get(t) / Parameters.numOfTrials);
					averageStorageCostT += (heuAlg.getStorageCostTrials().get(t) / Parameters.numOfTrials);
					averageUpdateCostT += (heuAlg.getUpdateCostTrials().get(t) / Parameters.numOfTrials);
					averageAccessCostT += (heuAlg.getAccessCostTrials().get(t) / Parameters.numOfTrials);
					averageProcessCostT += (heuAlg.getProcessCostTrials().get(t) / Parameters.numOfTrials);
				}
				
				aveCost[0][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[0][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[0][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[0][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[0][i] += (averageProcessCostT / Parameters.roundNum);
			}
		}
		
		System.out.println("total costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("storage costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveStorageCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("update costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveUpdateCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("access costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveAccessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("process costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveProcessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
	}
	
	public static void performanceApproximation() {
		int numOfAlgs = 1; 
		int [] network_sizes = {20, 30}; // 100, 150, 200
		double [][] aveCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveStorageCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveUpdateCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveAccessCost = new double [numOfAlgs][network_sizes.length];
		double [][] aveProcessCost = new double [numOfAlgs][network_sizes.length];
		
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			
			for(int round = 0; round < Parameters.roundNum; round ++) {// different network toplolgies.
				String networkIndexPostFix = "";
				if (round > 0) 
					networkIndexPostFix = "-" + round;
			
				SamplePlacementSimulator simulator = new SamplePlacementSimulator();
				simulator.InitializeDataCenterNetwork(simulator.getDatacenterNetwork(), networkIndexPostFix, network_size);//get the data center network (cloud network)			
				simulator.InitializeDatasetsAndSamples();
				simulator.InitializeQueries();
				ProposedApproximationAlg approAlg = new ProposedApproximationAlg(simulator);
				approAlg.run();
				
				double averageCostT = 0d;
				double averageStorageCostT = 0d;
				double averageUpdateCostT = 0d;
				double averageAccessCostT = 0d;
				double averageProcessCostT = 0d;
				int numOfInvalidTrials = 0;				
				for (int t = 0; t < Parameters.numOfTrials; t++) {
					if (approAlg.getCostTrials().get(t) == 0d) numOfInvalidTrials ++; 
					averageCostT += (approAlg.getCostTrials().get(t));
					averageStorageCostT += (approAlg.getStorageCostTrials().get(t));
					averageUpdateCostT += (approAlg.getUpdateCostTrials().get(t));
					averageAccessCostT += (approAlg.getAccessCostTrials().get(t));
					averageProcessCostT += (approAlg.getProcessCostTrials().get(t));
				}
				
				averageCostT /= (Parameters.numOfTrials - numOfInvalidTrials);
				averageStorageCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageUpdateCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageAccessCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				averageProcessCostT /= (Parameters.numOfTrials - numOfInvalidTrials); 
				
				aveCost[0][i] += (averageCostT / Parameters.roundNum);
				aveStorageCost[0][i] += (averageStorageCostT / Parameters.roundNum);
				aveUpdateCost[0][i] += (averageUpdateCostT / Parameters.roundNum);
				aveAccessCost[0][i] += (averageAccessCostT / Parameters.roundNum);
				aveProcessCost[0][i] += (averageProcessCostT / Parameters.roundNum);
			}
		}
		
		System.out.println("total costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("storage costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveStorageCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("update costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveUpdateCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("access costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveAccessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
		
		System.out.println("process costs");
		for (int i = 0; i < network_sizes.length; i ++) {
			int network_size = network_sizes[i];
			String out = "";
			for (int j = 0; j < numOfAlgs; j ++){
				out += aveProcessCost[j][i] + " ";
			}
			System.out.println("" + network_size + " " + out);
		}
	}

	/************************************
	 * Initialization functions
	 * 
	 * @return
	 ************************************/
	
	public SimpleWeightedGraph<Node, InternetLink> InitializeDataCenterNetwork(
			SimpleWeightedGraph<Node, InternetLink> dcNet, String networkIndexPostfix, int numOfDCs) {

		if (null == dcNet) {

			ConnectivityInspector<Node, InternetLink> connect = null;
			do {
				// Initialize the data center network
				datacenterNetwork = new SimpleWeightedGraph<Node, InternetLink>(InternetLink.class);
				NetworkGenerator<Node, InternetLink> networkGenerator = new NetworkGenerator<Node, InternetLink>();

				networkGenerator.setSize(numOfDCs);
				networkGenerator.setGenerateType(0);// generate data center
													// networks
				networkGenerator.setNetworkIndexPostFix(networkIndexPostfix);
				networkGenerator.generateGraph(datacenterNetwork, null, null);
				// displayGraph(substrateNetwork);
				connect = new ConnectivityInspector<Node, InternetLink>(datacenterNetwork);
			} while (!connect.isGraphConnected());

			List<DataCenter> dcs = this.getDataCenters();

			for (DataCenter dc1 : dcs) {
				for (DataCenter dc2 : dcs) {
					InternetLink il = datacenterNetwork.getEdge(dc1, dc2);
					if (null != il) {
						il.setEdgeSource(dc1);
						il.setEdgeTarget(dc2);
					}
				}
			}
		} else {// clear some parameters of dcNet.
			// Initialize parameters of sNet.
			// modify.
			for (Node node : dcNet.vertexSet()) {
				if (node instanceof DataCenter) {
					DataCenter dc = (DataCenter) node;
					dc.clearAdmittedSamples();
				} // TODO double check whether no need for front end servers.
			}
			for (InternetLink il : dcNet.edgeSet()) {
				il.clearUsers();
			}

			this.datacenterNetwork = dcNet;
		}
		return datacenterNetwork;
	}
	
	public Map<Integer, List<Dataset>> InitializeDatasetsAndSamples() {
		
		if (this.datasets.isEmpty()) {
			for(int i = 0; i < Parameters.numOfTrials; i ++) {
				int numOfDatasetsPerTS = RanNum.getRandomIntRange(Parameters.maxNumOfDatasetsPerTS, Parameters.minNumOfDatasetsPerTS);
				List<Dataset> dss = new ArrayList<Dataset>();
				List<Sample> sams = new ArrayList<Sample>();
				for(int j = 0; j < numOfDatasetsPerTS; j ++) {
					Dataset ds = new Dataset(this.getDataCenters());
					dss.add(ds);
					for (int s = 0; s < Parameters.errorBounds.length; s ++){
						Sample sample = new Sample(ds, s);
						sams.add(sample);
						ds.getSamples().add(sample);
					}
				}
				this.datasets.put(i, dss);
				this.samples.put(i, sams);
			}
		}
		
		return this.datasets;
	}
	
	public Map<Integer, List<Query>> InitializeQueries() {
		
		//String numOfQueries = "";
		if (this.queries.isEmpty()){// generate queries
			for(int i = 0; i < Parameters.numOfTrials; i ++) {
				
				int numOfQueriesPerTS = RanNum.getRandomIntRange(Parameters.maxNumOfQueriesPerTS, Parameters.minNumOfQueriesPerTS);
				//numOfQueries += numOfQueriesPerTS + " ";
				List<Query> qus = new ArrayList<Query>();
				Set<Sample> samplesToBePlaced = new HashSet<Sample>();
				double targetError = Parameters.errorBounds[1];
				
				for(int j = 0; j < numOfQueriesPerTS; j ++) {
					Query quer = new Query(this.getDataCenters(), this.getDatasets().get(i));
					qus.add(quer);
					
					for (Dataset ds : quer.getDatasets()) {
						samplesToBePlaced.add(ds.getSample(targetError));
					}
				}
				queries.put(i, qus);
				
				// make sure the total computing resource is more than that required by all queries. 
				double totalComputingDemand = 0d;
				for (Sample sample : samplesToBePlaced) {
					totalComputingDemand += sample.getVolume() * Parameters.computingAllocatedToUnitData; 
				}
				
				double totalComputingAvailable = 0d; 
				for (DataCenter dc : this.dataCenterList){
					totalComputingAvailable += dc.getAvailableComputing();
				}
				
				if (totalComputingDemand > totalComputingAvailable){
					double ratio = totalComputingAvailable / (totalComputingDemand + 10d);
					
					System.out.println("scalling down the size of all datasets by ratio: " + ratio);

					for (Dataset ds : this.datasets.get(i)) {
						ds.setVolume(ds.getVolume() * ratio);
						for (Sample sam : ds.getSamples()){
							sam.setVolume(sam.getVolume() * ratio);
						}
					}
				}
			}
			//System.out.println(numOfQueries);
		}
		return this.queries;
	}

	public void modifyCosts() {

		for (Node node : this.datacenterNetwork.vertexSet()) {

			if (node instanceof DataCenter) {
				double processingCost = RanNum.getRandomDoubleRange(Parameters.processCostUnitDataMax,
						Parameters.processCostUnitDataMin);
				double storageCost = RanNum.getRandomDoubleRange(Parameters.storageCostUnitDataMax,
						Parameters.storageCostUnitDataMin);

				DataCenter dc = (DataCenter) node;
				dc.setProcessingCost(processingCost);
				dc.setStorageCost(storageCost);
			}
		}

		for (InternetLink il : this.datacenterNetwork.edgeSet()) {
			il.setCost(RanNum.getRandomDoubleRange(Parameters.maxBandwidthCost, Parameters.minBandwidthCost));
		}
	}
	
	public List<Sample> samplesOfDataset(int trial, Dataset ds){
		
		List<Sample> foundSamples = new ArrayList<Sample>();
		for (Sample sample : this.getSamples().get(trial)){
			if (sample.getParentDataset().equals(ds))
				foundSamples.add(sample);
		}
		
		return foundSamples;
	}

	public List<DataCenter> getDataCenters() {
		if (!this.dataCenterList.isEmpty())
			return this.dataCenterList;

		for (Node dc : this.datacenterNetwork.vertexSet()) {
			if (dc instanceof DataCenter)
				this.dataCenterList.add((DataCenter) dc);
		}

		return this.dataCenterList;
	}

	public int getNumberOfDataCenters() {
		int numDC = 0;
		for (Node dc : this.datacenterNetwork.vertexSet()) {
			if (dc instanceof DataCenter)
				numDC++;
		}
		return numDC;
	}

	/******** getter and setter ************/
	public SimpleWeightedGraph<Node, InternetLink> getDatacenterNetwork() {
		return datacenterNetwork;
	}

	public void setDatacenterNetwork(SimpleWeightedGraph<Node, InternetLink> datacenterNetwork) {
		this.datacenterNetwork = datacenterNetwork;
	}

	public Map<Integer, List<Dataset>> getDatasets() {
		return datasets;
	}

	public void setDatasets(Map<Integer, List<Dataset>> datasets) {
		this.datasets = datasets;
	}

	public Map<Integer, List<Sample>> getSamples() {
		return samples;
	}

	public void setSamples(Map<Integer, List<Sample>> samples) {
		this.samples = samples;
	}

	public Map<Integer, List<Query>> getQueries() {
		return queries;
	}

	public void setQueries(Map<Integer, List<Query>> queries) {
		this.queries = queries;
	}

}
