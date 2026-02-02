/*
 * COMP2230 - Algorithms
 * Assignment 1 - Answers
 * @author  Amanda Foxley - c3137540
 * @version 1.1
 * Edited 26/8/24
 */

// -----------------------------------------------
// The TrafficAnalyser class loads a string of road names, intersection names and average travel times into a HashMap which represents a graph of the city.
//
// It uses a number of global variables containing information relevant to the city, including inner city edges and road speeds, as well as a number of look up HashMaps to convert
// between the index where the node is stored (in the relevant HashMap) and the road/intersection names.
// It also contains a number of methods to query and search this the graph information. e.g. isInInnerCity() cityBottleneckRoads() etc
// -----------------------------------------------

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

public class TrafficAnalyser {
	private final MapGenerator mapGen;
	private HashMap<String, GraphNode> graph; // Stores all nodes within the cityMap (Intersection name and graphNode object)
	private HashMap<String, GraphNode> innerCity; //Inner city, subgraph with the highest number of nodes (Intersection name and graphNode object)
	private HashMap<GraphEdge, Double> innerCityEdgeMap; //Inner city road names and their average travel times (Road name and avg travel time)
	private ArrayList<Double> innerCityRoadTimes;	//Sorted array of all inner city road/edge travel times

	//Lookup tables - for converting between string (names) to integer (indexes)
	private HashMap<Integer, String> indexToName;  //Index value to intersection name
	private HashMap<String, Integer> nameToIndex;  //Intersection name to index value

	//For disjoint sets
	private int[] parent;			
	private int[] rank;

	public String cityMap = null;  // We make this public in order to access it for testing, you would not normally do this 

	public TrafficAnalyser(int seed){
		mapGen = new MapGenerator(seed); // Pass a seed to the random number generator, allows for reproducibility
	}

	/*
	 * Precondition: None.
	 * Postcondition: The city map has been loaded into a graph HashMap.
	 * 		          Global lookup HashMaps have also been populated
	 */
	public void loadMap(){
		if (cityMap == null) {
			cityMap = mapGen.generateMap();
		}

		//Initialise global variables
		graph = new HashMap<String, GraphNode>();
		innerCity = new HashMap<String, GraphNode>();
		innerCityEdgeMap = new HashMap<GraphEdge, Double>();
		innerCityRoadTimes = new ArrayList<Double>();
		
		String map = cityMap; //Create another copy of the map so we can remove some of the symbols from the string

		//Remove bracket symbols
		map = map.replace("{", "");
		map = map.replace("}", "");
		String[] mapCopy = map.split(","); //spilt the map into parts based on "," symbol

		for(int i = 0; i < mapCopy.length; i=i+4) {//Add 4 each iteration since each edge has 4 parts

			//Remove any leading spaces on front of road an intersection names
			if(mapCopy[i].startsWith(" ")) { //Not all road names contain a space at the front
				mapCopy[i] = mapCopy[i].substring(1, mapCopy[i].length()); //Remove " " at the start of any road names
			}        	
			mapCopy[i+1] = mapCopy[i+1].substring(1, mapCopy[i+1].length()); //Remove " " at the start of any intersection 1 names
			mapCopy[i+2] = mapCopy[i+2].substring(1, mapCopy[i+2].length()); //Remove " " at the start of any intersection 2 names

			//Node/intersection names
			String source = mapCopy[i+1];	
			String destination = mapCopy[i+2];	

			GraphEdge newRoad = new GraphEdge(mapCopy[i], source, destination, Double.parseDouble(mapCopy[i+3])); 

			//Call method to add the nodes and edges to the graph
			addEdge(source, destination, newRoad);
		}  

		populateLookupHashMaps();	
		populateInnerCityGraph(); 	//Populate inner city using disjoint sets
		populateInnerCityEdgeMap(); //BFS for all inner city edges
	}	
	


	/**
	 * Precondition:  None
	 * Postcondition: The graph has been populated with the relevant nodes
	 * 				  road information has been added to the sourceNode and destinationNode (the road that runs between them) this is added to each node's edges List
	 * 		          neighbouringNode information has been added to the sourceNode and destinationNode (nodes that are 1 hop away) this is added to each node's neighbouringNode HashMap
	 * 
	 * @param source - name of the origin intersection, destination - name of the destination intersection, GraphEdge newRoad - the road you can traverse between the source and destination intersections
	 * @return None
	 */
	private void addEdge(String source, String destination, GraphEdge newRoad) {

		//Retrieve the node's information from the graph HashMap, if they don't exist in the HashMap then create new GraphNode objects instead
		GraphNode sourceNode = graph.getOrDefault(source, new GraphNode(source));
		GraphNode destinationNode = graph.getOrDefault(destination, new GraphNode(destination));

		// Update the nodes in the graph HashMap
		graph.put(source, sourceNode);
		graph.put(destination, destinationNode);

		// Add the edge/road to the each nodes edge list
		sourceNode.edges.add(newRoad);
		destinationNode.edges.add(newRoad);
		
		//Add nodes that are one hop away to each node's neighbouringNodes HashMap
		sourceNode.neighbourNodes.put(newRoad.name, destinationNode);
		destinationNode.neighbourNodes.put(newRoad.name, sourceNode);
	}
	
	/**
	 * Precondition:  graph HashMap is not empty
	 * Postcondition: indexToName, nameToIndex, have been populated with their relevant lookup values
	 * 
	 * @param None
	 * @return None
	 */
	private void populateLookupHashMaps() {
		indexToName = new HashMap<Integer, String>();
		nameToIndex = new HashMap<String, Integer>();
		
		//Populate the lookup HashMaps       
		int j =0;
		for (String key : graph.keySet()) {
			indexToName.put(j, key);     
			nameToIndex.put(key, j++);
		}
	}
			
	/**
	 * Precondition:  graph HashMap is not empty
	 * Postcondition: innerCity HashMap has been populated with relevant node values.
	 * 			      This method uses disjoint sets with ranks to identify the connected components of the map of the city.
	 * 			      The largest connected component is the inner city of the map.
	 * 				  Call the make() method to make each node in the graph it's own parent. Then call the unionNodes() method to join all nodes of the graph with other nodes that are linked to them via roads.
	 * 
	 * @param None
	 * @return None
	 */
	private void populateInnerCityGraph() {
		make(); //Call method to make each node it's own parent, in global parent array
		
		//Union all nodes of the graph
		Set<GraphNode> visited = new HashSet<>();
		for (GraphNode node : graph.values()) {
			if (!visited.contains(node)) { 		//If the node isn't in visited then it is part of a subgraph that hasn't been traversed. Send to unionNodes() method to traverse and find all linked nodes of this subgraph.
				unionNodesBFS(node.name, visited);
			}
		}

		//Find the highest ranked node in the graph - This is the parent of the innerCity
		int highestRankIndex = 0; 
		for(int i = 0; i < rank.length; i++) {

			if(rank[i] > rank[highestRankIndex]) { 
				highestRankIndex = i;
			}
		}
		
		//Using the index of parent of innerCity (highestRankIndex), add all nodes where their parent's index is the same i.e. if node 2's parent = highestRankIndex, then node 2 is part of inner city
		for(int i = 0; i < parent.length; i++) { //Populate the innerCity HashMap

			if(parent[i] == highestRankIndex) {// if the parent = highestRankedNode then they are in the same set
				String nodeName = indexToName.get(i);
				innerCity.put(nodeName, graph.get(nodeName));
			}
		}
	}
	
	/*
	 * Precondition: graph HashMap is not null
	 * Postcondition: All nodes have been set to be their own parent in the parent's array. All rank's have been set to 0.
     * 
	 * @param None
	 * @return None
	 */
	private void make() {
		parent = new int[graph.size()];			
		rank = new int[graph.size()];	
		
		//Make all nodes singletons and rank of all nodes is 0
		for(int i = 0; i < graph.size(); i++) {
			parent[i] = i;
			rank[i] = 0;
		}
	}
	
	/*
	 * Precondition: graph HashMap is not null
	 * Postcondition: The graph has been traversed and all nodes linked to each other through a mutual edge have been grouped using union() method
     * 			      A group of nodes is a subgraph of the map. A group of nodes all have the same parent in the global parent array.
     * 			      The values stored in the global array rank indicates how many nodes are in that group.
     * 
	 * @param startNode - name of the first node to traverse using iterative BFS algorithm. visited - Set to store nodes visited in BFS search
	 * @return None
	 */
	private void unionNodesBFS(String startNode, Set<GraphNode> visited) {
		Queue<GraphNode> openQueue = new LinkedList<>();
		GraphNode node = graph.get(startNode);

		openQueue.offer(node);
		
		//If we want the rank array to indicate how many nodes are in a subgraph, we need to add 1 to the parent's rank, i.e. the startNode (Since it counts as a node in the subgraph as well)
		int startNodeIndex = nameToIndex.get(startNode);
		rank[startNodeIndex] = rank[startNodeIndex] + 1; 
		
		//Iterate through graph, when a neighbour is reached through a mutual edge, the union method is called to form a relationship using disjoint sets (parent and rank arrays)
		while (!openQueue.isEmpty()) {
			GraphNode current = openQueue.poll();
			visited.add(current);
			
			for (GraphEdge edge : current.edges) {
				
				GraphNode neighbour = edge.intersection1.equals(current.name) ? graph.get(edge.intersection2) : graph.get(edge.intersection1); 
				
				if (!visited.contains(neighbour) && !openQueue.contains(neighbour)) {
					openQueue.offer(neighbour);
					
					//Get node indexes
					int currentNode = nameToIndex.get(current.name);
					int neighbourNode = nameToIndex.get(neighbour.name);
					
					union(neighbourNode, currentNode); //Current node will be the parent of neighbour in parent array
				}
			}
		}
	}
	
	/**
	 * Precondition:  i and j are valid integers
	 * Postcondition: The height of the tree has been flattened by calling find (find each node's parent)
	 * 				  Each node has been changed to the value that was stored in it's parent
	 * 				  The node with the higher rank will be set to be the parent of the node with the lower rank
	 * 				  The rank of the parent will be increased since it has had a node added to it's set
	 * 
	 * @param int i and int j - two integers representing nodes in separate sets
	 * @return None
	 */
	private void union(int i,int j) {
		i = find(i);
		j = find(j);

		if(rank[i] < rank[j]) { //i is smaller attach to j (or j's parent)
			parent[i] = j;
			rank[j] = rank[j] + 1; 	 	//Increase rank of j since we've added a member to j's set 
		}else if(rank[i] > rank[j]) {
			parent[j] = i;
			rank[i] = rank[i] + 1;		//Increase rank of i since we've added a member to i's set 
		}
		else {
			parent[i] = j;
			rank[j] = rank[j] + 1;		//Increase rank of j since we've added a member to j's set 
		}
	}
	
	/**
	 * Precondition:  parent array is not null
	 * Postcondition: The root/parent of the set that contains k has been returned
	 * 				  This helps with flattening the structure of the disjoint set. Instead of attaching nodeA to another nodeB, you attach a nodeA to the parent of nodeB.
	 * 				  This method finds that parent.
	 * 
	 * @param int k - A member of a set
	 * @return int r - returns the parent of k
	 */
	private int find(int k) {
		int r = k;

		while(r != parent[r]) {
			r = parent[r];
		}

		while(parent[k] != r) {
			int j = parent[k];
			parent[k] = r;
			k = j;
		}
		return r;  
	}
		
	/*
	 * Precondition: innerCity HashMap is not empty
	 * Postcondition: The inner city subgraph has been traversed (using BFS) and all roads in the component have been added to innerCityEdgeMap (road name and avg speed) HashMap and 
	 *                innerCityRoadSpeeds (avg road speeds) ArrayList. innerCityRoadSpeeds Arraylist has been sorted.
	 * @param None
	 * @return None
	 */
	private void populateInnerCityEdgeMap() {

		Set<GraphNode> visited = new HashSet<>();
		GraphNode startNode = (GraphNode) innerCity.values().toArray()[0]; // Get the first node in the inner city list
		Queue<GraphNode> openQueue = new LinkedList<>();

		openQueue.add(startNode); //Insert the first inner city node into the queue

		while (!openQueue.isEmpty()) {
			GraphNode current = openQueue.poll(); //Get the next node to be traversed
			visited.add(current);

			for (GraphEdge edge : current.edges) { //Check all edges of current node
				GraphNode neighbour = edge.intersection1.equals(current.name) ? innerCity.get(edge.intersection2) : innerCity.get(edge.intersection1); //Using the road information, get the next node (i.e. the one that isn't the current node)

				//Add the road to the HashMap if we havn't already
				if (!innerCityEdgeMap.containsKey(edge)) { 
					innerCityEdgeMap.put(edge, edge.avgTravel);
					innerCityRoadTimes.add(edge.avgTravel); //Add the road speed to the roadSpeed arrayList
				}

				//Add neighbour to visited HashMap so we don't visit this again
				if (!visited.contains(neighbour)) {
					openQueue.add(neighbour);
				}
			}
		}

		Collections.sort(innerCityRoadTimes); // Sort the innerCityRoadSpeeds ArrayList (Sorted by travel time)
	}
	
	/**
	 * Precondition:  innerCity HashMap is not empty
	 * Postcondition: populateInnerCityGraph has been called to populate the innerCity subgraph (global HashMap) using disjoint sets data structure
	 * 
	 * 				  Note: The 'inner city' is defined as the largest connected component of intersections in the map.
	 * @param String intersectionName - The name of the intersection to check
	 * @return boolean - true if the intersection is in the inner city, false otherwise
	 */
	public boolean isInInnerCity(String intersectionName) {
		return innerCity.containsKey(intersectionName);
	}

	/*
	 * Precondition:  innerCityRoadSpeeds is not null
	 * Postcondition: The number of roads in the inner city with a speed greater than the threshold has been returned.
	 * 				  This method has called the binarySearch method to find the index of the next highest number from the threshold.
	 * 				  The index returned has been used to discard all values from the roadCount Array before and including the threshold. 
	 * 				  The count of the roads remaining in the array has been returned.
	 * 
	 * @param double threshold - The travel time threshold in minutes
	 * @return int - The number of roads in the 'inner city' that have a travel time strictly greater (>) than the threshold
	 */
	public int countInnerCitySlowRoads(double threshold) {

		List<Double> roadCount = new ArrayList<Double>();

		int foundItemIndex = binarySearch(innerCityRoadTimes, threshold); //Find where the threshold belongs in the road speed array
		roadCount = innerCityRoadTimes.subList(foundItemIndex, innerCityRoadTimes.size());

		return roadCount.size(); // Return the number of roads that have speeds above the threshold
	}    
	
	/**
	 * Precondition:  arr is not null.
	 * Postcondition: The arrayList arr has been searched using an adapted binary search algorithm. 
	 * 				  The index that is returned is the first instance of the next highest road speed value from the threshold.
	 * 				  i.e. if the threshold is 7 and the array has stored the following numbers: 1,2,5,6,7,7,8,9,9 the index of the value 8 is returned. 
	 * 				  Uses an adapted binary search algorithm, so as not to include an item that is the same as the threshold. We need items strictly larger.
	 * 
	 * @param ArrayList arr - this array contains the values of all inner city road speeds. threshold - the target speed, we want all values stored in the array higher than this value. 
	 * @return left - int used in countInnerCitySlowRoads method to return part of an arrayList
	 */
	private static int binarySearch(List<Double> arr, double threshold) {
		int left = 0;
		int right = arr.size() - 1;

		while (left <= right) {
			int mid = left + (right - left) / 2;

			if (arr.get(mid) < threshold) {
				left = mid + 1; 
			} else if (arr.get(mid) > threshold) {
				right = mid - 1; 
			} else { 		//mid == threshold
				right = mid; 
				left = mid + 1; 
			}

			//The item stored at the left index is the same as the threshold, we need items greater than this. Search the right part of the array again, excluding the leftmost item.
			if(arr.get(left) == threshold) {
				left++;
				right = arr.size() - 1;
			}
		}
		return left;
	}
	
	/**
	 * Precondition:  innerCity HashMap is not empty
	 * Postcondition: getBottleNeckIterativeDFS method has been called to traverse the graph and add relevant edges to the bottleEdgeRoads set
	 * 				  A list of road names representing roads that are in the 'inner city' which, if any single one is closed, would result in an intersection no longer being reachable from the 'inner city'	
	 * 				  has been returned		  
	 * @param None
	 * @return String Array - An array of road names 
	 */
	public String[] cityBottleneckRoads() {
		Set<String> bottleEdgeRoads = new HashSet<String>();
		
		for(GraphEdge edge : innerCityEdgeMap.keySet()) { //Pass in each edge in the city
			getBottleNeckIterativeDFS(edge, bottleEdgeRoads);
		}
		
		String[] bottleNeckRoads = bottleEdgeRoads.toArray(new String[0]);
		
		return bottleNeckRoads;
	}

	/**
	 * Precondition:  innerCity HashMap is not empty
	 * Postcondition: All reachable nodes in the city has been traversed
	 * 			      The edge that is passed into the method represents a closed road. We want to see if this road is a bottleneck road
	 * 			      If we can reach every node in the city then this edge is not a bottleneck, if there are some nodes that can't be reached then 
	 * 				  it is a bottleneck and the road name is added to the bottleEdgeRoads set 				  
	 * 				 
	 * @param GraphEdge closedEdge - closed road in the city, Set<String> bottleEdgeRoads - List of road names that create a bottleneck in the city if they are closed
	 * @return None
	 */
	public void getBottleNeckIterativeDFS(GraphEdge closedEdge, Set<String> bottleEdgeRoads) {
		Stack<GraphNode> openStack = new Stack<>();
		HashSet<GraphNode> closedSet = new HashSet<>();
		
		GraphNode startNode = (GraphNode) innerCity.values().toArray()[0]; //Start with first node of the inner city
		openStack.push(startNode); //First open node 
		
		while (!openStack.isEmpty()) { 	
		GraphNode current = openStack.pop(); 	//Retrieve next node and remove from open
		closedSet.add(current);					

			for (GraphEdge edge : current.edges) {	//Check each edge for children/neighbours of current
				
				if (closedEdge != null && closedEdge.equals(edge)) { // Skip the closed edge i.e. don't traverse it
	                continue; 
	            }
				GraphNode neighbour = edge.intersection1.equals(current.name) ? graph.get(edge.intersection2) : graph.get(edge.intersection1); //Get the node that isn't the same as current node
		
				if (!closedSet.contains(neighbour) && !openStack.contains(neighbour)) { //Only if neighbour node hasn't been explored yet
					openStack.push(neighbour);    				//Node to be traversed next
				}
			}
		}
		
		Set<String> innerCityNodes = new HashSet<String>(innerCity.keySet()); //Get all node names in the inner city
		
		if(closedSet.size() != innerCityNodes.size()) { //If sizes are different then we couldn't traverse all nodes with this edge closed, add this edge to the bottle neck edge set
			bottleEdgeRoads.add(closedEdge.name);
		}
	}
	
	/**
	 * Precondition:  graph is not empty
	 * Postcondition: A String array has been returned.
	 * 				  This method has called the BFS method roadsToCloseIterativeBFS 
	 * 				  roadsToCloseIterativeBFS() uses a customised Breadth first search algorithm to return a String array containing all edge (names) that are within the specified number of hops to the intersection passed in (intersectionName).
	 *
	 * @param String intersectionName - The name of the intersection to check
	 * @param int hops - The distance from the intersection to lockdown roads
	 * @return String array - An array of road names that have endpoints within 'hops' of the intersection
	 * */
	public String[] lockdownIntersection(String intersectionName, int hops) {

		String[] roadsToClose = roadsToCloseIterativeBFS(intersectionName, hops);

		return roadsToClose;
	}
		
	/**
	 * Precondition:  startNode is stored in the graph HashMap, hops is a valid integer
	 * Postcondition: An array of edge names has been returned. 
	 * 				  These names are used in lockdownIntersection() method to determine which roads should be closed for the minister's visit.
	 * 				  This method uses an adapted iterative Breadth first search algorithm.
	 * 				  It keeps track of the level of the search tree (number of hops from start node), as well as the number of nodes in this level to be searched.
	 * 				  It also tracks the number of nodes in the next level to be searched, if we didn't know this we wouldn't know when to increment for the next level.
	 * 				  Once the desired level/hops has been reached the algorithm stops and the string array is processed and returned.
	 * 
	 * @param String startNode - The name of the starting node for the BFS algorithm. 
	 * @param int hops - The number of levels that this algorithm must iterate through before stopping.
	 * @return String[] roadsToClose - returns an array of road/edge names.
	 */
	private String[] roadsToCloseIterativeBFS(String startNode, int hops) {
		Queue<GraphNode> openQueue = new LinkedList<GraphNode>();
		Set<GraphNode> visited = new HashSet<GraphNode>();
		
		HashSet<String> edgeToClose = new HashSet<String>();
		GraphNode node = graph.get(startNode);

		openQueue.offer(node);
		
		//Keep track of levels/nodes in the graph
		int currentLevel = 0; 			//Where we are in the search tree
	    int nodesInCurrentLevel = 1;	//How many nodes we need to search on this (while) loop
	    int nodesInNextLevel = 0; 		//Start with 0 since we havn't discovered any nodes for the next level yet
	    
		while (!openQueue.isEmpty()) {

			GraphNode current = openQueue.poll();
			visited.add(current);
			nodesInCurrentLevel--; //We've traversed this node, decrement 
			
			for (GraphEdge edge : current.edges) {
				GraphNode neighbour = edge.intersection1.equals(current.name) ? graph.get(edge.intersection2) : graph.get(edge.intersection1); 
				if (!visited.contains(neighbour)) {
					openQueue.offer(neighbour); //Get all nodes connected to this one
					edgeToClose.add(edge.name); //Add all edges from this node to the edgeToClose Set, these are the roads that need to be closed for the minister's speech
					
					nodesInNextLevel++; //When a new node is discovered, add to this count so we can search it
				}
			}

			//Searched all nodes in this level, reset variables for next pass
			if (nodesInCurrentLevel == 0) {
	            nodesInCurrentLevel = nodesInNextLevel;
	            nodesInNextLevel = 0;
	            currentLevel++;		//Go up a level, since we've searched all nodes in this level
	        }
	        
	      //When we've searched to the desired level in the graph, break from the loop
	        if(currentLevel > hops) { 
				break;
			}
		}
		
		//Transfer the road names from the set to a String array and return it
		String[] roadsToClose = new String[edgeToClose.size()];

		int i = 0;
		for (String ele : edgeToClose) { 
			roadsToClose[i] = ele;
			i++;
		}

		return roadsToClose;
	}
	
	
	/**
	 * Precondition:  graph is not null
	 * Postcondition: A simulation of protests in the city has been run. The number of roads that need to be closed in order to contain the protests has been returned.
	 * 				  The city moves first then the protesters move. 
	 * 				  On the cities' first move it will identify a connected component of intersections blocked by protesters, and close the minimum number of roads to separate it from the rest of the city.
	 * 				  On the protesters move, they will spread to all connected intersections with a road open to the current intersection. They will only travel a single 'hop' per time step.
	 * 
	 * 				  More specifically, a List of protest objects is created (protests that are adjacent to each other are merged) once that list has been created, the city determines which roads to close
	 * 			      to contain a protest. It will first choose the protest with the most number of frontier nodes, then if there is a tie it will choose the protest that is easiest to contain (smallest number of outbound roads)
	 * 				  Then it is the protest's turn to spread out. This repeats until all protests have been contained. 
	 * 
	 * 				  Notes: Once a protest has nowhere else to move (all adjoining nodes closed) it will have it's status changed to contained. 
	 * 				  It is possible for a protest to self contain by spreading out to all connected nodes and therefore run out of nodes to spread to. This will affect the final road closed count.
	 * 
	 * @param String intersectionNames - An array of intersection names where protests are occurring
	 * @return int - The number of roads that need to be closed to contain the protests, according to the rules given
	 */
	public int containProtests(String[] intersectionNames){
		List<Protest> protestRegions = new ArrayList<Protest>(); //List of the city's protest regions
		Set<String> closedOutboundRoads = new HashSet<String>();
		boolean allProtestContained = false;
		
		//Step 1. Find all ‘protest regions’ - Create new protest regions from intersectionNames passed into the method
		for(String nodeName : intersectionNames) {
			boolean merged = false;
			
			GraphNode node = graph.get(nodeName);
			Protest newProtest = new Protest(nodeName);
			
			//Check if we need to merge the new Protest with an existing Protest
			if(protestRegions.size() > 0) { 
				for(int i = 0; i < protestRegions.size(); i++) {
					Protest region = protestRegions.get(i); //Get from the existing protest list
					
					if(region.checkMergeProtest(node)) { //If the node is in the protest's neigbouring node list it should be merged
						region.mergeProtest(newProtest); //Merge the protests by adding the new Protest's protestRegion, frontierRegion and roadsInRegion 
						merged = true;
					}	
				}
			}
			
			//If not merged add new Protest to the list
			if(!merged) {
				protestRegions.add(newProtest);
			}
		}
		
		//Continue the simulation until all protest regions have been contained
		while(!allProtestContained) { 
			Protest protest = null;
			
			//Step 2. Isolate the most dangerous region (the one that would spread to the most peaceful intersections on the next step)	
			int highestCount = 0;
			int index = 0;
			int lowestIndex = 0;
			int lowestCount = 0;
			boolean tie = false;
			
			for(int i = 0; i < protestRegions.size(); i++) { //Check each protest 
				
				protest = protestRegions.get(i);
				
				if(!protest.protestContained) { //protest is not already contained

					int frontierNodes = protest.frontierRegion.size();
					int currentRoadCount = protest.roadsInProtest.size();
					
					if(lowestCount == 0) { //First loop through this will be 0 - change it to protest's value
						lowestCount = currentRoadCount; 
					}
					
					if(frontierNodes > highestCount) { 	//Get index of protest containing the highest number of frontier nodes
						
						highestCount = frontierNodes; 	//Found a protest with a larger number of frontier nodes
						index = i; 
						
					}else if(frontierNodes == highestCount && i > 0) {//If there is a tie, choose the one that requires the least number of roads to be closed in order to contain it.
						tie = true;
					}
					
					if(currentRoadCount < lowestCount) { //Keep track of the protest with the lowest number of outbound roads (easiest to contain)
						lowestIndex = i;
					}			
				}
			}
			
			if(tie) { //Tie index = protest with least number of roads to be closed in order to contain it
				index = lowestIndex;
			}
			
			protest = protestRegions.get(index);
			
			//Close this protest's roads by adding to outboundRoads
			closedOutboundRoads.addAll(protest.roadsInProtest.keySet()); //close that road

			//Close outbound roads found in ALL of the Protest regions - If a protest region no longer has any outbound roads that can reach a frontier node it has been contained
			for(int i = 0; i < protestRegions.size(); i++) {
				protestRegions.get(i).closeRoads(closedOutboundRoads.toArray());
			}	
			
			//3. Spread the protest in the remaining regions outward into any non-closed roads connected to the protest group
			for(int i = 0; i < protestRegions.size(); i++) {
				Protest region = protestRegions.get(i);
				if(!region.protestContained) { //If not contained continue
					region.spreadProtest(); //All non-contained protests will spread out onto any adjacent non-closed roads by 1 hop - if this protest spreads out to another protest area and has nowhere else to move, it will be contained
				}	
			}

			allProtestContained = allProtestsContained(protestRegions); //Check if we have any protests not contained
		}
		return closedOutboundRoads.size();
	}
	
	/*
	 * Private static helper method
	 * Precondition:  List<Protest> protestRegions is not null
	 * Postcondition: A boolean indicating whether all protest's have been contained, has been returned
	 * 
	 * @param  protestRegions - List containing all protest regions in the city
	 * @return boolean indicating whether all protest regions in the city have been contained
	 * */
	private static boolean allProtestsContained(List<Protest> protestRegions) {
		boolean contained = false;
		
		for(Protest p: protestRegions) {
			if(!p.protestContained) { //Protest not contained - return false if at least one is false
				return false;
			}else {
				contained = true;
			}
		}
		return contained;
	}

	/**
	 * 	Private inner helper class Protest
	 * 
	 * 	Used to hold Protest information. 
	 *  A protest will start with a single node, if two protests are adjacent to each other they are merged. Both nodes will then be in the protest region with their neighbour nodes allocated to the frontier region.
	 *  The roads leading to the frontier roads are contained in the roadInProtest HashMap.
	 *  A protest may spread 1 hop to their adjacent nodes (frontier nodes) on their turn to move. They will spread only if they have not been contained by road closures. The city (simulation run in containProtests() method)
	 *  will try and contain the protest by closing roads. Once there is no longer any road for the protest to move through the protest is marked as contained.
	 */
	private class Protest{
		private HashSet<GraphNode> protestRegion; 	//Current area of a protest
		private HashSet<GraphNode> frontierRegion; 	//Next area to be affected (nodes 1 hop away)
		private HashMap<String, GraphNode> roadsInProtest;	//Roads leading from protest area to frontier nodes
		private List<String> closedCityRoads; 		//List of all roads closed in the city
		private boolean protestContained = false;
		
		/*
		 * Protest constructor method
		 * Precondition:  None
		 * Postcondition: A new protest region has been created with it's protest region set as the start node and it's frontier region, set to any nodes that reach that node by 1 hop.
		 * 				  roadsInProtest has been set to the contain the roads leading from the protest region to the frontier nodes. 
		 * 
		 * @param String nodeName - The name of the starting node of the protest
		 * */
		public Protest(String nodeName) {
			protestRegion = new HashSet<GraphNode>();
			frontierRegion = new HashSet<GraphNode>();
			roadsInProtest = new HashMap<String, GraphNode>();
			closedCityRoads = new ArrayList<String>();
			
			//Add the start node
			GraphNode node = graph.get(nodeName);
			protestRegion.add(graph.get(node.name)); 
			
			roadsInProtest.putAll(node.neighbourNodes);			//Link protest area to frontier nodes
			addFrontierRegion(node); 	//frontierRegion Starts with node's neighbours (within 1 hop)
		}
		
		
		/*
		 * Precondition:  node is not null, neighbourNodes is not null
		 * Postcondition: Nodes that are 1 hop away from the current node have been added to the frontier region of the protest
		 * 
		 * @param GraphNode node - The starting node of the protest
		 * @return None
		 * */
		public void addFrontierRegion(GraphNode node) {
			for(GraphNode current : node.neighbourNodes.values()) {
				frontierRegion.add(current);
			}
		}
		
		/*
		 * Precondition:  node is not null
		 * Postcondition: A boolean has been returned. If a protest group is linked to another protest group (via a shared road) the node will be contained in this protest's roadsInRegion HashMap
		 * 
		 * @param GraphNode node - The starting node of the protest
		 * @return boolean - Returns true if the protest containing the node passed into the method should be merged with this protest. False otherwise.		    
		 * */
		public boolean checkMergeProtest(GraphNode node) {
			if(roadsInProtest.containsValue(node)){
					return true;
			}
			return false;			
		}
		
		/*
		 * Precondition:  protest2 is not null
		 * Postcondition: Merges two adjacent protests into one protest. This protest and protest2's protestRegion, frontierRegion and roadsInProtest have been combined.
		 * 
		 * @param Protest protest2 - The protest object to be merged with this protest.	
		 * @return None	    
		 * */
		public void mergeProtest(Protest protest2) {
			HashSet<GraphNode> tempProtestRegion = new HashSet<GraphNode>();
			HashSet<GraphNode> tempFrontierRegion = new HashSet<GraphNode>();
			HashMap<String, GraphNode> tempRoads = new HashMap<String, GraphNode>();
			
			//Combine the protest and frontier regions
			tempProtestRegion.addAll(protest2.protestRegion);
			tempProtestRegion.addAll(this.protestRegion);
			
			tempFrontierRegion.addAll(protest2.frontierRegion);
			tempFrontierRegion.addAll(this.frontierRegion);
			
			//Remove any nodes from frontier that are already in the protest region
			for(GraphNode node: tempProtestRegion) { 
				if(tempProtestRegion.contains(node)) {
					tempFrontierRegion.remove(node);
				}
			}
			
			//Only add roads that don't connect to nodes in the protest region since we've already spread to these intersections
			for(HashMap.Entry<String, GraphNode> entry : roadsInProtest.entrySet()) {
			    String streetName = entry.getKey();
				GraphNode node = entry.getValue();
				
				if(!tempProtestRegion.contains(node)) {
					tempRoads.put(streetName, node);
				}
			}
			
			for(HashMap.Entry<String, GraphNode> entry : protest2.roadsInProtest.entrySet()) {
			    String streetName = entry.getKey();
				GraphNode node = entry.getValue();
				
				if(!tempProtestRegion.contains(node)) {
					tempRoads.put(streetName, node);
				}
			}
			
			//Overwrite existing regions/roads
			this.protestRegion = tempProtestRegion;
			this.frontierRegion= tempFrontierRegion;
			this.roadsInProtest = tempRoads;
			
		}
		
		/*
		 * Precondition:  roads is not null
		 * Postcondition: Closed roads have been removed from this protest's roadsInProtest HashMap as they can no longer be used to travel.
		 * 			      They have also been added to the closedCityRoads List (All roads closed in the city)
		 * 				  A check has been performed to see if this protest has a linked road that it can travel down for it's next move. If not then protestContained is set to true.
		 * 
		 * @param Object[] roads - List of Strings representing road names that have been closed across the city.
		 * @return None		    
		 * */
		public void closeRoads(Object[] roads) {
			for(Object roadName : roads) {
				roadsInProtest.remove(roadName);
				closedCityRoads.add((String) roadName);
			}
			
			//When there is no where else for the protest to move (all roads are closed) then the protest has been contained
			if(roadsInProtest.size() == 0) {
				protestContained = true;
			}
		}
		
		/*
		 * Precondition:  None.
		 * Postcondition: This protest has spread out by one hop through a non-closed road. protestRegion now includes all neigbouringNodes. Any previously closed roads have been removed.
		 * 				  Relevant nodes have been added to the frontier region (i.e. neighbouring nodes not already in the protest's current region)
		 * 			      Any roads that will not reach a frontier node have been removed - these roads are now part of the protest area and can't be closed/already infected
		 * 				  Roads that will reach a frontier node has been added to roadsInProtest HashMap 
		 * 				  A check has been performed to see if this protest has a linked road that it can travel down for it's next move. If not then protestContained is set to true
		 * @param None
		 * @return None
		 * */
		public void spreadProtest() {
			protestRegion.addAll(this.frontierRegion); //Add the intersections that were in the frontier to the current protestRegion
			frontierRegion = new HashSet<GraphNode>(); //Clear the HashSet ready for the next frontier of nodes
			HashMap<String, GraphNode> tempRoads = new HashMap<String, GraphNode>();
			
			//Get all neighbouring nodes from nodes that have moved into the protestRegion
			for(GraphNode current : protestRegion) { 
				tempRoads.putAll(current.neighbourNodes);
			}
			
			//Remove any closed roads from the roads list
			for(String roadName: closedCityRoads) { 
				if(tempRoads.containsKey(roadName)) {
					tempRoads.remove(roadName);
				}
			}
			
			//Only add nodes to the frontier if they're not already in the protest region
			for(GraphNode node: tempRoads.values()) { 
				if(!protestRegion.contains(node)) {
					frontierRegion.add(node);
				}
			}
			
			roadsInProtest.putAll(tempRoads); //Add all the roads list to the protest (mainly so we know how many to iterate over)
			
			//Remove any roads that will not reach a frontier node - these roads are now part of the protest area and can't be closed/already infected
			for(HashMap.Entry<String, GraphNode> entry : roadsInProtest.entrySet()) {
			    String streetName = entry.getKey();
				GraphNode node = entry.getValue();
				
				if(!frontierRegion.contains(node)) {
					tempRoads.remove(streetName);
				}
			}
			
			roadsInProtest = new HashMap<String, GraphNode>(); //Make a new list
			roadsInProtest.putAll(tempRoads); 	//Add roads that will reach a frontier node to the roads map - these are the possible roads that can be closed
				
			if(roadsInProtest.size() == 0) { //Check that the protest still has somewhere to go, if not then we've merged with another protest area and all other streets are closed already
				protestContained = true;
			}
		}
	}
	
	/**
	 * 	Private inner helper class GraphNode
	 * 
	 * 	Used to hold node (i.e. intersection) information. 
	 * 	These nodes are used in the TrafficAnalyser class to to represent a map/graph.
	 *  A node may have 1 or many edges/roads.
	 *  It also contains neighbour nodes - these are the nodes that are linked to this node by an edge (1 hop away)
	 */
	private class GraphNode {
		private String name;
		private List<GraphEdge> edges;
		private HashMap<String, GraphNode> neighbourNodes; //String = Road name linking the nodes, GraphNode = node linked to this node via the road
		
		public GraphNode(String name) {
			this.name = name;
			this.edges = new ArrayList<>();
			this.neighbourNodes = new HashMap<String, GraphNode>();
		}

		/**
		 *  Overridden method of toString.
		 *  Used when debugging to easier read node names, rather than the object name (e.g. TrafficAnalyser$GraphNode@44f75083).
		 *  
		 *  @param None.
		 *  @return String - Returns the nodes name.
		 */
		@Override
		public String toString() {
			return name;
		}
	}

	/**
	 * 	Private inner helper class GraphEdge
	 * 
	 * 	Used to hold edge (i.e. road) information. 
	 * 	These edges are used in the TrafficAnalyser class to to represent a map/graph.
	 *  An edge must have a previous and next intersection. (intersection1 and interesection2 respectively).
	 *  The average time it takes to travel the road is also stored as a double.
	 */
	private class GraphEdge {
		private String name;
		private String intersection1;
		private String intersection2;
		private double avgTravel = 0;

		public GraphEdge(String name, String intersection1, String intersection2, double averageTravel) {
			this.name = name;
			this.intersection1 = intersection1;
			this.intersection2 = intersection2;
			this.avgTravel = averageTravel;
		}

		/**
		 *  Overridden method of toString.
		 *  Used when debugging to easier read road names, rather than the object name (e.g. TrafficAnalyser$GraphEdge@44f75083).
		 *  
		 *  @param None.
		 *  @return String - Returns the nodes name.
		 */
		@Override
		public String toString() {
			return name;
		}
	}
}
