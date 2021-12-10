
import com.sun.jdi.IntegerValue;

import java.util.*;
import java.io.*;
public class DirectedGraph  {
    ArrayList<DirectedNodeList> dGraph;
    int numVertex;
   boolean [] marked;
   int edgecount = 0;
    
    
    public DirectedGraph() {
        dGraph = new ArrayList<>();
        numVertex=0;
        
    }
    
    public DirectedGraph(int n) {
      numVertex =n;
      dGraph = new ArrayList<>(n);
      marked= new boolean[n];
      for (int i=0;i<numVertex;i++)
       dGraph.add(new DirectedNodeList());
    }
    
    public void addEdge(int u, int v) {
       // assume all vertices are created
       // directed edge u to v will cause outdegree of u to go up and indegree of v to go up.
       
       if (u>=0 && u<numVertex && v>=0 && v<numVertex) { 
	      if (u!=v) {
           getNeighborList(u).addToOutList(v);
           getNeighborList(v).addToInList(u);
           edgecount++;
		  }
        }
        else throw new IndexOutOfBoundsException();
    }
    
    public DirectedNodeList getNeighborList(int u) {
        return dGraph.get(u);
    }
    
    public void printAdjacency(int u) {
       DirectedNodeList dnl = getNeighborList(u);
       System.out.println ("vertices going into "+u+"  "+dnl.getInList());
       System.out.println ("vertices going out of "+u+"  "+dnl.getOutList());
       System.out.println();
    }
    
    public void postOrderDepthFirstTraversal() {
       for (int i=0;i<numVertex;i++) 
       if (!marked[i])
           postOrderDFT (i);
       
    }
    public void postOrderDFT(int v){
        
        marked[v]=true;
        
        for (Integer u:dGraph.get(v).getOutList())
        if (!marked[u]) postOrderDFT(u);
       System.out.println(v);
    }
    
    public void depthFirstTraversal() {
       for (int i=0;i<numVertex;i++) 
       if (!marked[i])
           dFT (i);
       
    }
    public void dFT(int v){
        System.out.println(v);
        marked[v]=true;
        
        for (Integer u:dGraph.get(v).getOutList())
        if (!marked[u]) dFT(u);
       
    }

    public static DirectedGraph readAndStoreGraph(String fileName) {
        int maxV = 0;

        ArrayList<Edge> edgeSet = new ArrayList<>();

        try {
            Scanner sc = new Scanner(new File(fileName));
            String s;

            // graph.add(new ArrayList<Integer>());
            while (sc.hasNextLine()) {
                s = sc.nextLine();
                if (s.isEmpty()) continue;
                String[] line = s.split("\\s");
                // String[] line= s.split(",");
                int v1 = Integer.parseInt(line[0]);
                int v2 = Integer.parseInt(line[1]);
                int p = Math.max(v1, v2);
                if (p > maxV) maxV = p;
                edgeSet.add(new Edge(v1, v2));
            }
        } catch (FileNotFoundException e) {
        }

//numEdges = edgeSet.size();
        DirectedGraph dgtemp = new DirectedGraph(maxV+1);
        for(Edge e: edgeSet)
            dgtemp.addEdge(e.v1, e.v2);
        return dgtemp;
    }

    public ArrayList<Integer> finishingOrder(){
        ArrayList<Integer> finishingList = new ArrayList<>();
        marked = new boolean[numVertex];
        Stack<Integer> stack = new Stack<>();
        for (int i = 0; i < numVertex; i++) {
            marked[i] = false;
        }

        for(int i=0; i<numVertex; i++){
            if(!marked[i]){

                stack = new Stack<>();
                stack.add(i);

                while(!stack.isEmpty()){
                    int curr = stack.peek();
                    marked[curr] = true;

                    ArrayList<Integer> inList = this.getNeighborList(curr).getInList();

                    int nextUp = -1;
//                    boolean newVertex = false;

                    for(int a=0; a < inList.size() && nextUp < 0; a++){
                        if(!marked[inList.get(a)]) {
                            nextUp = inList.get(a);
//                            newVertex = true;
                        }
                    }

                    if(nextUp == -1){
                        finishingList.add(0, curr);
                        stack.pop();
                    }
                    else{
                        stack.add(nextUp);
                    }
                }

            }

        }
//        System.out.println(finishingList.size());
        return finishingList;
    }


    public ArrayList<ArrayList<Integer>> kosarajusAlgorithm(){
        ArrayList<ArrayList<Integer>> scc = new ArrayList<>();
        ArrayList<Integer> finishing = finishingOrder();
        Stack<Integer> stack = new Stack<>();
        marked = new boolean[numVertex];
        for (int i = 0; i < numVertex; i++) {
            marked[i] = false;
        }

        for(int i: finishing){
            if(!marked[i]){
                ArrayList<Integer> temp = new ArrayList<Integer>();

                stack = new Stack<>();
                stack.add(i);



                while(!stack.isEmpty()){
                    int curr = stack.peek();
                    marked[curr] = true;

                    ArrayList<Integer> outList = this.getNeighborList(curr).getOutList();

                    int nextUp = -1;

                    for(int a=0; a < outList.size() && nextUp < 0; a++){
                        if(!marked[outList.get(a)]) {
                            nextUp = outList.get(a);
                        }
                    }
                    if(nextUp == -1){
                        temp.add(curr);
                        stack.pop();
                    }
                    else{
                        stack.add(nextUp);
                    }
                }
                scc.add(temp);
            }
        }
//        System.out.println("\n# of strongly connected components: " + scc.size());
//        int max = Integer.MIN_VALUE;
//        ArrayList<Integer> maxAr = new ArrayList<>();
//        for(ArrayList<Integer> vertex: scc){
//            if(max < vertex.size()){
//                max = vertex.size();
//                maxAr = vertex;
//            }
//        }
//        System.out.println("\nLargest strongly connected component: " + max);
//        System.out.println(maxAr);
        return scc;
    }

    public HashMap<Integer, DirectedNodeList> reduceToMap(){
        //Create a scc list
        ArrayList<ArrayList<Integer>> scc = kosarajusAlgorithm();
        //Create a HashMap
        HashMap<Integer, DirectedNodeList> rd = new HashMap<>(scc.size());


        //For loop for all items in SCC, adding index 0 of each item as the key
        //then create a outList array and a inList array and add each outList and
        //inList of each component to the DirectedNodeList and set that as the HashMap value for it's
        //respective leader
        for(ArrayList<Integer> comp: scc){
            int key = comp.get(0);
            HashSet<Integer> inListSum = new HashSet<>();
            HashSet<Integer> outListSum = new HashSet<>();
            for(int fol: comp){
                ArrayList<Integer> inList = this.getNeighborList(fol).getInList();
                ArrayList<Integer> outList = this.getNeighborList(fol).getOutList();
                inListSum.addAll(inList);
                outListSum.addAll(outList);
            }
            DirectedNodeList newNode = new DirectedNodeList();
            newNode.inList = new ArrayList<>();
            newNode.outList = new ArrayList<>();

            for(int in: inListSum){
                newNode.addToInList(in);
            }
            for(int out: outListSum){
                newNode.addToOutList(out);
            }
            rd.put(key, newNode);
        }
        return rd;
    }

    public DirectedGraph createReducedGraph(){
        HashMap<Integer, DirectedNodeList> dnl = reduceToMap();

        //place holder so can fit the number of vertices, only contains 10559 actual values
        DirectedGraph rg = new DirectedGraph(numVertex);
        Set<Integer> keySet = dnl.keySet();
        for(int k: keySet){
            ArrayList<Integer> leaderInList = dnl.get(k).getInList();
            ArrayList<Integer> leaderOutList = dnl.get(k).getOutList();
            for(int in: leaderInList){
                if(keySet.contains(in)) {
                    rg.addEdge(in, k);
                }
            }
            for(int out: leaderOutList){
                if(keySet.contains(out)){
                    rg.addEdge(k, out);
                }
            }
        }
        rg.numVertex = dnl.size();
        return rg;
    }


    
    public static void main(String[] args) {
        DirectedGraph slashdot = readAndStoreGraph("Slashdot0902.txt");
        ArrayList<ArrayList<Integer>> StrongCC = slashdot.kosarajusAlgorithm();

        System.out.println("\nSize of Original Graph: " + slashdot.numVertex);
        System.out.println("Number of Edges in Original Graph: " + slashdot.edgecount);

        System.out.println("\nNumber of Strongly Connected Components: " + StrongCC.size());

        int maxSCC = Integer.MIN_VALUE;
        for(ArrayList<Integer> comp: StrongCC){
            if(maxSCC < comp.size()){
                maxSCC = comp.size();
            }
        }
        System.out.println("Largest Sized Strongly Connected Components: " + maxSCC);

        DirectedGraph reduced = slashdot.createReducedGraph();

        System.out.println("\nNumber of Vertices in Reduced Graph: " + reduced.numVertex);
        System.out.println("Number of Connections in Reduced Graph: " + reduced.edgecount);


        
        
        
    }
    

   
}
