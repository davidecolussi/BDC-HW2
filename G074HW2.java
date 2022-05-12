
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;


import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;


public class G074HW2 {

    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    // Input reading methods
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static Vector strToVector(String str) {
        String[] tokens = str.split(",");
        double[] data = new double[tokens.length];
        for (int i=0; i<tokens.length; i++) {
            data[i] = Double.parseDouble(tokens[i]);
        }
        return Vectors.dense(data);
    }

    public static ArrayList<Vector> readVectorsSeq(String filename) throws IOException {
        if (Files.isDirectory(Paths.get(filename))) {
            throw new IllegalArgumentException("readVectorsSeq is meant to read a single file.");
        }
        ArrayList<Vector> result = new ArrayList<>();
        Files.lines(Paths.get(filename))
                .map(str -> strToVector(str))
                .forEach(e -> result.add(e));
        return result;
    }

    /**
     * DUMMY SOLUTION FOR SEQWEIGHTOUTLIERS
     * @param P
     * @param W
     * @param k
     * @param z
     * @param alpha
     * @return
     */
    public static ArrayList<Vector> DummySeqWeightedOutliers(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {

        System.out.println(P);
        //dummy example
        //dummy fixed solution just for testing the objective function in main:
        ArrayList<Vector> dummySolution = new ArrayList<>();
        double[] dummyCenters = new double[2];
        dummyCenters[0] = 1.0;
        dummyCenters[1] = 1.0;
        Vector dummyCenter = Vectors.dense(dummyCenters);
        dummySolution.add(dummyCenter);

        double[] dummyCenters2 = new double[2];
        dummyCenters2[0] = 2.0;
        dummyCenters2[1] = 2.0;
        dummyCenter = Vectors.dense(dummyCenters2);
        dummySolution.add(dummyCenter);

        double[] dummyCenters3 = new double[2];
        dummyCenters3[0] = 3.0;
        dummyCenters3[1] = 3.0;
        dummyCenter = Vectors.dense(dummyCenters3);
        dummySolution.add(dummyCenter);
        return dummySolution;

    }

    /**
     * This method implements the weighted variant of KcenterOUT algorithm
     * @param P the set of points
     * @param W the set of weights
     * @param k the number of centers
     * @param z the number of outliers
     * @param alpha coefficient used by the algorithm
     * @return a solution of the K-center with z outliers problem
     */
    public static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {


        ArrayList<Vector> Z = new ArrayList<>();
        ArrayList<Vector> S;
        double max = 0, dist, r=-1, Wz, ballWeight;
        Vector newCenter = null, y;

        System.out.println("Calculating r as the min distance between the first k+z+1="+(k+z+1)+" points");
        double minDist = Math.sqrt(Vectors.sqdist(P.get(0), P.get(1))); // set the distance
        //calculate minDistance from the firt k+z+1 point
        for(int i = 0; i <(k+z+1); i++){
            for (int j = 0; j<(k+z+1); j++){
                if(i!=j){ //they must not be the same point
                    dist = Math.sqrt(Vectors.sqdist(P.get(i),P.get(j)));
                    //System.out.println("Distance between " + P.get(i) + " and " + P.get(j) + " is " + dist);
                    if(dist<minDist){
                        minDist = dist;
                    }
                }
            }
        }
        r = (minDist)/2;
        System.out.println("Calculated r = " + r);



        while(true) {
            System.out.println("outer while");
            for (Vector originalPoint : P) { //deep copy of P into Z (notice that the reference of each originalPoint is copied)
                System.out.println("copying " + originalPoint + "from P into Z");
                Z.add(originalPoint);
            }

            S = new ArrayList<>(); //clear S
            Wz = 0;
            for (double w : W) { // Calculate the sum of all the weight vector's weights
                System.out.println("summing to Wz the weight: " + w);
                Wz = Wz + w;
            }
            System.out.println("\nFinal Wz = " + Wz);


            while ((S.size() < k) && (Wz > 0)) {
                System.out.println("inner while [|S|=" + S.size() + ", Wz=" + Wz + "]");
                max = 0;
                newCenter = null; //todo: remove line

                //We have to compute the sum of all the weights of the points which are in the ball Bz(x,(1+2a)r)
                for (Vector x : P) { //for each point x in P (ball center)
                    System.out.println("Considering point of P x: " + x);
                    ballWeight = 0;
                    //calculate ballWeight by summing all weights of the points in the ball Bz with center x and radius (1+2a)r.
                    //a point y nof the set Z is considered in the ball Bz if its distance from the ball center x is less or equal than (1+2a)r
                    for (int i = 0; i < Z.size(); i++) { //for each point y in Z
                        y = Z.get(i);
                        System.out.println("considering point of Z y: " + y);
                        dist = Vectors.sqdist(x, y);
                        System.out.println("calculated distance of x to y " + dist);

                        if (dist <= (1 + 2 * alpha) * r) { //check if the point y is in the ball Bz(x,(1+2a)r)
                            //y is in the ball
                            System.out.println("The point y is in the ball. Adding its weight (" + W.get(i) + ")to the partial ballWeight sum");
                            ballWeight += W.get(i); //sum up the weight of y
                        }
                    }
                    System.out.println("\n ballWeight: " + ballWeight);

                    if (ballWeight > max) { //if i found a bigger ballWaight than the max found until now
                        max = ballWeight; //update max with current ballWeight
                        newCenter = x; //update newCenter with current x
                        System.out.println("*****NEW CENTER FOUND: " + newCenter + "******");
                    }

                }//end of for each point x in P



                if (newCenter != null) {
                    //todo: comment all this IF body for debugging
                    /*
                    S.add(newCenter); //add the newCenter found into the solution set S //todo: this line gives ERROR in run

                    System.out.println("*Adding newCenter (" + newCenter + ") to S\nCreating the ball Bz with center newCenter and radius (3+4a)r");

                    //find a new ball Bz with center "newCenter" and radius (3+4a)r.
                    //a point y of the set Z is considered in the ball Bz if his distance from the ball center "newCenter" is less or equal than (3+4a)r

                    for (int i = 0; i < Z.size(); i++) { //for each point y in Z
                        y = Z.get(i);
                        System.out.println("considering point of Z y: " + y);
                        dist = Vectors.sqdist(newCenter, y);
                        System.out.println("calculated distance of newCenter to y " + dist);

                        if (dist <= (3 + 4 * alpha) * r) { //check if the point y is in the ball Bz(newCenter,(1+2a)r)
                            System.out.println("The point y is in the ball Bz(newCenter,(3+4a)r)");
                            //remove y from Z
                            Z.remove(i);
                            Wz -= W.get(i);
                        }
                    }//end of for each point y in Z


                    if (Wz <= z) {
                        System.out.println("end of the algorithm, returning: " + S);
                        return S;
                    } else {
                        r = 2 * r;
                    }
                    //todo: comment all the IF body until here for debugging
                     */

                }//end of if(newCenter!=null)
            }//end of inner while
        }//end of outer while


    }

    /**
     * This method computes the value of the objective function of the K-center with z outliers problem.
     * Obj(P-Zs,S) = max(d(p,S)) for each p in P-Zs, where Zs is the largest set of points from S of total weight <= z
     * @param P the set of points
     * @param S the solution set
     * @param z the number of outliers
     * @return the value of the objective function
     */
    public static double ComputeObjective(ArrayList<Vector> P, ArrayList<Vector> S, int z) {
        double result = -1;
        double aux = -1;
        ArrayList<Double> distancesVector = new ArrayList<>();
        //compute a vector of distances of all point x in P to all points c in S, multiplied by weight vector

        System.out.println("Calculating objective function for the solution given by the centers: " + S);
        for (Vector point: P) { //for each point of the input set P
            //take the point
            System.out.println("Considering point: " + point);
            //compute all the distances (L2) from the point to every center c in S
            ArrayList<Double> distancesFromCenters = new ArrayList<>();


            for(Vector c: S){
                aux = Math.sqrt(Vectors.sqdist(point,c)); //L2 distance
                System.out.println("Distance from point " + point + " to center " + c + " is: " + aux);
                distancesFromCenters.add(aux); //add the L2 distance of point to c to the vector
            }
            //take only the minimum (which is the definition of d(point,S))
            double distanceFromS = Collections.min(distancesFromCenters);
            System.out.println("Distance from S is: "+ distanceFromS);
            //add the result to the distances Vector
            distancesVector.add(distanceFromS);

        }

        System.out.println("\n\ndistancesVector size: " + distancesVector.size() + " and content is:\n" + distancesVector);
        //sort the distances vector in an ascending way
        distancesVector.sort((d1,d2) -> Double.compare(d1,d2));
        System.out.println("\n\ndistancesVector SORTED:\n" + distancesVector);

        //exclude the last z elements in the sorted vector
        for(int i=0;i<z;i++){
            distancesVector.remove(distancesVector.size()-1);
        }

        System.out.println("\n\ndistancesVector WITH LAST " + z + " ELEMENTS REMOVED:\n" + distancesVector);

        //take the last element (which will represent the value of the objective function)
        result = distancesVector.get(distancesVector.size()-1);

        return result;
    }

    public static void main(String[] args) throws IOException {
        ArrayList<Vector> inputPoints = readVectorsSeq(args[0]); //File reading
        int k = Integer.parseInt(args[1]); //Number of centers
        int z = Integer.parseInt(args[2]); //Number of allowed outliers

        //create the vector of weight (for homework 2 all 1s)
        ArrayList<Long> weights = new ArrayList<Long>();;
        for(int i = 0; i<inputPoints.size(); i++) {
            weights.add(1L);
        }

        long startTime = System.currentTimeMillis();
        //ArrayList<Vector> solution = DummySeqWeightedOutliers(inputPoints,weights,k,z,0);
        ArrayList<Vector> solution = SeqWeightedOutliers(inputPoints,weights,k,z,0);


        long execTime = System.currentTimeMillis() - startTime;
        Double objective = null;
        //objective = ComputeObjective(inputPoints,solution,z);

        //Output printing
        System.out.println("Input size n = " + inputPoints.size());
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);
        /*System.out.println("Initial guess = ");
        System.out.println("Final guess = ");
        System.out.println("Number of guesses = ");
        System.out.println("Objective function = " + objective);
        System.out.println("Time of SeqWeightedOutliers = " + execTime);*/

    }
}