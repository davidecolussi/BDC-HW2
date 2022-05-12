
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;


import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;




public class G074HW2 {

    public static double initialGuess = 0;
    public static double finalGuess = 0;
    public static double guess = 1;


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

        double minDist = Math.sqrt(Vectors.sqdist(P.get(0), P.get(1))); // set the distance
        //calculate minDistance from the firt k+z+1 point
        for(int i = 0; i <(k+z+1); i++){
            for (int j = 0; j<(k+z+1); j++){
                if(i!=j){ //they must not be the same point
                    dist = Math.sqrt(Vectors.sqdist(P.get(i),P.get(j)));
                    if(dist<minDist){
                        minDist = dist;
                    }
                }
            }
        }
        r = (minDist)/2;




        initialGuess = r;
        finalGuess = r;

        while(true) {
            Z = new ArrayList<>();
            for (Vector originalPoint : P) { //deep copy of P into Z (notice that the reference of each originalPoint is copied)
                Z.add(originalPoint);
            }

            S = new ArrayList<>(); //clear S
            Wz = 0;
            for (double w : W) { // Calculate the sum of all the weight vector's weights
                Wz = Wz + w;
            }


            while ((S.size() < k) && (Wz > 0)) {
                max = 0;

                //We have to compute the sum of all the weights of the points which are in the ball Bz(x,(1+2a)r)
                for (Vector x : P) { //for each point x in P (ball center)
                    ballWeight = 0;
                    //calculate ballWeight by summing all weights of the points in the ball Bz with center x and radius (1+2a)r.
                    //a point y nof the set Z is considered in the ball Bz if its distance from the ball center x is less or equal than (1+2a)r
                    for (int i = 0; i < Z.size(); i++) { //for each point y in Z
                        y = Z.get(i);
                        dist = Math.sqrt(Vectors.sqdist(x, y));

                        if (dist <= (1 + 2 * alpha) * r) { //check if the point y is in the ball Bz(x,(1+2a)r)
                            //y is in the ball
                            ballWeight += W.get(i); //sum up the weight of y
                        }
                    }


                    if (ballWeight > max) { //if i found a bigger ballWaight than the max found until now
                        max = ballWeight; //update max with current ballWeight
                        newCenter = x; //update newCenter with current x
                    }


                }//end of for each point x in P



                if (newCenter != null) {


                    S.add(newCenter); //add the newCenter found into the solution set S


                    //find a new ball Bz with center "newCenter" and radius (3+4a)r.
                    //a point y of the set Z is considered in the ball Bz if his distance from the ball center "newCenter" is less or equal than (3+4a)r

                    for (int i = 0; i < Z.size(); i++) { //for each point y in Z
                        y = Z.get(i);
                        dist = Math.sqrt(Vectors.sqdist(newCenter, y));

                        if ((dist <= (3 + 4 * alpha) * r)) { //check if the point y is in the ball Bz(newCenter,(1+2a)r)
                            //remove y from Z
                            Z.remove(i);
                            Wz -= W.get(i);
                            i--;
                        }
                    }//end of for each point y in Z
                }//end of if(newCenter!=null)


            }//end of inner while
            if (Wz <= z) {
                return S;
            } else {
                r = 2 * r;
                finalGuess = r;
            }

            guess++;

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

        for (Vector point: P) { //for each point of the input set P
            //take the point
            //compute all the distances (L2) from the point to every center c in S
            ArrayList<Double> distancesFromCenters = new ArrayList<>();


            for(Vector c: S){
                aux = Math.sqrt(Vectors.sqdist(point,c)); //L2 distance
                distancesFromCenters.add(aux); //add the L2 distance of point to c to the vector
            }
            //take only the minimum (which is the definition of d(point,S))
            double distanceFromS = Collections.min(distancesFromCenters);
            //add the result to the distances Vector
            distancesVector.add(distanceFromS);

        }

        //sort the distances vector in an ascending way
        distancesVector.sort((d1,d2) -> Double.compare(d1,d2));

        //exclude the last z elements in the sorted vector
        for(int i=0;i<z;i++){
            distancesVector.remove(distancesVector.size()-1);
        }


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
        ArrayList<Vector> solution = SeqWeightedOutliers(inputPoints,weights,k,z,0);


        long execTime = System.currentTimeMillis() - startTime;
        Double objective = null;
        objective = ComputeObjective(inputPoints,solution,z);

        //Output printing
        System.out.println("Input size n = " + inputPoints.size());
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);
        System.out.println("Initial guess = "+ initialGuess);
        System.out.println("Final guess = " + finalGuess);
        System.out.println("Number of guesses = "+ guess);
        System.out.println("Objective function = " + objective);
        System.out.println("Time of SeqWeightedOutliers = " + execTime);

    }
}