
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
     * This method implements the weighted variant of KcenterOUT algorithm
     * @param P the set of points
     * @param W the set of weights
     * @param k the number of centers
     * @param z the number of outliers
     * @param alpha coefficient used by the algorithm
     * @return a solution of the K-center with z outliers problem
     */
    public static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {

        double minDist = Math.sqrt(Vectors.sqdist(P.get(0), P.get(1))); // set the distance
        //calculate minDistance from the firt k+z+1 point
        for(int i = 0; i <(k+z+1); i++){
            for (int j = 0; j<(k+z+1); j++){
                if(true) {
                    double dist = Math.sqrt(Vectors.sqdist(P.get(i),P.get(j)));
                    if(dist<minDist){
                        minDist = dist;
                    }
                }
            }
        }
        double r = (minDist)/2;

        while(true){
            ArrayList<Vector> Z = P;
            ArrayList<Vector> S = new ArrayList<>();
            double Wz = 0; // Calculate the sum of all weight
            for (double w: W) {
                Wz = Wz + w;
            }

            while(( S.size() < k )&&( Wz > 0)){
                double max = 0;
                Vector newCenter = null;
                for( int i = 0; i<P.size(); i++){

                    ArrayList<Vector> ball = new ArrayList<>();
                    double ballWeight = 0;
                    for( int j = 0; j< P.size(); j++){
                        if(i!=j){
                            double dist = Math.sqrt(Vectors.sqdist(P.get(i), P.get(j)));
                            if(dist<(1+2*alpha)*r){
                                ball.add(P.get(j));
                                ballWeight += W.get(j);
                            }
                        }

                    }

                    if(ballWeight > max){
                        max = ballWeight;
                        newCenter = P.get(i);
                    }
                }


                if(newCenter != null) {
                    S.add(newCenter);
                    ArrayList<Vector> Bz = new ArrayList<>();
                    for (int j = 0; j < P.size(); j++) {

                        double dist = Math.sqrt(Vectors.sqdist(P.get(j), newCenter));
                        if (dist < (3 + 4 * alpha) * r) {
                            Z.remove(P.get(j));
                            Wz -= W.get(j);
                        }
                    }
                }

            }
            if(Wz < z){
                System.out.println(S.toArray());
                return S;
            }else {
                r = 2*r;
            }

        }




/*

        System.out.println(P);
        //todo
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

        */
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
        ArrayList<Vector> solution = SeqWeightedOutliers(inputPoints,weights,k,z,0);

        long execTime = System.currentTimeMillis() - startTime;

        double objective = ComputeObjective(inputPoints,solution,z);

        //Output printing
        System.out.println("Input size n = " + inputPoints.size());
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);
        System.out.println("Initial guess = ");
        System.out.println("Final guess = ");
        System.out.println("Number of guesses = ");
        System.out.println("Objective function = " + objective);
        System.out.println("Time of SeqWeightedOutliers = " + execTime);

    }
}