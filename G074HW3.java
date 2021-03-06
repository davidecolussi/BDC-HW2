
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaDoubleRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.BLAS;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.codehaus.janino.Java;
import scala.Array;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;

public class G074HW3
{

    public static double initialGuess = 0;
    public static double finalGuess = 0;
    public static double guess = 1;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// MAIN PROGRAM
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static void main(String[] args) throws Exception {

        if (args.length != 4) {
            throw new IllegalArgumentException("USAGE: filepath k z L");
        }

        // ----- Initialize variables
        String filename = args[0];
        int k = Integer.parseInt(args[1]);
        int z = Integer.parseInt(args[2]);
        int L = Integer.parseInt(args[3]);
        long start, end; // variables for time measurements

        // ----- Set Spark Configuration
        Logger.getLogger("org").setLevel(Level.OFF);
        Logger.getLogger("akka").setLevel(Level.OFF);
        SparkConf conf = new SparkConf(true).setAppName("MR k-center with outliers");
        JavaSparkContext sc = new JavaSparkContext(conf);
        sc.setLogLevel("WARN");

        // ----- Read points from file
        start = System.currentTimeMillis();
        JavaRDD<Vector> inputPoints = sc.textFile(args[0], L)
                .map(x-> strToVector(x))
                .repartition(L)
                .cache();
        long N = inputPoints.count();
        end = System.currentTimeMillis();

        // ----- Pring input parameters
        System.out.println("File : " + filename);
        System.out.println("Number of points N = " + N);
        System.out.println("Number of centers k = " + k);
        System.out.println("Number of outliers z = " + z);

        System.out.println("Number of partitions L = " + L);
        System.out.println("Time to read from file: " + (end-start) + " ms");

        // ---- Solve the problem
        ArrayList<Vector> solution = MR_kCenterOutliers(inputPoints, k, z, L);


        // ---- Compute the value of the objective function
        start = System.currentTimeMillis();
        double objective = computeObjective(inputPoints, solution, z, sc);
        end = System.currentTimeMillis();
        System.out.println("Initial guess = "+ initialGuess);
        System.out.println("Final guess = " + finalGuess);
        System.out.println("Number of guesses = "+ guess);
        System.out.println("Objective function = " + objective);
        System.out.println("Time to compute objective function: " + (end-start) + " ms");

    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// AUXILIARY METHODS
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method strToVector: input reading
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static Vector strToVector(String str) {
        String[] tokens = str.split(",");
        double[] data = new double[tokens.length];
        for (int i=0; i<tokens.length; i++) {
            data[i] = Double.parseDouble(tokens[i]);
        }
        return Vectors.dense(data);
    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method euclidean: distance function
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static double euclidean(Vector a, Vector b) {
        return Math.sqrt(Vectors.sqdist(a, b));
    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method MR_kCenterOutliers: MR algorithm for k-center with outliers
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static ArrayList<Vector> MR_kCenterOutliers (JavaRDD<Vector> points, int k, int z, int L)
    {
        long start, end;
        start = System.currentTimeMillis();

        //------------- ROUND 1 ---------------------------

        JavaRDD<Tuple2<Vector,Long>> coreset = points.mapPartitions(x ->
        {
            ArrayList<Vector> partition = new ArrayList<>();
            while (x.hasNext()) partition.add(x.next());
            ArrayList<Vector> centers = kCenterFFT(partition, k+z+1);
            ArrayList<Long> weights = computeWeights(partition, centers);
            ArrayList<Tuple2<Vector,Long>> c_w = new ArrayList<>();
            for(int i =0; i < centers.size(); ++i)
            {
                Tuple2<Vector, Long> entry = new Tuple2<>(centers.get(i), weights.get(i));
                c_w.add(i,entry);
            }
            return c_w.iterator();
        }); // END OF ROUND 1


        //------------- ROUND 2 ---------------------------

        ArrayList<Tuple2<Vector, Long>> elems = new ArrayList<>((k+z)*L);
        elems.addAll(coreset.collect());

        end = System.currentTimeMillis();

        System.out.println("Time Round 1: "+ (end-start) + " ms");


        start = System.currentTimeMillis();
        //
        // ****** ADD YOUR CODE
        // ****** Compute the final solution (run SeqWeightedOutliers with alpha=2)
        // ****** Measure and print times taken by Round 1 and Round 2, separately
        // ****** Return the final solution
        //

        ArrayList<Vector> P = new ArrayList<>();
        ArrayList<Long> W = new ArrayList<>();

        for (Tuple2<Vector,Long> elem : elems) {
            P.add(elem._1);
            W.add(elem._2);
        }

        ArrayList<Vector> S = SeqWeightedOutliers1(P,W,k,z,2);

        end = System.currentTimeMillis();
        System.out.println("Time Round 2: "+ (end-start) + " ms");

        return S;
    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method kCenterFFT: Farthest-First Traversal
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static ArrayList<Vector> kCenterFFT (ArrayList<Vector> points, int k) {

        final int n = points.size();
        double[] minDistances = new double[n];
        Arrays.fill(minDistances, Double.POSITIVE_INFINITY);

        ArrayList<Vector> centers = new ArrayList<>(k);

        Vector lastCenter = points.get(0);
        centers.add(lastCenter);
        double radius =0;

        for (int iter=1; iter<k; iter++) {
            int maxIdx = 0;
            double maxDist = 0;

            for (int i = 0; i < n; i++) {
                double d = euclidean(points.get(i), lastCenter);
                if (d < minDistances[i]) {
                    minDistances[i] = d;
                }

                if (minDistances[i] > maxDist) {
                    maxDist = minDistances[i];
                    maxIdx = i;
                }
            }

            lastCenter = points.get(maxIdx);
            centers.add(lastCenter);
        }
        return centers;
    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method computeWeights: compute weights of coreset points
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static ArrayList<Long> computeWeights(ArrayList<Vector> points, ArrayList<Vector> centers)
    {
        Long weights[] = new Long[centers.size()];
        Arrays.fill(weights, 0L);
        for(int i = 0; i < points.size(); ++i)
        {
            double tmp = euclidean(points.get(i), centers.get(0));
            int mycenter = 0;
            for(int j = 1; j < centers.size(); ++j)
            {
                if(euclidean(points.get(i),centers.get(j)) < tmp)
                {
                    mycenter = j;
                    tmp = euclidean(points.get(i), centers.get(j));
                }
            }
            // System.out.println("Point = " + points.get(i) + " Center = " + centers.get(mycenter));
            weights[mycenter] += 1L;
        }
        ArrayList<Long> fin_weights = new ArrayList<>(Arrays.asList(weights));
        return fin_weights;
    }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method SeqWeightedOutliers: sequential k-center with outliers
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    /**
     * This method implements the weighted variant of KcenterOUT algorithm
     * @param P the set of points
     * @param W the set of weights
     * @param k the number of centers
     * @param z the number of outliers
     * @param alpha coefficient used by the algorithm
     * @return a solution of the K-center with z outliers problem
     */
    public static ArrayList<Vector> SeqWeightedOutliers1(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {


        ArrayList<Vector> Z = new ArrayList<>();
        ArrayList<Long> Weight = new ArrayList<>();

        ArrayList<Vector> S;
        double max = 0, dist, r=-1;
        long Wz, ballWeight;
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
            Weight = new ArrayList<>();

            for (Vector originalPoint : P) { //deep copy of P into Z (notice that the reference of each originalPoint is copied)
                Z.add(originalPoint);
            }

            for (long weight : W) { //deep copy of P into Z (notice that the reference of each originalPoint is copied)
                Weight.add(weight);
            }


            S = new ArrayList<>(); //clear S
            Wz = 0;
            for (long w : Weight) { // Calculate the sum of all the weight vector's weights
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

                        if (dist <= ((1 + (2 * alpha)) * r)) { //check if the point y is in the ball Bz(x,(1+2a)r)
                            //y is in the ball
                            ballWeight += Weight.get(i); //sum up the weight of y
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

                        if ((dist <= ((3 + (4 * alpha)) * r))) { //check if the point y is in the ball Bz(newCenter,(1+2a)r)
                            //remove y from Z
                            Z.remove(i);
                            Wz -= Weight.get(i);
                            Weight.remove(i);
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


    public static ArrayList<Vector> SeqWeightedOutliers2(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {



        ArrayList<Integer> indexs = new ArrayList<>();
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

            indexs = new ArrayList<>();

            for(int i = 0; i<P.size(); i++){
                indexs.add(i);
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
                    for (int i = 0; i < indexs.size(); i++) { //for each point y in Z
                        y = P.get(indexs.get(i));
                        dist = Math.sqrt(Vectors.sqdist(x, y));

                        if (dist <= (1 + 2 * alpha) * r) { //check if the point y is in the ball Bz(x,(1+2a)r)
                            //y is in the ball
                            ballWeight += W.get(indexs.get(i)); //sum up the weight of y
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

                    for (int i = 0; i < indexs.size(); i++) { //for each point y in Z
                        y = P.get(indexs.get(i));
                        dist = Math.sqrt(Vectors.sqdist(newCenter, y));

                        if ((dist <= (3 + 4 * alpha) * r)) { //check if the point y is in the ball Bz(newCenter,(1+2a)r)
                            //remove y from Z

                            Wz -= W.get(indexs.get(i));
                            indexs.remove(i);
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

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Method computeObjective: computes objective function
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    public static double computeObjective(JavaRDD<Vector> points, ArrayList<Vector> centers, int z, JavaSparkContext sc)
    {

        //
        // ****** ADD THE CODE FOR computeObjective
        //

        List <Double> highestDistances = points.mapPartitions(element -> {
            ArrayList <Double> distancesList = new ArrayList<>();
            while(element.hasNext()) {
                Vector point = element.next();
                Double min = null;
                for(Vector center : centers) {
                    double distance = euclidean(center,point);
                    if(min == null || distance<min)
                        min = distance;

                }
                distancesList.add(min);
            }
            return distancesList.iterator();
        }).top(z+1);

        return highestDistances.get(highestDistances.size()-1);
    }
}