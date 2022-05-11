
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
    public  static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {
        return null;
    }

    /**
     * This method computes the value of the objective function of the K-center with z outliers problem
     * @param P the set of points
     * @param S the solution set
     * @param z the number of outliers
     * @return the value of the objective function
     */
    public static double ComputeObjective(ArrayList<Vector> P, ArrayList<Vector> S, int z) {
        return 0;
    }

    public static void main(String[] args) throws IOException {
        System.out.println("Hello world");
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