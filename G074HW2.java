
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

    public  static ArrayList<Vector> SeqWeightedOutliers(ArrayList<Vector> P, ArrayList<Long> W, int k, int z, int alpha) {
        return null;
    }

    public static double ComputeObjective(ArrayList<Vector> P, ArrayList<Vector> S, int z) {
        return 0;
    }

    public static void main(String[] args) throws IOException {

        ArrayList<Vector> inputPoints = readVectorsSeq(args[0]); //File reading
        int k = Integer.parseInt(args[1]); //Number of centers
        int z = Integer.parseInt(args[2]); //Number of allowed outliers

        ArrayList<Long> weights = new ArrayList<Long>();;
        for(int i = 0; i<inputPoints.size(); i++) {
            weights.add(1L);
        }

        ArrayList<Vector> solution = SeqWeightedOutliers(inputPoints,weights,k,z,0);
        double objective = ComputeObjective(inputPoints,solution,z);

    }
}