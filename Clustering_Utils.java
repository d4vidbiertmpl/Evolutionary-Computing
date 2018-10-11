import java.lang.reflect.Array;
import java.util.*;
import java.util.Comparator;

class Clustering_Utils {
    /*
     *   Class that holds all functions needed for our clustering approach.
     *
     */
    private Parameters parameters;

    private Random rnd_;

    private EA_Utils ea_utils;

    public Clustering_Utils() {
        parameters = new Parameters();
        ea_utils = new EA_Utils();
        rnd_ = new Random();
        rnd_.setSeed(16438320L);
    }

    public double[][] calcElitistCluster(ArrayList<Individual> elitist, double distance_threshold) {
        /*
         * Performs hierarchical clustering with the single linkage technique to determine the number of clusters in the
         * elitist group
         *
         * @param elitist: The elitist group of the current generation
         *
         * @returns cluster_centroids: Returns the cluster of each found centroid
         *
         * */

        // perform hierarchical clustering with single linkage
        // check bottom up each hierarchy level and test if the distance threshold t is exceeded
        // The first level where this is not violated determines the found clusters

        double[][] distance_matrix = new double[parameters.elitist_size][parameters.elitist_size];

        // fill distance matrix of elitist individuals
        for (int i = 0; i < parameters.elitist_size; i++) {
            Individual current_elite = elitist.get(i);
            for (int j = 0; j < parameters.elitist_size; j++) {
                Individual other_elitist = elitist.get(j);
                if (i == j) {
                    distance_matrix[i][j] = Double.MAX_VALUE;
                } else {
                    double distance = euclideanDistance(current_elite.getValues(), other_elitist.getValues());
                    distance_matrix[i][j] = distance;
                }
            }
        }

        // print distance matrix for testing
//        for (int i = 0; i < parameters.elitist_size; i++) {
//            double[] distance_ind_i = distance_matrix[i];
//            for (int j = 0; j < parameters.elitist_size; j++) {
//                System.out.print(distance_ind_i[j]);
//                System.out.print(" ++ ");
//            }
//            System.out.println("\n------------------------------------------------------------------------------\n");
//        }

        Map<Integer, Integer> hierarchical_clusters = new HashMap<Integer, Integer>();

        // initialize clustering => each individual is its own cluster
        for (int i = 0; i < parameters.elitist_size; i++) {
            hierarchical_clusters.put(i, i);
        }

        double current_distance = 0.0;

        while (distance_threshold > current_distance) {
            int[] indices = minimumIndices(distance_matrix);
            int i = indices[0];
            int j = indices[1];

            // min value => current_distance
            current_distance = distance_matrix[i][j];

            if (distance_threshold > current_distance) {
                int corres_cluster = hierarchical_clusters.get(i);
                hierarchical_clusters.put(j, corres_cluster);

                distance_matrix[i][j] = Double.MAX_VALUE;
                distance_matrix[j][i] = Double.MAX_VALUE;
            }
        }

//        for (int i = 0; i < hierarchical_clusters.size(); i++) {
//            int cluster_no = hierarchical_clusters.get(i);
//            System.out.print(i);
//            System.out.print(" MAPPES TO CLUSTER ");
//            System.out.print(cluster_no);
//            System.out.println("-----------------------------------");
//        }

        // convert integer cluster representation to individuals, to calculate centroid and radius
        Map<Integer, ArrayList<Individual>> cluster_to_ind = new HashMap<Integer, ArrayList<Individual>>();
        int current_cluster_count = 0;
        for (int i = 0; i < parameters.elitist_size; i++) {
            int no_current_cluster = hierarchical_clusters.get(i);
            if (cluster_to_ind.containsKey(no_current_cluster)) {
                cluster_to_ind.get(no_current_cluster).add(elitist.get(i));
            } else {
                ArrayList<Individual> new_cluster = new ArrayList<Individual>();
                new_cluster.add(elitist.get(i));
                cluster_to_ind.put(current_cluster_count, new_cluster);
                current_cluster_count++;
            }
        }

//        for (int i = 0; i < cluster_to_ind.size(); i++) {
//            if (cluster_to_ind.containsKey(i)) {
//                ArrayList<Individual> ind_list = cluster_to_ind.get(i);
//                System.out.println(i);
//                System.out.println(ind_list.size());
//                for (int j = 0; j < sind_list.size(); j++) {
//                    printIndividual(ind_list.get(j));
//                }
//                System.out.println("\n+++++++++++++++++++++++++++++++++++++++++++++++++\n");
//            }
//        }

        // Calculate centroid of each cluster
        double[][] cluster_centroids = new double[cluster_to_ind.size()][parameters.individual_size];
        for (int i = 0; i < cluster_to_ind.size(); i++) {
            ArrayList<Individual> cluster = cluster_to_ind.get(i);
            double[] centroid = calcCentroid(cluster);
            cluster_centroids[i] = centroid;
        }
        return cluster_centroids;
    }


    public ArrayList<Individual> generateRandomProletarians(double[][] cluster, boolean adaptive) {
        /*
            @param cluster: the cluster of the elitist
            @return retun: random inittialized individuals outside the clusters determined by the individuals of the elitist group.
         */
        ArrayList<Individual> proletarians = new ArrayList<Individual>();

        int current_proletarians = 0;
        while (current_proletarians < parameters.proletarian_size) {
            Individual proletarian;
            if (adaptive) {
                proletarian = ea_utils.createRandomIndividual(true);
            } else {
                proletarian = ea_utils.createRandomIndividual(false);
            }
            if (!isIndividualInCluster(cluster, proletarian)) {
                proletarians.add(proletarian);
                current_proletarians++;
            }
        }
        return proletarians;
    }


    public ArrayList<Individual> getElitistGroup(ArrayList<Individual> population) {
        /*
         * Finds the elitist group of the population.
         *
         * @param population: the current population
         *
         * @return elitist: the elitist group
         *
         * */
        ArrayList<Individual> elitist = new ArrayList<Individual>();
        Collections.sort(population);
        for (int i = 0; i < parameters.elitist_size; i++) {
            elitist.add(population.get(i));
        }
        return elitist;
    }

    // --------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------
    private boolean isIndividualInCluster(double[][] clusters, Individual individual) {
        /*
         * Tests if an individual is in the spatial range of a cluster.
         *
         * @param clusters:
         * @param individual:
         *
         * @return: true if individual in spatial range of cluster, else false
         * */
        for (int i = 0; i < clusters.length; i++) {
            double[] current_cluster = clusters[i];
            double distance = euclideanDistance(current_cluster, individual.getValues());

//            System.out.println(distance);

            if (distance < parameters.cluster_distance_thresh) {
//                System.out.println("In Cluster");
                return true;
            }
        }
//        System.out.println("Not in Cluster");
        return false;
    }


    public double[] calcCentroid(ArrayList<Individual> cluster) {
        /*
         * Calculates centroid of one cluster
         *
         * @param: Individuals in an ArrayList representing a cluster.
         * */
        double[] sum_centroid = new double[parameters.individual_size];

        for (int i = 0; i < cluster.size(); i++) {
            double[] ind_vec = cluster.get(i).getValues();

            for (int j = 0; j < parameters.individual_size; j++) {
                double dim_value = ind_vec[j];
                sum_centroid[j] += dim_value;
            }
        }

        double[] centroid = new double[parameters.individual_size];
        for (int i = 0; i < parameters.individual_size; i++) {
            centroid[i] = sum_centroid[i] / cluster.size();
        }
        return centroid;
    }

    private int[] minimumIndices(double[][] distance_matrix) {
        /*
         * Finds the minimum value in a two dimensional double array.
         *
         * @param distance_matrix: distance matrix representing the euclidean distance between the individuals.
         *
         * @return double[][]: the coordinates (i, j) of the minimum in the table.
         * */
        int min_i = 0;
        int min_j = 0;

        double min_value = Double.MAX_VALUE;

        for (int i = 0; i < parameters.elitist_size; i++) {
            for (int j = 0; j < parameters.elitist_size; j++) {
                double current_value = distance_matrix[i][j];
                if (current_value < min_value) {
                    min_value = current_value;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        return new int[]{min_i, min_j};
    }

    public double euclideanDistance(double[] a, double[] b) {
        /*
         * Calculates the euclidean distance between two individuals
         *
         * @param a, b: the position of an individual in the search space
         *
         * @return: The euclidean distance
         * */
        double sum = 0;
        for (int i = 0; i < parameters.individual_size; i++) {
            sum = sum + Math.pow((a[i] - b[i]), 2.0);
        }
        return Math.sqrt(sum);
    }

}
