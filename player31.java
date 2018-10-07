import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import javax.naming.directory.InitialDirContext;
import java.util.*;
import java.util.Comparator;

public class player31 implements ContestSubmission {

    // individuals are always 10-dim vectors
    private static final int individual_size_ = 10;
    // all single values of the individuals are in [-5,5]
    private static final double values_min_ = -5.0;
    private static final double values_max_ = 5.0;

    // object for creating random numbers
    Random rnd_;

    // object for evaluating created individuals
    ContestEvaluation evaluation_;

    private int evaluations_limit_;   // number of allowed evaluations (defined by framework)
    private int evaluations_counter_ = 0; // counter for performed evaluations, starts at 0

    private int population_size_ = 100; // size of one population
    private int elitist_size_ = 10;
    private int proletarian_size_ = 20;

    private double cluster_distance_thresh_ = 2.0;


    // PARAMETERS FOR PARENT SELECTION
    private int tournament_size_ = 5;

    //PARAMETERS FOR MUTATION
    private double mutation_prop = 0.05;

    public player31() {
        rnd_ = new Random();
    }

    public void setSeed(long seed) {
        // Set seed of algortihm's random process
        rnd_.setSeed(seed);
    }

    public void setEvaluation(ContestEvaluation evaluation) {
        /**
         Set evaluation function used in the run.
         */
        evaluation_ = evaluation;

        // Get evaluation properties
        Properties props = evaluation.getProperties();

        // Get evaluation limit (Different for different evaluation functions)
        // Katsuura:   1.000.000
        // Schaffers:    100.000
        // BentCigar:     10.000
        // Sphere:        10.000
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));

        // Property keys depend on specific evaluation
        // E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

        // Do sth with property values, e.g. specify relevant settings of your algorithm

        // TODO: Following, we can adjust our algorithm to different evaluation functions.
        if (isMultimodal) {
            // Multimodal ->  KatsuuraEvaluation

        } else if (hasStructure) {
            // Regular ->   SchaffersEvaluation

        } else if (isSeparable) {
            // Seperable -> SphereEvaluation  (aka. Dummy-Evaluation/ Hillclimber)

        } else {
            // None -> BentCigarFunction

        }
    }

    private ArrayList<Individual> initialize_population() {
        /**
         Creates a random new population of individuals. The population is
         represented by a 2-dim array with doubles (size [population_size_][individual_size_]).
         The single values are from [-5,5].

         @return: double [population_size_][individual_size_]
         */
        ArrayList<Individual> population = new ArrayList<Individual>(population_size_);

        for (int i = 0; i < population_size_; i++) {
            double[] individual = new double[individual_size_];
            for (int j = 0; j < 10; j++) {
                double random_double = values_min_ + (values_max_ - values_min_) * rnd_.nextDouble();
                individual[j] = random_double;
            }
            //evaluate generated individual
            double fitness_i = (double) evaluation_.evaluate(individual);
            evaluations_counter_ += 1;
            // add individual to population
            population.add(new Individual(individual, fitness_i));
        }
        return population;
    }

    private void printPopulationFitness(ArrayList<Individual> population) {
      /*
      Prints the fitness values of the current population.

      @param population: ArrayList<Individual> of all individuals
                         whose fitness we want to print.
      */
        System.out.println("--------------------------------");
        System.out.println("Fitness: ");
        for (Individual individual : population) {
            double x = individual.getFitness();
            System.out.println(x);
        }
        System.out.println("--------------------------------");
    }

    private void printClusters(double[][] clusters) {
        for (int i = 0; i < clusters.length; i++) {
            double[] cluster = clusters[i];
            for (int j = 0; j < individual_size_; j++) {
                String s = Double.toString(cluster[j]);
                System.out.print(s.substring(0, 5).concat("\t"));
            }
            System.out.println("\n------------------------------------------------------------------------------\n");
        }
        System.out.println("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    }

    private void printIndividual(Individual in) {
        double[] values = in.getValues();
        double fitness = in.getFitness();
        for (int i = 0; i < values.length; i++) {
            String s = Double.toString(values[i]);
            System.out.print(s.substring(0, 5).concat("\t"));
        }
        System.out.println(fitness);
        System.out.println("\n------------------------------------------------------------------------------\n");
    }

    private ArrayList<Individual> tournamentSelection(ArrayList<Individual> pool) {
      /*
      Select two parents from pool of Individuals based on tournament selection.

      @param pool: ArrayList<Individual> of all individuals we can choose from
      @return: ArrayList<Individual> with two selected parents
      */
        ArrayList<Individual> mating_pool = new ArrayList<Individual>(tournament_size_);

        // Select tournament_size_ individuals from pool.
        while (mating_pool.size() < tournament_size_) {
            Individual i = pool.get(rnd_.nextInt(pool.size()));
            mating_pool.add(i);
        }

        // return the two fittest
        Collections.sort(mating_pool);
        return new ArrayList<Individual>(mating_pool.subList(0, 2));
    }

    private ArrayList<Individual> blendCrossover(ArrayList<Individual> parents) {
      /*
      Crossover two floatingpoint parents with Blend Crossover.

      @param parents: Should be of size two, we use the first two individuals in the list
      as parents.

      @returns: ArrayList<Individual> of size two, containing two children
      */

        double a = 0.5; // parameter, mentioned in book p.67
        // get parents value arrays
        double p0[] = parents.get(0).getValues();
        double p1[] = parents.get(1).getValues();
        // Initialize children value arrays
        double c0[] = new double[individual_size_];
        double c1[] = new double[individual_size_];

        double u = rnd_.nextDouble();
        double gamma = ((1 + (2 * a)) * u) - a;


        ArrayList<Individual> children = new ArrayList<Individual>(2);
        for (int i = 0; i < individual_size_; i++) {
            c0[i] = (1 - gamma) * p0[i] + gamma * p1[i];
            c1[i] = (1 - gamma) * p1[i] + gamma * p0[i];
        }

        children.add(new Individual(c0));
        children.add(new Individual(c1));
        return children;
    }

    private ArrayList<Individual> onePointCrossover(ArrayList<Individual> parents) {
      /*
      Crossover two floatingpoint parents with One Point Crossover.

      @param parents: Should be of size two, we use the first two individuals in the list
      as parents.

      @returns: ArrayList<Individual> of size two, containing two children
      */
        double p0[] = parents.get(0).getValues();
        double p1[] = parents.get(1).getValues();
        // Initialize children value arrays
        double c0[] = new double[individual_size_];
        double c1[] = new double[individual_size_];

        // choose random point for crossover
        int crossoverpoint = rnd_.nextInt(individual_size_);

        ArrayList<Individual> children = new ArrayList<Individual>(2);
        for (int i = 0; i < individual_size_; i++) {
            if (i < crossoverpoint) {
                c0[i] = p0[i];
                c1[i] = p1[i];
            } else {
                c0[i] = p1[i];
                c1[i] = p0[i];
            }
        }

        children.add(new Individual(c0));
        children.add(new Individual(c1));
        return children;
    }

    private void uniformMutation(ArrayList<Individual> individuals) {
      /*
      Performs uniform mutation on each of the individuals.

      */
        for (Individual i : individuals) {
            double values[] = i.getValues();
            for (int j = 0; j < values.length; j++) {
                if (rnd_.nextDouble() < mutation_prop) { // p = mutation_prop that this happens
                    double random_double = values_min_ + (values_max_ - values_min_) * rnd_.nextDouble();
                    values[j] = (double) random_double;
                }
            }
        }
    }

    private double[][] calcElitistCluster(ArrayList<Individual> elitist, double distance_threshold) {
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

        double[][] distance_matrix = new double[elitist_size_][elitist_size_];

        // fill distance matrix of elitist individuals
        for (int i = 0; i < elitist_size_; i++) {
            Individual current_elite = elitist.get(i);
            for (int j = 0; j < elitist_size_; j++) {
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
//        for (int i = 0; i < elitist_size_; i++) {
//            double[] distance_ind_i = distance_matrix[i];
//            for (int j = 0; j < elitist_size_; j++) {
//                System.out.print(distance_ind_i[j]);
//                System.out.print(" ++ ");
//            }
//            System.out.println("\n------------------------------------------------------------------------------\n");
//        }

        Map<Integer, Integer> hierarchical_clusters = new HashMap<Integer, Integer>();

        // initialize clustering => each individual is its own cluster
        for (int i = 0; i < elitist_size_; i++) {
            hierarchical_clusters.put(i, i);
        }

        double current_distance = 0.0;

        while (distance_threshold > current_distance) {
            int[] indices = findMin(distance_matrix);
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
        for (int i = 0; i < elitist_size_; i++) {
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
//                for (int j = 0; j < ind_list.size(); j++) {
//                    printIndividual(ind_list.get(j));
//                }
//                System.out.println("\n+++++++++++++++++++++++++++++++++++++++++++++++++\n");
//            }
//        }

        // Calculate centroid of each cluster
        double[][] cluster_centroids = new double[cluster_to_ind.size()][individual_size_];
        for (int i = 0; i < cluster_to_ind.size(); i++) {
            ArrayList<Individual> cluster = cluster_to_ind.get(i);
            double[] centroid = calcCentroid(cluster);
            cluster_centroids[i] = centroid;
        }

        return cluster_centroids;

    }

    /*
     * Calculates centroid of one cluster
     *
     * @param: Individuals in an ArrayList representing a cluster.
     * */
    private double[] calcCentroid(ArrayList<Individual> cluster) {
        double[] sum_centroid = new double[individual_size_];

        for (int i = 0; i < cluster.size(); i++) {
            double[] ind_vec = cluster.get(i).getValues();

            for (int j = 0; j < individual_size_; j++) {
                double dim_value = ind_vec[j];
                sum_centroid[j] += dim_value;
            }
        }

        double[] centroid = new double[individual_size_];
        for (int i = 0; i < individual_size_; i++) {
            centroid[i] = sum_centroid[i] / cluster.size();
        }

        return centroid;
    }


    private int[] findMin(double[][] distance_matrix) {
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

        for (int i = 0; i < elitist_size_; i++) {
            for (int j = 0; j < elitist_size_; j++) {
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


    private double euclideanDistance(double[] a, double[] b) {
        /*
         * Calculates the euclidean distance between two individuals
         *
         * @param a, b: the position of an individual in the search space
         *
         * @return: The euclidean distance
         * */
        double sum = 0;
        for (int i = 0; i < individual_size_; i++) {
            sum = sum + Math.pow((a[i] - b[i]), 2.0);
        }
        return Math.sqrt(sum);
    }

    private ArrayList<Individual> getElitistGroup(ArrayList<Individual> population) {
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

        for (int i = 0; i < elitist_size_; i++) {
            elitist.add(population.get(i));
        }

        return elitist;
    }


    private Individual createRandomIndividual() {
        /*
         * Creates a random individual in the search space and returns it.
         * */

        double[] individual = new double[individual_size_];
        for (int i = 0; i < individual_size_; i++) {
            double random_double = values_min_ + (values_max_ - values_min_) * rnd_.nextDouble();
            individual[i] = random_double;
        }
        //evaluate generated individual
        double fitness_i = (double) evaluation_.evaluate(individual);
        evaluations_counter_ += 1;

        return new Individual(individual, fitness_i);
    }

    private boolean isIndividualInCluster(double[][] clusters, Individual individual) {
        /*
         * Tests if an individual is in the spatial range of a cluster.
         * */

        for (int i = 0; i < clusters.length; i++) {
            double[] current_cluster = clusters[i];
            double distance = euclideanDistance(current_cluster, individual.getValues());

            if (distance < cluster_distance_thresh_) {
                return true;
            }
        }
        return false;
    }

    private ArrayList<Individual> generateRandomProletarians(ArrayList<Individual> population, double[][] cluster) {
        /*
        * Deletes the worst individuals in a population depending on the proletarian_size_ and refills the population
        * with random individuals outside the clusters determined by the individuals of the elitist group.
        * */

        Collections.sort(population);
        population = new ArrayList<Individual>(population.subList(0, (population_size_ - proletarian_size_)));

        int current_proletarians = 0;
        while (current_proletarians < proletarian_size_) {
            Individual proletarian = createRandomIndividual();
            if (!isIndividualInCluster(cluster, proletarian)) {
                population.add(proletarian);
                current_proletarians++;
            }

        }
        return population;
    }


/*
    TODO: Implement first basic approach => maybe at least non-uniform mutation for that?
    TODO: Improve parent selection? => what selection pressure do we want to have?
    TODO: Implement more sophisticated approach with uncorrelated mutation with n step size
    TODO: Extend both approaches (simple and sophisticated) with the proletarian group approach
    TODO: Implement function to save data to csv files for visualization
    TODO: Implement visualization => reason about data ...
*/

    public void run() {
        // Run your algorithm here

        // Initialize random population
        ArrayList<Individual> population = initialize_population();

        // Rank population in fitness
        Collections.sort(population);

        // print initial population for testing
//        for (int i = 0; i < population_size_; i++) {
//            printIndividual(population.get(i));
//        }
//        printPopulationFitness(population); // just for testing

        while (evaluations_counter_ < evaluations_limit_) {

            ArrayList<Individual> elitist = getElitistGroup(population);
            double[][] elitist_clusters = calcElitistCluster(elitist, cluster_distance_thresh_);
//            printClusters(elitist_clusters);
            System.out.println(elitist_clusters.length);


            // Create next generation (offspring)
            ArrayList<Individual> offspring = new ArrayList<Individual>(population_size_);
            while (offspring.size() < population_size_) {

                // PARENT SELECTION (parents has size 2)
                ArrayList<Individual> parents = tournamentSelection(population);

                // RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // (1) One Point Crossover
                // ArrayList<Individual> children = onePointCrossover(parents);

                // (2) Blend Crossover
                ArrayList<Individual> children = blendCrossover(parents);

                // MUTATE children
                // (1) Uniform mutation
                uniformMutation(children);

                // Evaluate final children's fitness
                for (Individual child : children) {
                    child.setFitness(((double) evaluation_.evaluate(child.getValues())));
                    evaluations_counter_ += 1;
                }
                // Add children to offspring
                offspring.addAll(children);
            }

            // SURVIVOR SELECTION
            // (1) Just replace whole population by offspring
            //  population = new ArrayList<Individual>(offspring);

            // (2) Rank population+offspring and select population_size_ best.
            population.addAll(offspring);
            Collections.sort(population);
            population = new ArrayList<Individual>(population.subList(0, population_size_));

// Function that introduces the random individuals => just needs to be properly embedded into the evaluations (causes Null pointer at the moment)
//            population = generateRandomProletarians(population, elitist_clusters);

        }
    }
}


class Individual implements Comparable<Individual> {
    /*
    Class that wraps a single individual and its fitness.
    Enables easy comparison among individuals.
    */
    private double[] values_;
    private double fitness_;

    public Individual(double values[]) {
        values_ = values;
    }

    public Individual(double values[], double fitness) {
        fitness_ = fitness;
        values_ = values;
    }

    public double getFitness() {
        return fitness_;
    }

    public double[] getValues() {
        return values_;
    }

    public void setFitness(double fitness) {
        fitness_ = fitness;
    }

    public void setValues(double values[]) {
        values_ = values;
    }

    @Override
    public int compareTo(Individual i) {
        return this.fitness_ < i.getFitness() ? 1 : this.fitness_ > i.getFitness() ? -1 : 0;
    }
}
