import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import javax.naming.directory.InitialDirContext;
import java.util.*;
import java.util.Comparator;

public class player31 implements ContestSubmission {

    // object for creating random numbers
    Random rnd_;

    // object for evaluating created individuals
    ContestEvaluation evaluation_;

    // parameters
    Parameters parameters;

    // EA utils (= EA-functions applied on individuals)
    EA_Utils ea_utils;

    // Clustering Utils (= Clustering-functions applied on individuals)
    Clustering_Utils clustering_utils;


    private int evaluations_limit_;   // number of allowed evaluations
    private int evaluations_counter_ = 0; // counter for performed evaluation


    public player31() {
        rnd_ = new Random();
        parameters = new Parameters();
        ea_utils = new EA_Utils();
        clustering_utils = new Clustering_Utils();
    }

    public void setSeed(long seed) {
        rnd_.setSeed(seed);
    }

    public void setEvaluation(ContestEvaluation evaluation) {
        //Set evaluation function used in the run.
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

/*
    TODO: Implement first basic approach => maybe at least non-uniform mutation for that?

    TODO: Improve parent selection? => what selection pressure do we want to have?

    TODO: Implement more sophisticated approach with uncorrelated mutation with n step size

    TODO: Extend both approaches (simple and sophisticated) with the proletarian group approach

    TODO: Implement function to save data to csv files for visualization

    TODO: Implement visualization => reason about data ...
*/

    public void run() {

        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population();
        evaluateIndividuals(population);

        // Rank population in fitness
        Collections.sort(population);

        while (evaluations_counter_ < evaluations_limit_) {  // until no exception

            // Get top Individuals
            // ArrayList<Individual> elitist = getElitistGroup(population

            // Calculate Clusters
            // double[][] elitist_clusters = calcElitistCluster(elitist, parameters.cluster_distance_thresh);
            // printClusters(elitist_clusters); // for test purpose
            // System.out.println(elitist_clusters.length);


            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.offspring_size) {

                // (1) PARENT SELECTION (parents has size 2)
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // (2.1) One Point Crossover
                // ArrayList<Individual> children = onePointCrossover(parents);

                // (2.2) Blend Crossover
                ArrayList<Individual> children = ea_utils.blendCrossover(parents);

                // (3) MUTATE children
                // (3.1) Uniform mutation
                for (Individual child: children){
		    ea_utils.uniformMutation(child);
		}

		// (3.2) Non-uniform mutation
                //for (Individual child: children){
		//  ea_utils.nonUniformMutation(child);
		//}

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) SURVIVOR SELECTION
            // (4.1) Just replace whole population by offspring
            // population = new ArrayList<Individual>(offspring);

            // (4.2) Rank population+offspring and select parameters.population_size best.
            // population.addAll(offspring);
            // Collections.sort(population);
            // population = new ArrayList<Individual>(population.subList(0, parameters.population_size));

            // (4.3) tournamentSelection
            population.addAll(offspring);
            ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
            while (new_population.size() < parameters.population_size) {
              new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
            }

        population = new ArrayList<>(new_population); // copy new_population to population

          // Function that introduces the random individuals => just needs to be properly embedded into the evaluations (causes Null pointer at the moment)
          // population = generateRandomProletarians(population, elitist_clusters);

        }
    }

    private void evaluateIndividuals(ArrayList<Individual> individuals){
      /*
       * Evaluate + assign the fitness of every individual in individuals that
       * has not yet been assigned a fitness.
       */
       for(Individual individual : individuals){
         if (!individual.isEvaluated()){ // if individual was not evaluated yet
           double fitness = (double) evaluation_.evaluate(individual.getValues());
           individual.setFitness(fitness);
           evaluations_counter_ +=1;
         }
       }
    }

    private double measureDiversity(ArrayList<Individual> individuals){
      /*
       * Returns a measure of how diverse the individuals are.
       * This measure is the average euclidian distance to the
       * average individual.
       */
       double diff = 0.0;
       double average[] =  clustering_utils.calcCentroid(individuals);

       for (Individual individual : individuals){
         double values[] = individual.getValues();
         for(int i=0; i<parameters.individual_size; i++){
           diff += clustering_utils.euclideanDistance(average, values);
         }
       }
       return diff/individuals.size();
    }


    // --------------------------------------------------------------------------
    // Following functions are just needed for testing purposes.
    // -------------------------------------------------------------------------
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
            for (int j = 0; j < parameters.individual_size; j++) {
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
}
