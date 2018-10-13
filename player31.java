import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import javax.naming.directory.InitialDirContext;
import java.lang.reflect.Array;
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

    private int generation_;


    public player31() {
        rnd_ = new Random();
        parameters = new Parameters();
        ea_utils = new EA_Utils(parameters);
        clustering_utils = new Clustering_Utils(parameters);

        setParams();
    }

    public void setSeed(long seed) {
        rnd_.setSeed(seed);
        ea_utils.setSeed(rnd_.nextLong());
        clustering_utils.setSeed(rnd_.nextLong());
    }

    public void setParams() {

        // parse parameters from command line
        if (System.getProperty("non_uniform_mutation_step_size") != null) {
            parameters.non_uniform_mutation_step_size = Double.parseDouble(System.getProperty("non_uniform_mutation_step_size"));
        }
        if (System.getProperty("offspring_size") != null) {
          double d = Double.parseDouble(System.getProperty("offspring_size"));
          int i = (int) d;
          parameters.offspring_size = i;
        }
        if (System.getProperty("parent_tournament_size") != null) {
          double d = Double.parseDouble(System.getProperty("parent_tournament_size"));
          int i = (int) d;
          parameters.parent_tournament_size = i;
        }
        if (System.getProperty("survivor_tournament_size") != null) {
          double d = Double.parseDouble(System.getProperty("survivor_tournament_size"));
          int i = (int) d;
          parameters.survivor_tournament_size = i;
        }

        // Parameters for the clustering
        if (System.getProperty("elitist_size") != null) {
          double d = Double.parseDouble(System.getProperty("elitist_size"));
          int i = (int) d;
          parameters.elitist_size = i;
        }
        if (System.getProperty("proletarian_size") != null) {
          double d = Double.parseDouble(System.getProperty("proletarian_size"));
          int i = (int) d;
          parameters.proletarian_size = i;
        }
        if (System.getProperty("cluster_distance_thresh") != null) {
          parameters.cluster_distance_thresh = Double.parseDouble(System.getProperty("cluster_distance_thresh"));
        }
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

        // Which experiment do we want to execute?
        // simple_approach, simple_approach with own, sophisticated approach, sophisticated approach with own

//        simple_approach();
//        simple_approach_with_own();
       sophisticated_approach();
        //sophisticated_approach_with_own();

    }

    private void evaluateIndividuals(ArrayList<Individual> individuals) {
        /*
         * Evaluate + assign the fitness of every individual in individuals that
         * has not yet been assigned a fitness.
         */
        for (Individual individual : individuals) {
            if (!individual.isEvaluated()) { // if individual was not evaluated yet
                double fitness = (double) evaluation_.evaluate(individual.getValues());
                individual.setFitness(fitness);
                evaluations_counter_ += 1;
            }
        }
    }

    private double measureDiversity(ArrayList<Individual> individuals) {
        /*
         * Returns a measure of how diverse the individuals are.
         * This measure is the average euclidian distance to the
         * average individual.
         */
        double diff = 0.0;
        double average[] = clustering_utils.calcCentroid(individuals);

        for (Individual individual : individuals) {
            double values[] = individual.getValues();
            for (int i = 0; i < parameters.individual_size; i++) {
                diff += clustering_utils.euclideanDistance(average, values);
            }
        }
        return diff / individuals.size();
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

    private void printPopulationCSV(ArrayList<Individual> population) {
        String csv_population = "";

        for (int i = 0; i < parameters.population_size; i++) {
            double[] v = population.get(i).getValues();
            String csv_line = "";
            for (int j = 0; j < parameters.individual_size - 1; j++) {
                String value = Double.toString(v[j]);
                String value_comma = value.concat(",");
                csv_line = csv_line.concat(value_comma);
            }
            csv_line = csv_line.concat(Double.toString(v[9]));
            csv_line = csv_line.concat("\n");

            csv_population = csv_population.concat(csv_line);
        }
        csv_population = csv_population.concat("+");
        System.out.println(csv_population);
    }

    
    // --------------------------------------------------------------------------
    // Functions for Hybridisation
    // -------------------------------------------------------------------------

    public Individual hillClimb(Individual proletarian, int evaluations_left) {

	if (evaluations_left == 0) {
	    return proletarian;
	}

	double original_values[] = proletarian.getValues();
	double values[] = original_values.clone();
        for (int i = 0; i < original_values.length; i++) {
            double random_gauss = rnd_.nextGaussian() * parameters.non_uniform_mutation_step_size + original_values[i];
            if (random_gauss < parameters.values_min) {
                random_gauss = parameters.values_min;
            } else if (random_gauss > parameters.values_max) {
                random_gauss = parameters.values_max;
            }

            values[i] = random_gauss;
	}       	
	
	double possible_impr_fitness = (double) evaluation_.evaluate(values);
        evaluations_counter_ += 1;
	if (possible_impr_fitness > proletarian.getFitness()) {
	    Individual possible_improvement = new Individual(values, proletarian.isAdaptive_());
	    possible_improvement.setFitness(possible_impr_fitness);
	    return hillClimb(possible_improvement, evaluations_left -1);
	} else {
	    return hillClimb(proletarian, evaluations_left - 1);
	}
    }

    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------

    private void simple_approach() {

        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(false);
        evaluateIndividuals(population);

        // Rank population in fitness
        Collections.sort(population);

        while (evaluations_counter_ < evaluations_limit_) {  // until no exception

            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.offspring_size) {

                // (1) PARENT SELECTION (parents has size 2)
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // (2.1) One Point Crossover
                // ArrayList<Individual> children = onePointCrossover(parents);

                // (2.2) Blend Crossover
//                 ArrayList<Individual> children = ea_utils.blendCrossover(parents);

                // (2.3) Whole Arithmetic Recombination
                ArrayList<Individual> children = ea_utils.wholeArithmeticRecombination(parents);

                // (3) MUTATE children
                // (3.1) Uniform mutation
                for (Individual child : children) {
//                     ea_utils.uniformMutation(child);
                }

                // (3.2) Non-uniform mutation
                for (Individual child : children) {
                    ea_utils.nonUniformMutation(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) SURVIVOR SELECTION
            // (4.1) Just replace whole population by offspring
//             population = new ArrayList<Individual>(offspring);

            // (4.2) Rank population+offspring and select parameters.population_size best.
//             population.addAll(offspring);
//             Collections.sort(population);
//             population = new ArrayList<Individual>(population.subList(0, parameters.population_size));

            // (4.3) tournamentSelection
            population.addAll(offspring);
            ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
            while (new_population.size() < parameters.population_size) {
                new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
            }
            population = new ArrayList<Individual>(new_population); // copy new_population to population

        }


    }

    private void simple_approach_with_own() {

        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(false);
        evaluateIndividuals(population);

        // Rank population in fitness
        // Collections.sort(population);

        while (evaluations_counter_ < evaluations_limit_) {  // until no exception

            // Test for diversity
            //System.out.println(measureDiversity(population));

            // Get top Individuals
            ArrayList<Individual> elitist = clustering_utils.getElitistGroup(population);

            // Calculate Clusters
            double[][] elitist_clusters = clustering_utils.calcElitistCluster(elitist, parameters.cluster_distance_thresh);
            // printClusters(elitist_clusters); // for test purpose
            System.out.println(elitist_clusters.length);

            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.offspring_size) {

                // (1) PARENT SELECTION (parents has size 2)
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // (2.1) One Point Crossover
                // ArrayList<Individual> children = onePointCrossover(parents);

                // (2.2) Blend Crossover
//                 ArrayList<Individual> children = ea_utils.blendCrossover(parents);

                // (2.3) Whole Arithmetic Recombination
                ArrayList<Individual> children = ea_utils.wholeArithmeticRecombination(parents);

                // (3) MUTATE children

                for (Individual child : children) {

                    // (3.1) Uniform mutation
                    //ea_utils.uniformMutation(child);

                    // (3.2) Non-Uniform Mutation
                    ea_utils.nonUniformMutation(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) SURVIVOR SELECTION
            // (4.1) Just replace whole population by offspring
//             population = new ArrayList<Individual>(offspring);

            // (4.2) Rank population+offspring and select parameters.population_size best.
//             population.addAll(offspring);
//             Collections.sort(population);
//             population = new ArrayList<Individual>(population.subList(0, parameters.population_size));

            // (4.3) tournamentSelection
            population.addAll(offspring);

            ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
            while (new_population.size() < parameters.population_size) {
                new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
            }
            population = new ArrayList<Individual>(new_population); // copy new_population to population

            ArrayList<Individual> proletarians = clustering_utils.generateRandomProletarians(elitist_clusters, false);
            Collections.sort(population);
            population = new ArrayList<Individual>(population.subList(0, (parameters.population_size - parameters.proletarian_size)));

            for (int i = 0; i < parameters.proletarian_size; i++) {
                Individual current_proletarian = proletarians.get(i);
                double prol_fitness = (double) evaluation_.evaluate(current_proletarian.getValues());
                current_proletarian.setFitness(prol_fitness);
                evaluations_counter_ += 1;
                population.add(current_proletarian);
            }

        }
    }

    private void sophisticated_approach() {

        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(true);
        evaluateIndividuals(population);

        // Rank population in fitness
        Collections.sort(population);

        while (evaluations_counter_ < evaluations_limit_) {  // until no exception


            printPopulationCSV(population);

            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.offspring_size) {

                // (1) PARENT SELECTION (parents has size 2)
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // (2.2) Blend Crossover
                ArrayList<Individual> children = ea_utils.blendCrossover(parents);
//                ArrayList<Individual> children = ea_utils.wholeArithmeticRecombination(parents);


                // (3) MUTATE children
                // (3.1) Uniform mutation
                for (Individual child : children) {
                    ea_utils.adaptiveMutationNStSz(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) SURVIVOR SELECTION
            // (4.1) Just replace whole population by offspring
            Collections.sort(offspring);
            population = new ArrayList<Individual>(offspring.subList(0, parameters.population_size));

            // (4.2) Rank population+offspring and select parameters.population_size best.
//            population.addAll(offspring);
//            Collections.sort(population);
//            population = new ArrayList<Individual>(population.subList(0, parameters.population_size));


            // (4.3) tournamentSelection
//            population.addAll(offspring);
//            ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
//            while (new_population.size() < parameters.population_size) {
//                new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
//            }
//            population = new ArrayList<Individual>(new_population); // copy new_population to population

        }
    }

    private void sophisticated_approach_with_own() {

        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(true);
        evaluateIndividuals(population);

        // Rank population in fitness
        Collections.sort(population);

        while (evaluations_counter_ < evaluations_limit_) {  // until no exception
            generation_++;

            printPopulationCSV(population);

            // Get top Individuals
            ArrayList<Individual> elitist = clustering_utils.getElitistGroup(population);

            // Calculate Clusters
            double[][] elitist_clusters = clustering_utils.calcElitistCluster(elitist, parameters.cluster_distance_thresh);
            // printClusters(elitist_clusters); // for test purpose
//            System.out.println(elitist_clusters.length);


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
                for (Individual child : children) {
                    ea_utils.adaptiveMutationNStSz(child);
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
            Collections.sort(offspring);
            population = new ArrayList<Individual>(offspring.subList(0, parameters.population_size));

            // (4.2) Rank population+offspring and select parameters.population_size best.
//             population.addAll(offspring);
//             Collections.sort(population);
//             population = new ArrayList<Individual>(population.subList(0, parameters.population_size));

            ArrayList<Individual> proletarians = clustering_utils.generateRandomProletarians(elitist_clusters, true);
            Collections.sort(population);
            population = new ArrayList<Individual>(population.subList(0, (parameters.population_size - parameters.proletarian_size)));

            for (int i = 0; i < parameters.proletarian_size; i++) {
                Individual current_proletarian = proletarians.get(i);
                double prol_fitness = (double) evaluation_.evaluate(current_proletarian.getValues());
                current_proletarian.setFitness(prol_fitness);
                evaluations_counter_ += 1;

		Properties props = evaluation_.getProperties();
		if (parameters.use_hybridisation && Boolean.parseBoolean(props.getProperty("Multimodal"))
		    && ((evaluations_limit_ - evaluations_counter_) > parameters.evaluations_per_proletarian)) {
		    current_proletarian = hillClimb(current_proletarian, parameters.evaluations_per_proletarian);
		}

                population.add(current_proletarian);
            }

            // Function that introduces the random individuals => just needs to be properly embedded into the evaluations (causes Null pointer at the moment)
            // population = generateRandomProletarians(population, elitist_clusters);

        }
    }
}
