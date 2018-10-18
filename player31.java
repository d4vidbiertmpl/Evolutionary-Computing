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
        if (System.getProperty("population_size") != null) {
            double d = Double.parseDouble(System.getProperty("population_size"));
            int i = (int) d;
            parameters.population_size = i;

            parameters.elitist_size = (int) 10*i/100;
            parameters.proletarian_size = (int) 15*i/100;
        }
        // Parameters for the clustering
        if (System.getProperty("cluster_distance_thresh") != null) {
            parameters.cluster_distance_thresh = Double.parseDouble(System.getProperty("cluster_distance_thresh"));
        }
        if (System.getProperty("hill_climb_step_size") != null) {
            parameters.hill_climb_step_size = Double.parseDouble(System.getProperty("hill_climb_step_size"));
        }
    }

    public void setEvaluation(ContestEvaluation evaluation) {
        //Set evaluation function used in the run.
        evaluation_ = evaluation;

        // Get evaluation properties
        Properties props = evaluation.getProperties();

        // Get evaluation limit (Different for different evaluation functions)
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));

        // Property keys depend on specific evaluation
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

        // TODO: Following, we can adjust our algorithm to different evaluation functions.
        if (isMultimodal) {

        } else if (hasStructure) {
            // Regular ->   SchaffersEvaluation

        } else if (isSeparable) {
            // Seperable -> SphereEvaluation  (aka. Dummy-Evaluation/ Hillclimber)

        } else {
            // None -> BentCigarFunction
        }
    }

    public void run() {
        // Which experiment do we want to execute?

        // simple_approach();
        // simple_approach_with_own();
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
                if (evaluations_counter_ < evaluations_limit_) {
                    double fitness = (double) evaluation_.evaluate(individual.getValues());
                    individual.setFitness(fitness);
                    evaluations_counter_ += 1;
                } else {
                    individual.setFitness(0.0);
                }
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
        double res = diff / individuals.size();
        System.out.println("Diversity: ".concat(Double.toString(res)));
        return res;
    }


    // --------------------------------------------------------------------------
    // Following functions are just needed for testing purposes.
    // -------------------------------------------------------------------------
    private void printMaxFitness(ArrayList<Individual> population){
      double max_f = 0.0;
      double av_f = 0.0;
      for (Individual i: population){
          double f_i = i.getFitness();
          if (f_i > max_f){
            max_f = f_i;
          }
          av_f += f_i;
      }
      av_f = av_f/population.size();
      System.out.println("MaxFitness: ".concat(Double.toString(max_f)));
      System.out.println("AvFitness: ".concat(Double.toString(av_f)));
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
    // --------------------------------------------------------------------------



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
            double random_gauss = rnd_.nextGaussian() * parameters.hill_climb_step_size + original_values[i];
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
            return hillClimb(possible_improvement, evaluations_left - 1);
        } else {
            return hillClimb(proletarian, evaluations_left - 1);
        }
    }

    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------



    // --------------------------------------------------------------------------
    // Approaches
    // --------------------------------------------------------------------------

    private void simple_approach() {
        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(false);
        evaluateIndividuals(population);

        while (evaluations_counter_ < evaluations_limit_) {

            // Create offspring
            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);

            while (offspring.size() < parameters.offspring_size) {

                // (1) Parent Selection : Tournament Selection
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) Recombine Parents to receive 2 children (fitness not evaluated)
                //  Whole Arithmetic Recombination
                ArrayList<Individual> children = ea_utils.wholeArithmeticRecombination(parents);

                // (3) Mutate children
                // Non-uniform mutation
                for (Individual child : children) {
                    ea_utils.nonUniformMutation(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) Survivor Selection: TournamentSelection
            population.addAll(offspring);
            ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
            while (new_population.size() < parameters.population_size) {
                new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
            }
            population = new ArrayList<Individual>(new_population.subList(0, parameters.population_size));
        }
    }


    private void simple_approach_with_own() {

      // Initialize random population
      ArrayList<Individual> population = ea_utils.initialize_population(false);
      evaluateIndividuals(population);

      while (evaluations_counter_ < evaluations_limit_) {

          // Create offspring
          ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
          while (offspring.size() < parameters.offspring_size) {

              // (1) Parent Selection : Tournament Selection
              ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

              // (2) Recombine Parents to receive 2 children (fitness not evaluated)
              //  Whole Arithmetic Recombination
              ArrayList<Individual> children = ea_utils.wholeArithmeticRecombination(parents);

              // (3) Mutate children
              // Non-uniform mutation
              for (Individual child : children) {
                  ea_utils.nonUniformMutation(child);
              }

              // Evaluate final children's fitness
              evaluateIndividuals(children);

              // Add children to offspring
              offspring.addAll(children);
          }

          // (4) Survivor Selection: TournamentSelection
          population.addAll(offspring);
          ArrayList<Individual> new_population = new ArrayList<Individual>(parameters.population_size);
          while (new_population.size() < parameters.population_size) {
              new_population.addAll(ea_utils.tournamentSelection(population, parameters.survivor_tournament_size));
          }
          population = new ArrayList<Individual>(new_population.subList(0, parameters.population_size));


          // OWN APPROACH
          // Get Top Individuals
          ArrayList<Individual> elitist = clustering_utils.getElitistGroup(population);

          // Delete worst individuals
          Collections.sort(population);
          population = new ArrayList<Individual>(offspring.subList(0, (parameters.population_size - parameters.proletarian_size)));

          // Calculate Clusters + Generate Random Proletarians
          double[][] elitist_clusters = clustering_utils.calcElitistCluster(elitist, parameters.cluster_distance_thresh);
          ArrayList<Individual> proletarians = clustering_utils.generateRandomProletarians(elitist_clusters, false);
          evaluateIndividuals(proletarians);

          //Hillclimber
          // Only if function is multimodal
          for (Individual current_proletarian: proletarians){
             Properties props = evaluation_.getProperties();
              if (parameters.use_hybridisation && Boolean.parseBoolean(props.getProperty("Multimodal"))
                      && ((evaluations_limit_ - evaluations_counter_) > parameters.evaluations_per_proletarian)) {
                  current_proletarian = hillClimb(current_proletarian, parameters.evaluations_per_proletarian);
              }
          }

          // Add proletarians to population
          population.addAll(new ArrayList<Individual>(proletarians));
        }
    }

    private void sophisticated_approach() {
        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(true);
        evaluateIndividuals(population);

        while (evaluations_counter_ < evaluations_limit_) {

            //printPopulationCSV(population);
            //printMaxFitness(population);

            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.population_size) {

                // (1) Parent Selection : Tournament Selection
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) Recombine 2 Parents to receive 2 children
                // Blend Crossover
                ArrayList<Individual> children = ea_utils.blendCrossover(parents);

                // (3) Mutate Children with adaptive mutation step size
                for (Individual child : children) {
                    ea_utils.adaptiveMutationNStSz(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);
                // Add children to offspring
                offspring.addAll(children);
            }

            // (4) Instead of Survivor Selection
            // Just replace whole population by offspring
            population = new ArrayList<Individual>(offspring.subList(0, parameters.population_size));
        }
    }


    private void sophisticated_approach_with_own() {
        // Initialize random population
        ArrayList<Individual> population = ea_utils.initialize_population(true);
        evaluateIndividuals(population);

        while (evaluations_counter_ < evaluations_limit_) {

            // Create Offspring from population
            ArrayList<Individual> offspring = new ArrayList<Individual>(parameters.population_size);
            while (offspring.size() < parameters.population_size) {

                // (1) PARENT SELECTION (parents has size
                ArrayList<Individual> parents = ea_utils.tournamentSelection(population, parameters.parent_tournament_size);

                // (2) RECOMBINE Parents to receive 2 children (fitness not evaluated)
                // Blend Crossover
                ArrayList<Individual> children = ea_utils.blendCrossover(parents);

                // (3) MUTATE children
                for (Individual child : children) {
                    ea_utils.adaptiveMutationNStSz(child);
                }

                // Evaluate final children's fitness
                evaluateIndividuals(children);

                // Add children to offspring
                offspring.addAll(children);
            }
            // (4) Instead of Survivor Selection
            // Just replace whole population by offspring
            population = new ArrayList<Individual>(offspring.subList(0, parameters.population_size));


            // OWN APPROACH
            // Get Top Individuals
            ArrayList<Individual> elitist = clustering_utils.getElitistGroup(population);

            // Delete worst individuals
            Collections.sort(population);
            population = new ArrayList<Individual>(offspring.subList(0, (parameters.population_size - parameters.proletarian_size)));

            // Calculate Clusters + Generate Random Proletarians
            double[][] elitist_clusters = clustering_utils.calcElitistCluster(elitist, parameters.cluster_distance_thresh);
            ArrayList<Individual> proletarians = clustering_utils.generateRandomProletarians(elitist_clusters, true);

            // mutate + evaluate proletarians to also test their mutation step size
            for (Individual current_proletarian: proletarians){
                ea_utils.adaptiveMutationNStSz(current_proletarian);
            }
            evaluateIndividuals(proletarians);

            //Hillclimber
            // Only if function is multimodal
            for (Individual current_proletarian: proletarians){
              Properties props = evaluation_.getProperties();
               if (parameters.use_hybridisation && Boolean.parseBoolean(props.getProperty("Multimodal"))
                       && ((evaluations_limit_ - evaluations_counter_) > parameters.evaluations_per_proletarian)) {
                   current_proletarian = hillClimb(current_proletarian, parameters.evaluations_per_proletarian);
               }
            }
              // Add proletarians to population
            population.addAll(new ArrayList<Individual>(proletarians));
        }
    }
}
