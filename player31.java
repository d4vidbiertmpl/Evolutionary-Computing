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
    private int evaluations_counter_ = 0 ; // counter for performed evaluations, starts at 0

    private int population_size_ = 10; // size of one population

    // PARAMETERS FOR PARENT SELECTION
    private int tournament_size_ = 3;

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
      Prints the fitnesses of the current population.

      @param population: ArrayList<Individual> of all individuals
                         whose fitness we want to print.
      */
        System.out.println("--------------------------------");
        System.out.println("Fitness: ");
        for (Individual individual : population){
          double x = individual.getFitness();
          System.out.println(x);
        }
        System.out.println("--------------------------------");
    }

    private void printIndividual(Individual in){
      double values[] = in.getValues();
      for(int i=0; i<values.length; i++){
        String s = Double.toString(values[i]);
        System.out.print(s.substring(0,5).concat("\t"));
      }
      System.out.println("\n------------------------------------------------------------------------------\n");
    }

    private ArrayList<Individual> tournamentSelection(ArrayList<Individual> pool){
      /*
      Select two parents from pool of Individuals based on tournament selection.

      @param pool: ArrayList<Individual> of all individuals we can choose from
      @return: ArrayList<Individual> with two selected parents
      */
      ArrayList<Individual> mating_pool = new ArrayList<Individual>(tournament_size_);

      // Select tournament_size_ individuals from pool.
      while(mating_pool.size()<tournament_size_){
        Individual i = pool.get(rnd_.nextInt(pool.size()));
        mating_pool.add(i);
      }

      // return the two fittest
      Collections.sort(mating_pool);
      return new ArrayList<Individual>(mating_pool.subList(0,2));
    }

    private ArrayList<Individual> blendCrossover(ArrayList<Individual> parents){
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
      double gamma = ((1+(2*a))*u)-a;


      ArrayList<Individual> children  = new ArrayList<Individual>(2);
      for(int i=0; i<individual_size_; i++){
        c0[i] = (1-gamma)*p0[i]+gamma*p1[i];
        c1[i] = (1-gamma)*p1[i]+gamma*p0[i];
      }

      children.add(new Individual(c0));
      children.add(new Individual(c1));
      return children;
    }

    private ArrayList<Individual> onePointCrossover(ArrayList<Individual> parents){
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

      ArrayList<Individual> children  = new ArrayList<Individual>(2);
      for(int i=0; i<individual_size_; i++){
        if (i<crossoverpoint){
          c0[i] = p0[i];
          c1[i] = p1[i];
        } else{
          c0[i] = p1[i];
          c1[i] = p0[i];
        }
      }

      children.add(new Individual(c0));
      children.add(new Individual(c1));
      return children;
    }

    private void uniformMutation(ArrayList<Individual>  individuals){
      /*
      Performs uniform mutation on each of the individuals.

      */
      for(Individual i : individuals){
        double values[] = i.getValues();
        for(int j=0; j<values.length; j++){
          if (rnd_.nextDouble()<mutation_prop){ // p = mutation_prop that this happens
            double random_double = values_min_ + (values_max_ - values_min_) * rnd_.nextDouble();
            values[j] = (double) random_double;
          }
        }
      }
    }

    public void run() {
        // Run your algorithm here

        // Initialize random population
        ArrayList<Individual> population = initialize_population();

        // Rank population in fitness
        Collections.sort(population);

        //printPopulationFitness(population); // just for testing

        while (evaluations_counter_ < evaluations_limit_) {


            // Create next generation (offspring)
            ArrayList<Individual> offspring = new ArrayList<Individual>(population_size_);
            while (offspring.size() < population_size_){

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
              for (Individual child: children){
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
          population =  new ArrayList<Individual>(population.subList(0,population_size_));


        }
      }
}


class Individual implements Comparable<Individual> {
  /*
  Class that wraps a single individual and its fitness.
  Enables easy comparison among individuals.
  */
    private double values_[];
    private double fitness_;

    public Individual (double values[])
    {
        values_ = values;
    }

    public Individual (double values[], double fitness)
    {
        fitness_  = fitness;
        values_ = values;
    }

    public double getFitness()
    {
        return fitness_;
    }

    public double[] getValues()
    {
        return values_;
    }

    public void setFitness(double fitness){
        fitness_  = fitness;
    }

    public void setValues(double values[]){
        values_ = values;
    }

    @Override
    public int compareTo(Individual i) {
        return this.fitness_ < i.getFitness() ? 1 : this.fitness_ > i.getFitness() ? -1 : 0;
    }
}
