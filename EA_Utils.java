
import java.util.*;
import java.util.Comparator;

class EA_Utils{
  /*
   *   Class that holds all EA functions that can be applied on Individuals, such as
   *
   *   - Variation :  Recombination, Mutation
   *   - Selection :  Parent Selection, Survivor Selection
   *
   * */
  private Parameters parameters;

  private Random rnd_;

  public EA_Utils() {
      parameters = new Parameters();
      rnd_ = new Random();
      rnd_.setSeed(16438320L);
  }

  // --------------------------------------------------------------------------
  // (1) Functions for Intitalization
  // --------------------------------------------------------------------------

  public Individual createRandomIndividual() {
      /*
       * Creates a random individual in the search space and returns it.
       * */
      double[] individual = new double[parameters.individual_size];
      for (int i = 0; i < parameters.individual_size; i++) {
          double random_double = parameters.values_min + (parameters.values_max - parameters.values_min) * rnd_.nextDouble();
          individual[i] = random_double;
      }
      return new Individual(individual);
  }

  public ArrayList<Individual> initialize_population() {
      /*
       * Creates a random new population of individuals.
       * */
      ArrayList<Individual> population = new ArrayList<Individual>(parameters.population_size);
      for (int i = 0; i < parameters.population_size; i++) {
          Individual ramdom_individual = createRandomIndividual();
          population.add(ramdom_individual);
      }
      return population;
  }
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------



  // --------------------------------------------------------------------------
  // (2) Functions for Selection
  // --------------------------------------------------------------------------
  public ArrayList<Individual> tournamentSelection(ArrayList<Individual> pool, int tourn_size) {
    /*
     * Select two parents from pool of Individuals based on tournament selection.
     *
     * @param pool: ArrayList<Individual> of all individuals we can choose from
     * @param tourn_size: Size of tournament
     *
     *@return: ArrayList<Individual> with two selected parents
     * */
     ArrayList<Individual> mating_pool = new ArrayList<Individual>(tourn_size);

     //System.out.println(pool.size());
     // Select tourn_size individuals from pool.
     while (mating_pool.size() < tourn_size) {
          int index = rnd_.nextInt(pool.size());
          Individual i = pool.get(index);
          //System.out.println(index);
          mating_pool.add(i);
     }

     //System.out.println("---------");

     // return the two fittest
     Collections.sort(mating_pool);
     return new ArrayList<Individual>(mating_pool.subList(0, 2));
  }
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------



  // --------------------------------------------------------------------------
  // (3) Functions for Recombination
  // -------------------------------------------------------------------------
  public ArrayList<Individual> blendCrossover(ArrayList<Individual> parents) {
    /*
     * Crossover two floatingpoint parents with Blend Crossover.
     *
     * @param parents: Should be of size two, we use the first two individuals in the list
     * as parents.
     *
     * @returns: ArrayList<Individual> of size two, containing two children
     * */

      double a = 0.5; // parameter, mentioned in book p.67

      // get parents value arrays
      double p0[] = parents.get(0).getValues();
      double p1[] = parents.get(1).getValues();
      // Initialize children value arrays
      double c0[] = new double[parameters.individual_size];
      double c1[] = new double[parameters.individual_size];

      double u = rnd_.nextDouble();
      double gamma = ((1 + (2 * a)) * u) - a;

      ArrayList<Individual> children = new ArrayList<Individual>(2);
      for (int i = 0; i < parameters.individual_size; i++) {
          c0[i] = (1 - gamma) * p0[i] + gamma * p1[i];
          c1[i] = (1 - gamma) * p1[i] + gamma * p0[i];
      }
      children.add(new Individual(c0));
      children.add(new Individual(c1));
      return children;
  }

  public ArrayList<Individual> onePointCrossover(ArrayList<Individual> parents) {
    /*
     * Crossover two floatingpoint parents with One Point Crossover.
     *
     *   @param parents: Should be of size two, we use the first two individuals in the list
     *   as parents.
     *
     *   @returns: ArrayList<Individual> of size two, containing two children
     * */
      double p0[] = parents.get(0).getValues();
      double p1[] = parents.get(1).getValues();
      // Initialize children value arrays
      double c0[] = new double[parameters.individual_size];
      double c1[] = new double[parameters.individual_size];

      // choose random point for crossover
      int crossoverpoint = rnd_.nextInt(parameters.individual_size);

      ArrayList<Individual> children = new ArrayList<Individual>(2);
      for (int i = 0; i < parameters.individual_size; i++) {
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
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------



  // --------------------------------------------------------------------------
  // (4) Functions for Mutation
  // -------------------------------------------------------------------------
  public void uniformMutation(Individual individual) {
    /*
     * Performs uniform mutation on an individual.
     * */
    double values[] = individual.getValues();
    for (int j = 0; j < values.length; j++) {
        // p = parameters.uniform_mutation_prop that this happens
        if (rnd_.nextDouble() < parameters.uniform_mutation_prop) {
            double random_double = parameters.values_min + (parameters.values_max - parameters.values_min) * rnd_.nextDouble();
            values[j] = (double) random_double;
        }
    }
  }

    public void nonUniformMutation(Individual individual) {
	double values[] = individual.getValues();
	for (int j = 0; j < values.length; j++) {
	    double random_gauss = rnd_.nextGaussian() * parameters.non_uniform_mutation_step_size + values[j];
	    if (random_gauss < parameters.values_min) {
		random_gauss = parameters.values_min;
	    } else if (random_gauss > parameters.values_max) {
		random_gauss = parameters.values_max;
	    }

	    values[j] = random_gauss;
	}
    }
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
}
