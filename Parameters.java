class Parameters{
  /*
   * Class that holds all parameters of our Algorithm.
   * Consists of only public parameters.
   */

    // BASIC PROBLEM PARAMETERS (DO NOT CHANGE)
    public static final int individual_size = 10;
    public static final double values_min = -5.0;
    public static final double values_max =  5.0;
    // number of individuals of one population
    public static final int population_size = 100;


    // Parameters for the clustering
    public int elitist_size = 5;
    public int proletarian_size = 13;
    public double cluster_distance_thresh = 1.3007402405183437;


    // Parameters for the EA-components
    public int offspring_size = 100;  //TODO: Should that be more/less ?

    // (1) Parameters for parent selection
    public int parent_tournament_size = 14;

    // (2) Parameters for mutation
    public double uniform_mutation_prop = 0.05;
    public double non_uniform_mutation_step_size = 0.3;

    // Parameters regarding the adaptive mutation
    public double tau = 1 / Math.sqrt(2*individual_size);
    public double tau_prime = 1 / Math.sqrt(2*Math.sqrt(individual_size));
    public double boundary = 0.01;
    
    // (3) Parameters for survivor selection
    public int survivor_tournament_size = 4;

    // (4) Parameters for hybridisation
    public int evaluations_per_proletarian = 81;
    public boolean use_hybridisation = true; //If true, hybridisation will be used in the 'sophisticated approach with own' on multimodal functions only
    public double hill_climb_step_size = 0.6827029110776486;

}
