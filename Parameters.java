class Parameters {
    /*
     * Class that holds all parameters of our Algorithm.
     * Consists of only public parameters.
     */

    // BASIC PROBLEM PARAMETERS (DO NOT CHANGE)
    public static final int individual_size = 10;
    public static final double values_min = -5.0;
    public static final double values_max = 5.0;
    // number of individuals of one population
    public static int population_size = 500;


    // Parameters for the clustering
    public int elitist_size = 10 * population_size / 100;
    public int proletarian_size = 15 * population_size / 100;
    public double cluster_distance_thresh = 2.2690223010178983;


    // Parameters for the EA-components
    public int offspring_size = 650;

    // (1) Parameters for parent selection
//    public int parent_tournament_size = 47;

    // (2) Parameters for mutation
    public double uniform_mutation_prop = 0.05;
    public double non_uniform_mutation_step_size = 0.3;

    // Parameters regarding the adaptive mutation
    public double tau = 1 / Math.sqrt(2 * individual_size);
    public double tau_prime = 1 / Math.sqrt(2 * Math.sqrt(individual_size));
    public double boundary = 0.01;

    // (3) Parameters for survivor selection
//    public int survivor_tournament_size = 4;

    // (4) Parameters for hybridisation
    public int evaluations_per_proletarian = 10;
    public boolean use_hybridisation = true; //If true, hybridisation will be used in the 'sophisticated approach with own' on multimodal functions only
    public double hill_climb_step_size = 0.9601293606479896;

    // sophisticated approach Bentcigar
    public int parent_tournament_size = 39;

    //Simple approach Bentcigar
//    public int offspring_size=123;
    public int survivor_tournament_size=4;
    //    public double non_uniform_mutation_step_size=0.15143002448265783;

//Simple approach Schaffers
// offspring_percentage=130.12368486099646
// parent_tournament_size=39
// survivor_tournament_size=4
// non_uniform_mutation_step_size=0.15143002448265783

//Simple approach Katsuura
// offspring_percentage= 134.15988296420113
// parent_tournament_size=40
// survivor_tournament_size=16
// non_uniform_mutation_step_size=1.9916534678949145


//Simple approach Clustering bentcigar
//    public int offspring_percentage= 123.9789487679547
//    parent_tournament_size=10
//    survivor_tournament_size=17
//    non_uniform_mutation_step_size=0.06542202660344815
//    population_size=100
//    cluster_distance_thresh=1.3260420688384487
//    hill_climb_step_size=0.18175913741175742


}

// Soph with clustering bentcigar
// parent_tournament_size=44
// cluster_distance_thresh=1.6097353905781224
// hill_climb_step_size=0.19520320686518625


//Katsuura Simple approach without clustering
// offspring_percentage= 134.15988296420113
// parent_tournament_size=40
// survivor_tournament_size=16
// non_uniform_mutation_step_size=1.9916534678949145


// Schaffers simple approach clustering
// offspring_percentage=111.30764393411431
// parent_tournament_size=40
// survivor_tournament_size=5
// non_uniform_mutation_step_size=0.09650924903657412
// population_size=500
// cluster_distance_thresh=2.2690223010178983
// hill_climb_step_size=0.9601293606479896
