import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import javax.naming.directory.InitialDirContext;
import java.util.*;
import java.util.Comparator;

public class player31 implements ContestSubmission {
    Random rnd_;
    ContestEvaluation evaluation_;
    private int evaluations_limit_;

    public player31() {
        rnd_ = new Random();
    }

    public void setSeed(long seed) {
        // Set seed of algortihms random process
        rnd_.setSeed(seed);
    }

    public void setEvaluation(ContestEvaluation evaluation) {
        // Set evaluation problem used in the run
        evaluation_ = evaluation;

        // Get evaluation properties
        Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
        // Property keys depend on specific evaluation
        // E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

        // Do sth with property values, e.g. specify relevant settings of your algorithm
        if (isMultimodal) {
            // Do sth
        } else {
            // Do sth else
        }
    }

    // Initialization => Population => Parent selection => Parents
    // => Recombination (Crossover) => Mutation => Offspring
    // => Survivor Selection => Population

    // Pseudo-code
    // Initialize Population with random candidate solutions
    // Evaluate each candidate
    // REPEAT UNTIL T_COND is met
    //  SELECT Parents
    //  RECOMBINE Parents
    //  MUTATE Offspring
    // EVALUATE new candidates
    // SELECT individuals for next generation


    private double[][] initialize_population(int population_size) {

        Random r = new Random();
        double r_min = -5.0;
        double r_max = 5.0;

        double[][] population = new double[population_size][10];

        for (int i = 0; i < population_size; i++) {
            double[] individual = new double[10];
            for (int j = 0; j < 10; j++) {
                double random_double = r_min + (r_max - r_min) * r.nextDouble();
                individual[j] = random_double;
            }
            population[i] = individual;

        }
        return population;
    }

    public double[][] evaluate_population(double[][] population, int population_size, ContestEvaluation evaluation_) {

        double[][] fitness = new double[population_size][2];

        for (int i = 0; i < population_size; i++) {
            double[] current_individual = population[i];
            double current_fitness = (double) evaluation_.evaluate(current_individual);
            fitness[i][0] = current_fitness;
            fitness[i][1] = i;
        }
        return fitness;
    }


    private double getmaxValue(double array[]) {
        double max = Arrays.stream(array).max().getAsDouble();
        return max;
    }


    private void displayArray(double[][] fitness) {
        System.out.println("--------------------------------");
        System.out.println("Fitness\t\t Index");
        for (int i = 0; i < fitness.length; i++) {
            double[] itemRecord = fitness[i];
            System.out.println(itemRecord[0]);
            System.out.println(itemRecord[1]);
        }
        System.out.println("--------------------------------");
    }

    //        for (int i = 0; i < 10; i++) {
//            System.out.print(Double.toString(population[1][i]).concat(" "));
//        }

    public void run() {
        // Run your algorithm here

        MatrixComp matrixComp = new MatrixComp();

        int evals = 0;

        int population_size = 10;

        // Initialize random population
        double[][] population = initialize_population(population_size);

        // evaluate population
        double[][] pop_fittness = evaluate_population(population, population_size, evaluation_);
        evals += population_size;

        displayArray(pop_fittness);

        Arrays.sort(pop_fittness, matrixComp);

        displayArray(pop_fittness);


        // calculate fitness
        while (evals < evaluations_limit_) {

//            System.out.print("HELLO");

            evals += 1000;

            // find a way to only evaluate the children every generation
            // Evaluate whole population before loop
            // Select parents
            // Generate children with crossover and mutate
            // Evaluate children
            // concat children to parents and rank the
            // select survivors with specific method, alter population and fitness array respectively
            // Now new population with respective fittness values => only evaluate # of new children lambda


            // Select parents

            // Apply crossover / mutation operators
            //double child[] = {0.2, 0.8, 0.7, 0.0, 0.2, 0.6, 0.3, 0.3, 0.4, 0.6};


            // Check fitness of unknown function
            //Double fitness = (double) evaluation_.evaluate(child);
            //evals++;
            // Select survivors
        }

    }
}


class MatrixComp implements Comparator<double[]> {

    public int compare(double[] o1, double[] o2) {
        //get the item ids which are at index 0 of the array
        double itemIdOne = o1[0];
        double itemIdTwo = o2[0];
        // sort on item id
        return Double.compare(itemIdOne, itemIdTwo);
    }
}

