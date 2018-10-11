import java.util.Arrays;

class Individual implements Comparable<Individual> {
    /*
    Class that wraps a single individual and its fitness.
    Enables easy comparison among individuals.
    */

    private double[] values_;
    private double fitness_;
    private boolean evaluated_ = false;
    private double[] mutation_sts_;
    private boolean adaptive_;

    public Individual(double[] values, boolean adaptive) {
        values_ = values;
        mutation_sts_ = new double[10];
        Arrays.fill(mutation_sts_, 1.5);
        adaptive_ = adaptive;
    }

    public double getFitness() {
        return fitness_;
    }

    public double[] getValues() {
        return values_;
    }

    public double[] getStepSizes() {
        return mutation_sts_;
    }

    public void setFitness(double fitness) {
        fitness_ = fitness;
        evaluated_ = true;
    }

    public void setValues(double values[]) {
        values_ = values;
    }

    public void setStepSizes(double[] step_size) {
        mutation_sts_ = step_size;
    }

    public boolean isAdaptive_() {
        return adaptive_;
    }

    public boolean isEvaluated() {
        return evaluated_;
    }

    @Override
    public int compareTo(Individual i) {
        return this.fitness_ < i.getFitness() ? 1 : this.fitness_ > i.getFitness() ? -1 : 0;
    }

}