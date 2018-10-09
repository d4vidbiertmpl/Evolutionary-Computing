class Individual implements Comparable<Individual> {
    /*
    Class that wraps a single individual and its fitness.
    Enables easy comparison among individuals.
    */

    public double[] values_;
    public double fitness_;
    private boolean evaluated_ = false;

    public Individual(double values[]) {
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
        evaluated_ = true;
    }

    public void setValues(double values[]) {
        values_ = values;
    }

    public boolean isEvaluated(){
      return evaluated_;
    }

    @Override
    public int compareTo(Individual i) {
        return this.fitness_ < i.getFitness() ? 1 : this.fitness_ > i.getFitness() ? -1 : 0;
    }
}
