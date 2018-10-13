import os
import random
import numpy as np

"""
This implementation follows the decription of the REVAC optimizer in
https://www.ijcai.org/Proceedings/07/Papers/157.pdf .
"""

def build_command(parameter_values, parameter_names):
    """
    Builds the command to run the Java application according to given parameters.
    """
    # check if same amount of values and fields
    assert len(parameter_values) == len(parameter_names)

    command = "java "
    for name, value in zip(parameter_names, parameter_values):
        command+="-D"+name+"="+str(value)+" "
    command += "-jar testrun.jar -submission=player31 -evaluation=BentCigarFunction -seed="
    seed = random.randint(0,10000000)
    return command+str(seed);

def log_vector(vector, response, evaluation=0, path="revac_output.txt"):
    """
    Logs a parameter vector and its response to path.
    """
    with open(path, 'a') as f:
        f.write("["+str(evaluation)+"\t"+str(vector)+"\t"+str(response)+"] \n")

def evaluate_parameters(parameter_vector, parameter_names, num_runs):
    """
    Evaluates parameter_vector by computings its average performance over num_runs runs.
    """
    score =  0.0
    for r in range(num_runs):
        x = os.popen(build_command(parameter_vector, parameter_names)).readlines()
        score += float([i for i in x if 'Score:' in i][0].split()[1])
    return score/num_runs

def multi_parent_crossover(parents):
    """
    Performs multi-parent crossover with uniform scanning.
    """
    child=[]
    for i in range(len(parents[0])):
        child += [parents[random.randint(0, len(parents)-1)][i]]
    return child

def revac_mutation(child, parents, h, parameters_info):
    """
    Performs the mutation described in the REVAC paper.
    """
    mutated_child = []

    # for each parameter:
    for index, param_info in enumerate(parameters_info):

        # Sort the values of this parameter for all vectors into ascending order
        p_values = [p[index] for p in parents]
        p_values.sort()

        # compute upper and lower bounds of the mutation distribution
        c_ind = p_values.index(child[index])
        max_ind = c_ind+h
        min_ind = c_ind-h
        # catch cases in that c_ind is very in beginning/end of list:
        if min_ind<0:
            min_ind=0
        if max_ind>=len(p_values):
            max_ind=len(p_values)-1

        # drawing a new parameter value at uniform random
        if param_info[2]=="float":
            mutated_child += [np.random.uniform(p_values[min_ind], p_values[max_ind])]
        if param_info[2]=="int":
            mutated_child += [np.random.randint(p_values[min_ind], p_values[max_ind]+1)]

    return mutated_child


def tune(parameters, population_size=80, parent_size=35, h=5, runs_per_vector=1):
    """
    Tunes a given algorithm using REVAC optimization.
    """
    parameter_names =  [p[0] for p in parameters]

    # (1) Draw an initial set of parameter vectors at uniform random from
    # their initial distributions.
    population = []
    for i in range(population_size):
        vector = []
        for p in parameters:
            assert p[2]=="float"or p[2]=="int"
            if p[2]=="float":
                vector += [np.random.uniform(p[1][0], p[1][1])]
            if p[2]=="int":
                vector += [np.random.randint(p[1][0], p[1][1]+1)]
        population += [vector]

    # (2) Compute the utility of each parameter vector, before finding and recording the best parameter vector
    utility_table = [evaluate_parameters(parameter_vector, parameter_names,runs_per_vector) for parameter_vector in population]

    # log best parameters
    max_utility = max(utility_table)
    best_parameter_vector = population[utility_table.index(max_utility)]
    log_vector(best_parameter_vector, max_utility)

    # TODO: Different termination condition?
    evaluations=0
    oldest_index=0

    while(evaluations<1000):
        # (3) Select the parent_size-best vectors from the table as the parents
        # of the next parameter vector.
        parent_indices = np.argpartition(utility_table, -parent_size)[-parent_size:]
        parents = [population[i] for i in parent_indices]

        # (4) Perform multi-parent crossover on the N parents to create
        # a proto-child vector.
        child = multi_parent_crossover(parents)

        # (5) Mutate Child vector with regard to parents
        child = revac_mutation(child, parents, h, parameters)

        # Replace oldest member of population by child vector
        population[oldest_index] = child.copy()
        #Update utility table accoriding to new member
        child_utility = evaluate_parameters(child, parameter_names, runs_per_vector)
        utility_table[oldest_index] = child_utility

        # Update the best vector if the child vector is an improvement.
        if child_utility > max_utility:
            max_utility = child_utility
            best_parameter_vector = child.copy()

        # Update Statistics
        oldest_index = (oldest_index+1)%population_size
        evaluations += 1
        log_vector(child, child_utility, evaluation=evaluations)

    log_vector(best_parameter_vector, max_utility, evaluation="BEST")
    # Return optimal parameter vectore
    return best_parameter_vector


if __name__ == "__main__":

    # Set Parameters you want to tune here. If you don't want them to be tuned, don't put them in the parameters list
    # Define a reasonable range for the values of the parameters.
    parameters=[["elitist_size", [5,25], "int"],
                ["proletarian_size", [10,30], "int"],
                ["cluster_distance_thresh", [0.2,3], "float"],
                ["offspring_size", [50,100], "int"],
                ["parent_tournament_size",[2,50], "int"],
                ["survivor_tournament_size",[2,50], "int"],
                ["non_uniform_mutation_step_size", [0.05, 3], "float"],
                ["hill_climb_step_size", [0.1, 2], "float"],
                ["evaluations_per_proletarian", [40, 100], "int"]]

    #params of the revac opimizer
    population_size=100
    parent_size=50
    h=int(round(parent_size/10)) # mentioned in paper
    runs_per_vector = 5 # Not too big 

    # tune parameters with revac
    tune(parameters, population_size, parent_size, h,runs_per_vector)
