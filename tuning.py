import os
import random
import numpy as np

"""
This implementation follows the decription of the REVAC optimizer in
https://www.ijcai.org/Proceedings/07/Papers/157.pdf .
"""

def testDiversity():
    x = os.popen(build_command([],[])).readlines()
    diversities = [i.split()[1] for i in x if 'Diversity:' in i]
    for d in diversities:
        with open("divs.txt", 'a') as f:
            f.write(d+"\n")

def build_command(parameter_values, parameter_names):
    """
    Builds the command to run the Java application according to given parameters.
    """
    # check if same amount of values and fields
    assert len(parameter_values) == len(parameter_names)
    command = "java "
    for name, value in zip(parameter_names, parameter_values):
        if name=="offspring_percentage":
            pop_size = parameter_values[parameter_names.index("population_size")];
            command+="-Doffspring_size="+str(round(value*pop_size/100))+" "
        else:
            command+="-D"+name+"="+str(value)+" "
    command += "-jar testrun.jar -submission=player31 -evaluation=KatsuuraEvaluation -seed="
    seed = random.randint(0,10000000)
    return command+str(seed);

def log_vector(vector, response, evaluation=0, path="revac_output.txt"):
    """
    Logs a parameter vector and its response to path.
    """
    with open(path, 'a') as f:
        f.write("["+str(evaluation)+"\t"+str(vector)+"\t"+str(response)+"] \n")

def evaluate_parameters(parameter_vector, parameter_names, num_runs, ini=false):
    """
    Evaluates parameter_vector by computings its average performance over num_runs runs.
    """
    score =  0.0
    for r in range(num_runs):
        command = build_command(parameter_vector, parameter_names)
        x = os.popen(command).readlines()
        utility =  float([i for i in x if 'Score:' in i][0].split()[1])
        score += utility
    if ini:
        log_vector(parameter_vector, score/num_runs, evaluation="RANDOM:  ")
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


def tune(parameters, population_size=80, parent_size=35, h=5, runs_per_vector=1, total_evaluations=1000):
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
                vector += [random.randint(p[1][0], p[1][1])]
        population += [vector]
    #draw initial samples with latin hypercube sampling
    #population = latin_hypercube_sampling(parameters)
    #print("Initialized. Start Evaluation")

    # (2) Compute the utility of each parameter vector, before finding and recording the best parameter vector
    utility_table = [evaluate_parameters(parameter_vector, parameter_names,runs_per_vector, ini=true) for parameter_vector in population]

    # log best parameters
    max_utility = max(utility_table)
    best_parameter_vector = population[utility_table.index(max_utility)]
    log_vector(best_parameter_vector, max_utility)

    evaluations=0
    oldest_index=0

    while(evaluations<total_evaluations):
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

def latin_hypercube_sampling(parameters, k=100):
    # divide range of each dimension into k equally sized intervals
    ranges_parameters = []
    for p in parameters:
        rng = p[1]
        x = np.linspace(rng[0],rng[1], k+1)
        if p[2]=="int":
            x = [round(i) for i in x]
        # ranges for parameters
        y = []
        for i in range(len(x)-1):
            y += [[x[i], x[i+1]]]
        ranges_parameters += [y]

    # For each new design point: map a range in each dimension one–to–one + draw a random value within it
    [random.shuffle(i) for i in ranges_parameters] # shuffling once = drawing k times random
    new_vectors=[]
    for i in range(k):
        new_vector = []
        for j in range(len(parameters)):
            rng = ranges_parameters[j][i]
            if parameters[j][2]=="float":
                new_vector += [random.uniform(rng[0], rng[1])]
            if parameters[j][2]=="int":
                new_vector += [random.randint(rng[0], rng[1])]
        new_vectors += [new_vector]
    return new_vectors


if __name__ == "__main__":

    # Set Parameters you want to tune here. If you don't want them to be tuned, don't put them in the parameters list
    # Define a reasonable range for the values of the parameters.

    ## PARAMS TO TUNE:
    # (1) Simple approach
    # On Katsuura:
    """
    parameters=[["offspring_percentage", [100,150], "float"],
                ["parent_tournament_size", [4,50], "int"],
                ["survivor_tournament_size", [2,50], "int"],
                ["non_uniform_mutation_step_size", [0.05, 3], "float"],
                ["population_size", [1000,1000], "int"]]
    """
    #  Else:
    """
    parameters=[["offspring_percentage", [100,150], "float"],
                ["parent_tournament_size", [4,50], "int"],
                ["survivor_tournament_size", [2,50], "int"],
                ["non_uniform_mutation_step_size", [0.05, 3], "float"],
                ["population_size", [100,100], "int"]]
    """

    # (2) Simple approach + Ours
    # On Katsuura:
    """
    parameters=[["offspring_percentage", [100,150], "float"],
                ["parent_tournament_size", [4,50], "int"],
                ["survivor_tournament_size", [2,50], "int"],
                ["non_uniform_mutation_step_size", [0.05, 3], "float"],
                ["population_size", [1000,1000], "int"],
                ["cluster_distance_thresh", [1.2,3], "float"],
                ["hill_climb_step_size", [0.1, 1], "float"]]
    """
    # Else:
    """
    parameters=[["offspring_percentage", [100,150], "float"],
                ["parent_tournament_size", [4,50], "int"],
                ["survivor_tournament_size", [2,50], "int"],
                ["non_uniform_mutation_step_size", [0.05, 3], "float"],
                ["population_size", [100,100], "int"],
                ["cluster_distance_thresh", [1.2,3], "float"],
                ["hill_climb_step_size", [0.1, 1], "float"]]
    """

    # (3) Sophisticated Approach:
    # On Katsuura
    """
    parameters=[["parent_tournament_size", [4,50], "int"],
                ["population_size", [1000,1000], "int"]]
    """
    # Else
    """
    parameters=[["parent_tournament_size", [4,50], "int"],
                ["population_size", [100,100], "int"]]
    """

    # (4) Sophisticated Approach + Ours:
    # On Katsuura
    """
    parameters=[["parent_tournament_size", [4,50], "int"],
                ["population_size", [1000,1000], "int"],
                ["cluster_distance_thresh", [1.2,3], "float"],
                ["hill_climb_step_size", [0.1, 1], "float"]]
    """
    # Else
    """
    parameters=[["parent_tournament_size", [4,50], "int"],
                ["population_size", [100,100], "int"],
                ["cluster_distance_thresh", [1.2,3], "float"],
                ["hill_climb_step_size", [0.1, 1], "float"]]
    """

    parameters=[["parent_tournament_size", [4,50], "int"],
                ["population_size", [1000,1000], "int"],
                ["cluster_distance_thresh", [1.2,3], "float"],
                ["hill_climb_step_size", [0.1, 1], "float"]]


    ## REVAC PARAMS
    ##### (1) USE THESE IF YOU'RE OPTIMIZING ON KATSUURA
    population_size=75
    parent_size=30
    h=3 # about parent_size/10
    runs_per_vector = 5
    num_evaluations = 1000

    ##### (2) USE THESE ELSE:
    #population_size=100
    #parent_size=50
    #h=5 # about parent_size/10
    #runs_per_vector = 20
    #num_evaluations = 1000


    # tune parameters with revac
    tune(parameters, population_size, parent_size, h,runs_per_vector, num_evaluations)
