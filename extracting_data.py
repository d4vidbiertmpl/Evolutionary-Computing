import time
import os
import itertools

import numpy as np

# double because I do not trust the compiling :D
os.system(
    "jar cmf MainClass.txt submission.jar player31.class Parameters.class Individual.class EA_Utils.class Clustering_Utils.class")
os.system("javac -cp contest.jar player31.java Parameters.java Individual.java EA_Utils.java Clustering_Utils.java")
os.system(
    "jar cmf MainClass.txt submission.jar player31.class Parameters.class Individual.class EA_Utils.class Clustering_Utils.class")
os.system("javac -cp contest.jar player31.java Parameters.java Individual.java EA_Utils.java Clustering_Utils.java")
# os.system("java -jar testrun.jar -submission=player31  -evaluation=SchaffersEvaluation -seed=1")

print "HERE PYTHON SCRIPT"

output = os.popen('java -jar testrun.jar -submission=player31  -evaluation=SchaffersEvaluation -seed=1').readlines()

populations = [list(group) for k, group in itertools.groupby(output, lambda x: x.strip() == "+") if not k][:-1]

numpy_populations = np.zeros((0, 100, 10), np.float64)

for population in populations:
    numpy_population = np.zeros((0, 10), np.float64)

    for individual in population:
        individual_list = individual.split(',')
        individual_list[9] = individual_list[9][:-1]
        individual_list_float = [float(i) for i in individual_list]

        numpy_individual = np.reshape(np.asarray(individual_list_float, np.float64), (1, 10))

        numpy_population = np.append(numpy_population, numpy_individual, axis=0)

    numpy_populations = np.append(numpy_populations, np.reshape(numpy_population, (1, 100, 10)), axis=0)

np.save("saved_runs/01_Schaffers_11_10_parameters", np.float64(numpy_populations))

print "success"