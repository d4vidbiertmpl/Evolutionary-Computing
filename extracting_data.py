import time
import os
import itertools
import numpy as np
import random
import csv

import pickle


def extract_data(link):
    os.system(
        "jar cmf MainClass.txt submission.jar player31.class Parameters.class Individual.class EA_Utils.class Clustering_Utils.class")
    os.system("javac -cp contest.jar player31.java Parameters.java Individual.java EA_Utils.java Clustering_Utils.java")
    os.system(
        "jar cmf MainClass.txt submission.jar player31.class Parameters.class Individual.class EA_Utils.class Clustering_Utils.class")
    os.system("javac -cp contest.jar player31.java Parameters.java Individual.java EA_Utils.java Clustering_Utils.java")

    print "HERE PYTHON SCRIPT                         "

    entire_runs = []

    for i in range(0, 100):

        seed = random.randint(0, 10000000)

        # decide on which function to evaluate
        command = 'java -jar testrun.jar -submission=player31  -evaluation=BentCigarFunction -seed=' + str(seed)
        # command = 'java -jar testrun.jar -submission=player31  -evaluation=SchaffersEvaluation -seed=' + str(seed)
        # command = 'java -jar testrun.jar -submission=player31  -evaluation=KatsuuraEvaluation -seed=' + str(seed)

        output = os.popen(command).readlines()

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

        entire_runs.append(numpy_populations)

    with open(link, 'wb') as f:
        pickle.dump(entire_runs, f)

    # with open("saved_runs/NAME_OF_APPROACH_AND_FUNC.csv", "w") as output:
    #     writer = csv.writer(output, lineterminator='\n')
    #     writer.writerows(entire_runs)

    print entire_runs

    print "success"


# dummy function for loading the data, just copy it where you need
def load_data(link):
    # data = list(csv.reader(open("saved_runs/NAME_OF_APPROACH_AND_FUNC.csv", "rU"), delimiter=','))

    with open(link, 'rb') as f:
        data = pickle.load(f)

    print len(data)
    print data[0].shape


if __name__ == "__main__":

    # before running check where the data will be saved
    link = 'saved_runs/NAME_OF_APPROACH_AND_FUNC.pkl'
    extract_data(link)
    load_data(link)
