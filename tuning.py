import os
import random

seed = random.randint(0,10000000)
# run the java code with a new uniform_mutation_prop
os.system("java -Dvar1=0.5 -Dnon_uniform_mutation_step_size=2 -jar testrun.jar -submission=player31 -evaluation=BentCigarFunction -seed="+str(seed))
