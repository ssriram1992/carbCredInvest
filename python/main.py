from EPECinterface import *
from gurobipy import *
import itertools
import csv

class Sets:
	pass

sets = Sets()


sets.dirty = ["coal", "gas"]
sets.green = ["wind", "solar"]

sets.countries = ["c1", "c2"]
sets.producers = ["p1_1", "p1_2", "p2_1", "p2_2"]
sets.domesticity = tuplelist([("c1","p1_1"), ("c1","p1_2"),
               ("c2","p2_1"), ("c2","p2_2")
              ])

sets.scenario = ["s"+str(i+1) for i in range(2)]

mylist = readInputFromFile("inputs.csv")



solution = solveMPECs(sets, mylist)

with open("output2.txt","w") as file:
    for vals in solution:
        file.write(",".join([str(vv) for vv in vals]))
        file.write('\n')
