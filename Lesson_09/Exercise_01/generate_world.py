#!/usr/bin/python

import random
import math

number_of_cities=30
dim=1.

file1 = open("city_config.square","w+")
file2 = open("city_config.circle","w+")

file1.write(str(number_of_cities)+"\n")
file2.write(str(number_of_cities)+"\n")

for i in range(number_of_cities):
    x=random.uniform(-dim,dim)
    y=random.uniform(-dim,dim)
    file1.write(str(x)+"\t"+str(y)+"\n")

for i in range(number_of_cities):
    theta=random.uniform(0, 2.*math.pi)
    x=dim*math.cos(theta)
    y=dim*math.sin(theta)
    file2.write(str(x)+"\t"+str(y)+"\n")

file1.close()
file2.close()
