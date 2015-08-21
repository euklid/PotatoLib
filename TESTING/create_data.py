#!/usr/bin/python
import random

num_data = 1000
max_leaf = 10
exp_terms = 15
loc_terms = exp_terms


f = open('fmm_data.conf','w')
f.write("0\n")

for i in range(num_data):
    x = random.gauss(3,1)
    if x < 0:
        x = 0
    if x > 6:
        x = 0
    y = random.gauss(3,1)
    if y <0:
        y = 0
    if y > 6:
        y = 0
    init_val = random.random()
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)
    
for i in range(num_data):
    x = random.gauss(3,1)
    if x < 0:
        x = 0
    if x > 6:
        x = 0
    y = random.gauss(3,1)
    if y <0:
        y = 0
    if y > 6:
        y = 0
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " 0" + "\n"
    f.write(line) 
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();
