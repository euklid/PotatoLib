#!/usr/bin/python
import random
import math

num_data = 1000000
max_leaf = 30
exp_terms = 20
loc_terms = exp_terms

data_box_size = 6000000.0


f = open('fmm_data.conf','w')
f.write("0\n")
radius = data_box_size/2

for i in range(num_data):
    x = radius*math.cos(float(i)/num_data*2*math.pi)
    y = radius*math.sin(float(i)/num_data*2*math.pi)
    init_val = 0.1 #random.random()/num_data
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)

    
for i in range(num_data):
    x = radius*math.cos((float(i)+0.5)/num_data*2*math.pi)
    y = radius*math.sin((float(i)+0.5)/num_data*2*math.pi) 
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " 0" + "\n"
    f.write(line) 
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();
