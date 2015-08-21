#!/usr/bin/python
import random

num_data = 50
max_leaf = 12
exp_terms = 6
loc_terms = exp_terms

data_box_size = 60000.0


f = open('fmm_data.conf','w')
f.write("0\n")

for i in range(num_data):
    invalid = True
    x = 0
    y = 0
    while(invalid) :
        x = random.gauss(data_box_size/2,data_box_size/4)
        if x < 0:
            continue
        if x > data_box_size:
            continue
        y = random.gauss(data_box_size/2,data_box_size/4)
        if y <0:
            continue
        if y > data_box_size:
            continue
        invalid = False
    init_val = 0.01 #random.random()/num_data
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)

    
for i in range(num_data):
    invalid = True
    x = 0
    y = 0
    while(invalid) :
        x = random.gauss(data_box_size/2,data_box_size/4)
        if x < 0:
            continue
        if x > data_box_size:
            continue
        y = random.gauss(data_box_size/2,data_box_size/4)
        if y <0:
            continue
        if y > data_box_size:
            continue
        invalid = False
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " 0" + "\n"
    f.write(line) 
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();
