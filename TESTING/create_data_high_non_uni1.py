#!/usr/bin/python
import random
import math

num_data = 28000
max_leaf = 3000
exp_terms = 20
loc_terms = exp_terms

data_box_size = 1000.0


f = open('fmm_data.conf','w')
f.write("0\n")

# 1/5 uniformly in square

for i in range (1*num_data/5):
    x = random.random()*data_box_size
    y = random.random()*data_box_size
    init_val = 0.1 #random.random()/num_data
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)
    
#2/5 in a small circle
for i in range(2*num_data/5):
    invalid = True
    x = 0
    y = 0
    while(invalid) :
        x = random.random() - 0.5
        y = random.random() - 0.5
        if x*x + y*y > 0.25:
            continue
        invalid = False
        x *= 0.003*data_box_size
        y *= 0.003*data_box_size
        x += data_box_size/2;
        y += data_box_size/2;
    init_val = 0.1 #random.random()/num_data
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)        
        
#2/5 in a circle with 1/r^2 density
# thanks to http://www.twam.info/software/non-uniform-distributed-random-numbers
for i in range(2*num_data/5):
    invalid = True
    x = 0
    y = 0
    radius = data_box_size/2
    while(invalid) :
        x = 2*(random.random() - 0.5)
        y = 2*(random.random() - 0.5)
        if x*x + y*y > 1:
            continue
        invalid = False
        sample_len = x*x + y*y;
        #adjust sample_len to 1/r^2 distribution:
        smp_min = 2.0
        smp_max = 2.02
        smp_fac = 1/(smp_max-smp_min)        
        
        sample_len = radius*smp_fac*(1/(sample_len*(1/smp_max - 1/smp_min)+1/smp_min) - smp_min)
        #adjust coordinates to their length
        x *= sample_len
        y *= sample_len
        #move them to center of box
        x += data_box_size/2
        y += data_box_size/2
    init_val = 0.1 #random.random()/num_data
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)       
        
##targets
# 1/5 uniformly in square
for i in range (1*num_data/5):
    x = random.random()*data_box_size
    y = random.random()*data_box_size
    init_val = 0.1 #random.random()/num_data
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)
    
#2/5 in a small square
for i in range(2*num_data/5):
    invalid = True
    x = 0
    y = 0
    while(invalid) :
        x = random.random() - 0.5
        y = random.random() - 0.5
        if x*x + y*y > 0.25:
            continue
        invalid = False
        x *= 0.003*data_box_size
        y *= 0.003*data_box_size
        x += data_box_size/2;
        y += data_box_size/2;
    init_val = 0.1 #random.random()/num_data
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)        
        
#2/5 in a square with density proportional to 1/r^2
for i in range(2*num_data/5):
    invalid = True
    x = 0
    y = 0
    radius = data_box_size/2
    while(invalid) :
        x = 2*(random.random() - 0.5)
        y = 2*(random.random() - 0.5)
        if x*x + y*y > 1:
            continue
        invalid = False
        sample_len = x*x + y*y;
        #adjust sample_len to 1/r^2 distribution:
        smp_min = 2.0
        smp_max = 2.02
        smp_fac = 1/(smp_max-smp_min)        
        
        sample_len = radius*smp_fac*(1/(sample_len*(1/(smp_max) - 1/smp_min)+1/smp_min) - smp_min)
        #adjust coordinates to their length
        x *= sample_len
        y *= sample_len
        #move them to center of box
        x += data_box_size/2
        y += data_box_size/2
    init_val = 0.1 #random.random()/num_data
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)  
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();
