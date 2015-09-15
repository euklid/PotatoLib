#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 02:10:04 2015

@author: euklid
"""

import random
import math

num_data = 512000
max_leaf = 30
exp_terms = 17
loc_terms = exp_terms

data_box_size = 1.0


f = open('fmm_data_spiro_' + str(num_data) + '.conf','w')
f.write("0\n")
radius = data_box_size/3
rings = 7
shift = math.pi/rings;


for i in range(num_data):
    x = radius*math.cos(float(i)/num_data*2*math.pi) + radius/2*math.cos(rings*float(i)/num_data*2*math.pi) + data_box_size/2
    y = radius*math.sin(float(i)/num_data*2*math.pi) + radius/2*math.sin(rings*float(i)/num_data*2*math.pi) + data_box_size/2
    init_val = random.random()
    s = 1
    t = 0
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " " + str(init_val)+"\n"
    f.write(line)

    
for i in range(num_data):
    x = radius*math.cos((float(i)+0.5)/num_data*2*math.pi + shift) + radius/2*math.cos(rings*(float(i)+0.5)/num_data*2*math.pi) + data_box_size/2
    y = radius*math.sin((float(i)+0.5)/num_data*2*math.pi + shift) + radius/2*math.sin(rings*(float(i)+0.5)/num_data*2*math.pi) + data_box_size/2
    s = 0
    t = 1
    line = str(x) + " " + str(y) + " " + str(s) + " " + str(t) + " 0" + "\n"
    f.write(line) 
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();