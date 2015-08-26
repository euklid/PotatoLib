#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 04:51:34 2015

@author: euklid
"""


import math

num_el = 4000
max_leaf = 30
exp_terms = 15
loc_terms = exp_terms

data_box_size = 2000.0


f = open('fmm_data.conf','w')
f.write("1\n")
f.write(str(num_el)+"\n")

radius = data_box_size/2

for i in range(num_el):
    x = radius*math.cos(float(i)/num_el*2*math.pi)
    y = radius*math.sin(float(i)/num_el*2*math.pi)
    line = str(i) + " " + str(x) + " " + str(y) + "\n"
    f.write(line)

    
for i in range(num_el):
    end_node = (i+1)%num_el
    s = 1
    t = 1
    line = str(i) + " " + str(end_node) + " " + str(s) + " " + str(t) + " "+  "0.5" + "\n"
    f.write(line) 
    
f.close()


f = open('fmm_conf.conf','w')
line = str(exp_terms) + " " + str(loc_terms) + " " + str(max_leaf);
f.write(line);
f.close();
