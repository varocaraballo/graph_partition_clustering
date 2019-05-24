# -*- coding: utf-8 -*-
import numpy as np
import time
import main

n = 20
k = 3
s = np.random.rand(n,n)
for i in range(n):
    for j in range(i+1):
        s[i,j] = 0

start_time = time.time()
q, k_clustering = main.getclusters(s,k)
elapsed_time = time.time() - start_time
print((q, k_clustering, elapsed_time))


start_time = time.time()
q1, clustering = main.getclusters(s)
elapsed_time = time.time() - start_time
print((q1, clustering, elapsed_time))