# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import  pairwise_distances
from sklearn.cluster import spectral_clustering
import main

a = np.pi/4

points = [(np.cos(a*i)+4, np.sin(a*i)+7) for i in range(8)]

a = np.pi/40

points = points+[(np.cos(a*i)*6, np.sin(a*i)*6) for i in range(-40,-66, -1)]

a = np.pi/8
points = points + [(np.cos(a*i)*7-10, np.sin(a*i)*7+1) for i in range(0, -10,-1)]

D = pairwise_distances(points, points)
D = D / np.max(D)
s = 1 - D
s = np.triu(s)

k = 3
q, k_clustering = main.getclustering(s,k)

markers = ['o','>','s']
colors = [(0.2, 0.5, 0.9), 'k', 'r']

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
for i in range(len(k_clustering)):
    plt.scatter([points[j][0] for j in k_clustering[i]], [points[j][1] for j in k_clustering[i]], 
                marker = markers[i],
                edgecolor = colors[i],
                facecolor='white')
plt.title('proposed algorithm')
    
    
plt.subplot(1, 2, 2)
normcut_labels = spectral_clustering(s, n_clusters=k)
for i in range(len(normcut_labels)):
    plt.scatter(points[i][0], points[i][1],
                marker = markers[normcut_labels[i]],
                edgecolor = colors[normcut_labels[i]],
                facecolor='white')
plt.title('normalized cut')
plt.tight_layout()
plt.show()
