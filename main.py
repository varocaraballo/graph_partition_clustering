# -*- coding: utf-8 -*-
import numpy as np
import collections
from common import binarysearch
from common import bottom_up
import fixed_k
import estimate_k

def getclusters(sm: list, k: int = None)->(float, list):
    """Given an NxN normalized similarity matrix computes a clustering of the elements {0,...,N-1}.
    
        Parameters:        
        ----------
        sm : list of lists of floats
            NxN (upper triangular) normalized similarity matrix, i.e.: 0 < sm[i][j] < 1 is the similarity between the elements i and j. If sm[i][j] --> 0 then i and j are very different. If sm[i][j] --> 1 then i and j are very similar.        
        k : int, optional
            Number of required clusters, 2<k<=N.
        
        Returns:
        --------
        C : set
            C is a set of k sets of integers, where every set in C is a cluster, i.e.: assuming as input a 5x5 matrix and k=2, a possible returns is {{0,1,4},{2,3}}.
            If k is not given the number of sets in C is the value of k that minimizes P.
        """
    from scipy.sparse import csr_matrix, lil_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree 
    n = len(sm)
    M = csr_matrix(sm)
    s = -minimum_spanning_tree(-M)
    parent = [None]*n
    child_tree = [collections.deque() for i in range(n)]
    w_set = set(s.data)
    w_set.add(0)
    w_set.add(1)
    weights = sorted([i for i in w_set]) # weights
    edges = lil_matrix(np.zeros((n,n)), dtype=int) #edges 
    q = collections.deque([0])
    marks = [0]*n
    
    while len(q)>0:
        p = q.popleft()
        marks[p] = 1
        for i in range(n):
            e_weight = s[min(i,p), max(i,p)]
            if e_weight > 0 and marks[i] == 0:
                edges[i,p] = binarysearch(weights, e_weight)  
                edges[p,i] = edges[i,p]
                child_tree[p].append(i)
                parent[i] = p
                q.append(i)
    edges = csr_matrix(edges)
    child_tree = [[i for i in l] for l in child_tree]
    if k is not None:        
        trees_tables = {}
        trees_l_mu = {}
        for v in bottom_up(child_tree):
            if len(child_tree[v]) == 0:
                trees_tables[v] = {(1,len(weights)-1):(0,0, None, None)} # (l,mu):(M, b, depends on tree_key, key (l',mu') of the entry on which it depends)
                trees_l_mu[v] = {1:[len(weights)-1]}
            else:
                s_table = s_l_mu = None
                for i in range(len(child_tree[v])):
                    _q_table, _q_l_mu = fixed_k.up_to_parent(trees_tables[child_tree[v][i]], trees_l_mu[child_tree[v][i]], child_tree[v][i], weights, edges[v,child_tree[v][i]], k)
                    if i == 0:
                        s_table, s_l_mu = _q_table, _q_l_mu
                    else:
                        trees_tables[(v,0,i)], trees_l_mu[(v,0,i)] = s_table, s_l_mu
                        trees_tables[(v,i,1)], trees_l_mu[(v,i,1)] = _q_table, _q_l_mu
                        s_table, s_l_mu = fixed_k.add_child(s_table, s_l_mu, (v,0,i), _q_table, _q_l_mu, (v,i,1), weights, edges[v,child_tree[v][i]], k)
                trees_tables[v], trees_l_mu[v] = s_table, s_l_mu
        root_table = trees_tables[0]                   
        root_l_mu = trees_l_mu[0]
        best = float('inf')
        best_key = None 
        for mu in root_l_mu[k]:
            if root_table[(k,mu)][1]<best:
                best = root_table[(k,mu)][1]
                best_key = (k,mu)        
        return best, list(get_clustering(trees_tables, 0, best_key, len(weights)-1))
    else:
        trees_tables = {}
        trees_l_mu = {}
        for v in bottom_up(child_tree):
            if len(child_tree[v]) == 0:
                trees_tables[v] = {(0,len(weights)-1):(0,0, None, None)} # (l,mu):(M, b, depends on tree_key, key (l',mu') of the entry on which it depends)
                trees_l_mu[v] = {0:[len(weights)-1]}
            else:
                s_table = s_l_mu = None
                for i in range(len(child_tree[v])):
                    _q_table, _q_l_mu = estimate_k.up_to_parent(trees_tables[child_tree[v][i]], trees_l_mu[child_tree[v][i]], child_tree[v][i], weights, edges[v,child_tree[v][i]])
                    if i == 0:
                        s_table, s_l_mu = _q_table, _q_l_mu
                    else:
                        trees_tables[(v,0,i)], trees_l_mu[(v,0,i)] = s_table, s_l_mu
                        trees_tables[(v,i,1)], trees_l_mu[(v,i,1)] = _q_table, _q_l_mu
                        s_table, s_l_mu = estimate_k.add_child(s_table, s_l_mu, (v,0,i), _q_table, _q_l_mu, (v,i,1), weights, edges[v,child_tree[v][i]])
                trees_tables[v], trees_l_mu[v] = s_table, s_l_mu
        root_table = trees_tables[0]                   
        root_l_mu = trees_l_mu[0]
        best = float('inf')
        best_key = None 
        for mu in root_l_mu[1]:
            if root_table[(1,mu)][1]<best:
                best = root_table[(1,mu)][1]
                best_key = (1,mu)        
        return best, list(get_clustering(trees_tables, 0, best_key, len(weights)-1))
    

def get_clustering(trees_tables: dict, t_key: int or tuple, p_key: tuple, edge_1: int, s: set = None) -> collections.deque:
    l = collections.deque()
    if s is None:
        s = set()
        l.append(s)
    if type(t_key) is int:
        s.add(t_key)                
    tpl = trees_tables[t_key][p_key]
    if tpl[2] is not None:
        if tpl[2][0] == 1:
            if p_key[1] == edge_1:
                l.extend(get_clustering(trees_tables, tpl[2][1], tpl[3], edge_1))
            else:
                l.extend(get_clustering(trees_tables, tpl[2][1], tpl[3], edge_1, s))
        else:
            l.extend(get_clustering(trees_tables, tpl[2][1], tpl[3][0], edge_1, s))
            l.extend(get_clustering(trees_tables, tpl[2][2], tpl[3][1], edge_1, s))
    return l

    

    
    
    
    
    
        
        