# -*- coding: utf-8 -*-
import collections
from common import binarysearch_geq, min_M_b_pairs

default_pair = (float('inf'), float('inf'), None, None)

def up_to_parent(s_table: dict, s_l_mu: dict, s_key: int or tuple, weights: list, w_index: int, k: int) -> (dict, dict):
    _s_table = {} # Corresponding table to the resultant tree of adding the parent of the root node of S
    _s_l_mu = {}
    for l in range(1, k+1):
        used_mu = collections.deque()
        if l in s_l_mu:        
            for mu in s_l_mu[l]:
                tpl = s_table[(l,mu)]            
                if mu < w_index:   
                    # Computing p_tree_table[(l,mu)] when mu < w_index (by using Claim 19 of the draft)
                    _s_table[(l,mu)] = (tpl[0], tpl[1], (1, s_key), (l, mu))
                    used_mu.append(mu)                
                    # End of computing p_tree_table[(l,mu)] when mu < w_index
                else:                
                    # Computing p_tree_table[(l,w_index)] when mu >= w_index (by using Claim 18 of the draft)
                    b = max(weights[tpl[0]]/weights[w_index], tpl[1])
                    if b <=1:
                        _s_table[(l, w_index)] = (tpl[0], b, (1, s_key), (l, mu))
                        used_mu.append(w_index)
                    break
                    # End of computing p_tree_table[(l,w_index)]
        if l>1 and l-1 in s_l_mu:
            # Computing p_tree_table[(l,len(weigths)-1)] (by using Claim 15 of the draft)
            mus = s_l_mu[l-1]
            pos = binarysearch_geq(mus, w_index)
            if pos>-1:
                for i in range(pos, len(mus)):
                    mu = mus[i]
                    tpl = s_table[(l-1,mu)]   
                    b = max(weights[w_index], tpl[1], max(weights[w_index],weights[tpl[0]])/weights[mu])
                    if b<=1: 
                        _s_table[(l, len(weights)-1)] = min_M_b_pairs(_s_table.setdefault((l, len(weights)-1), default_pair), 
                                  (w_index, b, (1, s_key), (l-1, mu)))
                        if len(used_mu) == 0 or used_mu[-1] != len(weights)-1:
                            used_mu.append(len(weights)-1)
            # End of computing p_tree_table[(l,len(weigths)-1)]
        if len(used_mu)>0:
            _s_l_mu[l] = [mu for mu in used_mu]
    return _s_table, _s_l_mu
                
def add_child(s_table: dict, s_l_mu: dict, s_key: int or tuple, _q_table: dict, _q_l_mu: dict, _q_key: int or tuple, weights, w_index, k) -> (dict, dict):
    p_table = {} # Corresponding table to the resultant tree of adding the tree q to s
    p_l_mu = {}
    for l in range(1, k+1):
        used_mu = set()
        for x in range(1,l+1):
            if x not in s_l_mu or l-x+1 not in _q_l_mu:
                continue
            s_mus = s_l_mu[x]
            _q_mus = _q_l_mu[l-x+1]
            # Computing the values assuming that the minimum inside head cluster is apported by tree s
            _mu_index = -1 if len(s_mus)==0 else binarysearch_geq(_q_mus, s_mus[0])            
            for mu in s_mus:
                s_tpl = s_table[(x,mu)]
                if mu < w_index:
                    # Computing o1 when mu < w_index (by using Claim 26 of the draft)
                    while _mu_index != -1 and _mu_index<len(_q_mus) and _q_mus[_mu_index]<mu:
                        _mu_index += 1
                    if _mu_index == -1 or _mu_index >= len(_q_mus):
                        break
                    _mu = _q_mus[_mu_index]
                    _q_tpl = _q_table[(l-x+1, _mu)]
                    M = max(s_tpl[0], _q_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(l, mu)] = min_M_b_pairs(p_table.setdefault((l, mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((x,mu),(l-x+1,_mu))))
                        used_mu.add(mu)
                    # End of computing o1 when mu < w_index
                else:
                    # Computing o1 when mu >= w_index (by using Claims 22 and 23 of the draft)
                    if len(_q_mus)==0 or _q_mus[-1] != len(weights)-1:
                        break
                    _mu = _q_mus[-1]
                    _q_tpl = _q_table[(l-x+1, _mu)]
                    M = max(w_index, s_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(l, mu)] = min_M_b_pairs(p_table.setdefault((l, mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((x,mu),(l-x+1, _mu))))
                        used_mu.add(mu)
                    # End of computing o1 when mu >= w_index
            # End of computing the values assuming that the minimum inside head cluster is apported by tree s
            
            # Computing the values assuming that the minimum inside head cluster is apported by tree _q
            _mu_index = -1 if len(_q_mus)==0 else binarysearch_geq(s_mus, _q_mus[0])            
            for mu in _q_mus:
                _q_tpl = _q_table[(l-x+1, mu)]
                if mu <= w_index:
                    # Computing o2 when mu <= w_index (by using Claim 26 and 24 of the draft)
                    while _mu_index != -1 and _mu_index<len(s_mus) and s_mus[_mu_index]<mu:
                        _mu_index += 1
                    if _mu_index == -1 or _mu_index >= len(s_mus):
                        break
                    _mu = s_mus[_mu_index]
                    s_tpl = s_table[(x, _mu)]
                    M = max(s_tpl[0], _q_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(l, mu)] = min_M_b_pairs(p_table.setdefault((l, mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((x, _mu), (l-x+1, mu))))
                        used_mu.add(mu)
                    # End of computing o2 when mu <= w_index
                else:
                    break
            # End of computing the values assuming that the minimum inside head cluster is apported by tree _q
            
            # Computing p_table[(l, len(weights)-1)] (by using Claim 21 of the draft)
            if (x, len(weights)-1) in s_table and (l-x+1, len(weights)-1) in _q_table:
                s_tpl, _q_tpl = s_table[(x, len(weights)-1)], _q_table[(l-x+1, len(weights)-1)]
                b = max(s_tpl[1], _q_tpl[1])
                if b <= 1:                    
                    p_table[(l, len(weights)-1)] = min_M_b_pairs(p_table.setdefault((l, len(weights)-1), default_pair), 
                             (max(w_index, s_tpl[0]), b, (2, s_key, _q_key), ((x, len(weights)-1), (l-x+1, len(weights)-1))))                    
                    used_mu.add(len(weights)-1)
            # End of Computing p_table[(l, len(weights)-1)]
        if len(used_mu)>0:
            p_l_mu[l] = sorted([mu for mu in used_mu])
    return p_table, p_l_mu
                        
