# -*- coding: utf-8 -*-
from common import binarysearch_geq, min_M_b_pairs

default_pair = (float('inf'), float('inf'), None, None)

def up_to_parent(s_table: dict, s_l_mu: dict, s_key: int or tuple, weights: list, w_index: int) -> (dict, dict):
    _s_table = {} # Corresponding table to the resultant tree of adding the parent of the root node of S
    _s_l_mu = {}    
    used_mu = [set(), set()]
    for l in range(2):
        if l in s_l_mu:        
            for mu in s_l_mu[l]:
                tpl = s_table[(l,mu)]            
                if mu < w_index:   
                    # Computing p_tree_table[(l,mu)] when mu < w_index (by using Claim 19 of the draft)
                    _s_table[(l,mu)] = (tpl[0], tpl[1], (1, s_key), (l, mu))
                    used_mu[l].add(mu)                
                    # End of computing p_tree_table[(l,mu)] when mu < w_index
                else:                
                    # Computing p_tree_table[(l,w_index)] when mu >= w_index (by using Claim 18 of the draft)
                    b = max(weights[tpl[0]]/weights[w_index], tpl[1])
                    if b <=1:
                        _s_table[(l, w_index)] = (tpl[0], b, (1, s_key), (l, mu))
                        used_mu[l].add(w_index)
                    break
                    # End of computing p_tree_table[(l,w_index)]
        if l in s_l_mu:
            # Computing p_tree_table[(l,len(weigths)-1)] (by using Claim 15 of the draft)
            mus = s_l_mu[l]
            pos = binarysearch_geq(mus, w_index)
            if pos>-1:
                for i in range(pos, len(mus)):
                    mu = mus[i]
                    tpl = s_table[(l,mu)]   
                    b = max(weights[w_index], tpl[1], max(weights[w_index],weights[tpl[0]])/weights[mu])
                    if b<=1: 
                        _s_table[(1, len(weights)-1)] = min_M_b_pairs(_s_table.setdefault((1, len(weights)-1), default_pair), 
                                  (w_index, b, (1, s_key), (l, mu)))
                        used_mu[1].add(len(weights)-1)
            # End of computing p_tree_table[(l,len(weigths)-1)]
    for l in range(2):
        if len(used_mu[l])>0:
            _s_l_mu[l] = sorted(list(used_mu[l]))
    return _s_table, _s_l_mu
                
def add_child(s_table: dict, s_l_mu: dict, s_key: int or tuple, _q_table: dict, _q_l_mu: dict, _q_key: int or tuple, weights: list, w_index: int) -> (dict, dict):
    p_table = {} # Corresponding table to the resultant tree of adding the tree q to s
    p_l_mu = {}       
    used_mu = [set(), set()]
    for l_s in range(2):
        for l_q in range(2):
            if l_s not in s_l_mu or l_q not in _q_l_mu:
                continue
            s_mus = s_l_mu[l_s]
            _q_mus = _q_l_mu[l_q]
            # Computing the values assuming that the minimum inside head cluster is apported by tree s
            _mu_index = -1 if len(s_mus)==0 else binarysearch_geq(_q_mus, s_mus[0])            
            for mu in s_mus:
                s_tpl = s_table[(l_s,mu)]
                if mu < w_index:
                    # Computing o1 when mu < w_index (by using Claim 26 of the draft)
                    while _mu_index != -1 and _mu_index<len(_q_mus) and _q_mus[_mu_index]<mu:
                        _mu_index += 1
                    if _mu_index == -1 or _mu_index >= len(_q_mus):
                        break
                    _mu = _q_mus[_mu_index]
                    _q_tpl = _q_table[(l_q, _mu)]
                    M = max(s_tpl[0], _q_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(int(l_s+l_q > 0), mu)] = min_M_b_pairs(p_table.setdefault((int(l_s+l_q > 0), mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((l_s,mu),(l_q,_mu))))
                        used_mu[int(l_s+l_q > 0)].add(mu)
                    # End of computing o1 when mu < w_index
                else:
                    # Computing o1 when mu >= w_index (by using Claims 22 and 23 of the draft)
                    if len(_q_mus)==0 or _q_mus[-1] != len(weights)-1:
                        break
                    _mu = _q_mus[-1]
                    _q_tpl = _q_table[(l_q, _mu)]
                    M = max(w_index, s_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(int(l_s+l_q > 0), mu)] = min_M_b_pairs(p_table.setdefault((int(l_s+l_q > 0), mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((l_s,mu),(l_q, _mu))))
                        used_mu[int(l_s+l_q > 0)].add(mu)
                    # End of computing o1 when mu >= w_index
            # End of computing the values assuming that the minimum inside head cluster is apported by tree s
            
            # Computing the values assuming that the minimum inside head cluster is apported by tree _q
            _mu_index = -1 if len(_q_mus)==0 else binarysearch_geq(s_mus, _q_mus[0])            
            for mu in _q_mus:
                _q_tpl = _q_table[(l_q, mu)]
                if mu <= w_index:
                    # Computing o2 when mu <= w_index (by using Claim 26 and 24 of the draft)
                    while _mu_index != -1 and _mu_index<len(s_mus) and s_mus[_mu_index]<mu:
                        _mu_index += 1
                    if _mu_index == -1 or _mu_index >= len(s_mus):
                        break
                    _mu = s_mus[_mu_index]
                    s_tpl = s_table[(l_s, _mu)]
                    M = max(s_tpl[0], _q_tpl[0])
                    b = max(weights[M]/weights[mu], s_tpl[1], _q_tpl[1])
                    if b <= 1:
                        p_table[(int(l_s+l_q > 0), mu)] = min_M_b_pairs(p_table.setdefault((int(l_s+l_q > 0), mu), default_pair),
                               (M, b, (2, s_key, _q_key), ((l_s, _mu), (l_q, mu))))
                        used_mu[int(l_s+l_q > 0)].add(mu)
                    # End of computing o2 when mu <= w_index
                else:
                    break
            # End of computing the values assuming that the minimum inside head cluster is apported by tree _q
            
            # Computing p_table[(l, len(weights)-1)] (by using Claim 21 of the draft)
            if (l_s, len(weights)-1) in s_table and (l_q, len(weights)-1) in _q_table:
                s_tpl, _q_tpl = s_table[(l_s, len(weights)-1)], _q_table[(l_q, len(weights)-1)]
                b = max(s_tpl[1], _q_tpl[1])
                if b <= 1:                    
                    p_table[(int(l_s+l_q > 0), len(weights)-1)] = min_M_b_pairs(p_table.setdefault((int(l_s+l_q > 0), len(weights)-1), default_pair), 
                             (max(w_index, s_tpl[0]), b, (2, s_key, _q_key), ((l_s, len(weights)-1), (l_q, len(weights)-1))))                    
                    used_mu[int(l_s+l_q > 0)].add(len(weights)-1)
            # End of Computing p_table[(l, len(weights)-1)]
    for l in range(2):
        if len(used_mu[l])>0:
            p_l_mu[l] = sorted(list(used_mu[l]))
    return p_table, p_l_mu
                        
