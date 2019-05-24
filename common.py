# -*- coding: utf-8 -*-
import collections

def binarysearch(l: list, x: float, i: int = None, f : int = None) -> int:
    """Return the position of x in the ordered list l using binary search between the positions i and f of l"""
    if f is None and i is None:
        f = len(l)-1        
        i = 0
    if i>f:
        return -1
    m = (i+f)//2
    if l[m] == x:
        return m
    if l[m]<x:
        return binarysearch(l, x, m+1, f)
    return binarysearch(l, x, i, m-1)

def binarysearch_geq(l:list, x: float, i: int = None, f:int = None) -> int:
    """Return the position of the first element greater than or equal to x in the ordered list l using binary search between the positions i and f of l"""
    if f is None and i is None:
        f = len(l)-1        
        i = 0
    if i>f:
        return -1
    if l[f]<x:
        return -1
    if l[i]>=x:
        return i    
    m = (i+f)//2
    if l[m] == x:
        return m
    if l[m]<x:
        return binarysearch_geq(l, x, m+1, f)
    return binarysearch_geq(l, x, i, m)

def bottom_up(child_tree: list) -> list:    
    """Returns a list with the nodes sortered in bottom-up order"""
    q = collections.deque([[0,0]])
    r = collections.deque()
    while len(q)>0:
        v = q[-1]
        if v[1] == len(child_tree[v[0]]):
            r.append(v[0])
            q.pop()
        else:            
            q.append([child_tree[v[0]][v[1]],0])
            v[1] += 1
    return [v for v in r]

def min_M_b_pairs(p1: tuple, p2: tuple) -> tuple:
    if p1[1]<p2[1]:
        return p1
    if p1[1]>p2[1]:
        return p2
    if p1[0]<=p2[0]: # Notice that the indices of the weights are in increasing order
        return p1
    return p2