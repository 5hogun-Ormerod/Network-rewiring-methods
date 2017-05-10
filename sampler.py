# -*- coding: utf-8 -*-
"""
Created on Tue May  9 22:21:39 2017

@author: Chris
"""



def fat_model():
    ER = nx.gnm_random_graph(94,203)
    GScale = clust_rewire_3_node_c(ER,5000)
    return clust_rewire_3_loop_a(GScale,100)

def fat_model_b():
    ER = nx.gnm_random_graph(94,203)    
    return clust_rewire_3_loop_b(ER,100)

def mito_model():
    ER = nx.gnm_random_graph(747,2477)
    GScale = clust_rewire_3_node_c(ER,10000)
    return clust_rewire_3_loop_a(GScale,9000)
