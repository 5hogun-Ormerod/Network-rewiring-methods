import networkx as nx
import random

def Random_Path(G, length, path = []):  
    """
    Given a graph, G, this function returns a random path of a given number
    of edges (length =n). This path is represented by a list of nodes of G
    [v0,v1,..., vn] such that (vi,vi+1) is an edge. 
    
    This is done recursively by appending a random neighbor of the last
    node in the path then appending a path of length n-1. 
    
    """
    if path == []:
        rand_node = random.choice(nx.nodes(G))
        path.append(rand_node)
        print(rand_node)
    if length == 1:
        prev_node = list(path)[-1] 
        connected = list(set(nx.all_neighbors(G,prev_node)).difference(set(path)))
        if len(connected) > 0:
            random_new_node = random.choice(connected)
            path.append(random_new_node)
            return path
        else:
            return path
    else:
        prev_node = list(path)[-1] 
        connected = list(set(nx.all_neighbors(G,prev_node)).difference(set(path)))
        if len(connected) > 0:
            random_new_node = random.choice(connected)
            path.append(random_new_node)
            return Random_Path(G, length-1, path)
        else:
            return path
            
def rewire_scale_free(G,rewire_iter,max_iter = 10000000):
    """
    We create a scale free graph by choosing a random node, x, and a 
    random path [xij] and rewiring [xi] to [xj] sufficiently many times.   
    """
    G_new = nx.Graph(G)
    count_iter = 0
    count_rewire = 0
    while (count_iter < max_iter) & (count_rewire < rewire_iter):
        count_iter = count_iter + 1
        rand_node = random.choice(G_new.nodes())
        path = Random_Path(G_new, 2 ,[rand_node])      
        if len(path) ==  3:
            if path[2] in nx.non_neighbors(G_new,path[0]):
                G_new.remove_edge(path[0],path[1])
                G_new.add_edge(path[0],path[2])
                count_rewire = count_rewire +1
    return G_new
    
def rewire_cluster_binomial(G,rewire_iter,max_iter = 10000000):
    G_new = nx.Graph(G) 
    count_iter = 0
    count_rewire = 0
    while (count_iter < max_iter) & (count_rewire < rewire_iter):
        count_iter = count_iter + 1
        rand_node = random.choice(G_new.nodes())
        if len(list(nx.all_neighbors(G_new,rand_node))) >= 2:
            path = Random_Path(G_new, 3,[rand_node])
            if len(path)== 4:
                if path[0] in nx.non_neighbors(G_new,path[2]):
                    G_new.remove_edge(path[2],path[3])
                    G_new.add_edge(path[0],path[2])
                    count_rewire = count_rewire +1
    return G_new 

def rewire_cluster_scale_free(G,rewire_iter,max_iter = 10000000):
    G_new = nx.Graph(G) 
    count_iter = 0
    count_rewire = 0
    while (count_iter < max_iter) & (count_rewire < rewire_iter):
        count_iter = count_iter + 1
        rand_node = random.choice(G_new.nodes())
        if len(list(nx.all_neighbors(G_new,rand_node))) >= 2:
            path = Random_Path(G_new, 2,Random_Path(G_new,2,[rand_node])[::-1])
            if len(path)== 5:
                if path[2] in nx.non_neighbors(G_new,path[4]):
                    G_new.remove_edge(path[1],path[2])
                    G_new.add_edge(path[2],path[4])
                    count_rewire = count_rewire +1
    return G_new 
    

    
def C_elegans_metabolic(m,n):
    return rewire_cluster_scale_free(rewire_scale_free(nx.gnm_random_graph(435,2026),n),m)
    
def Saccharomyces_metabolic():
    return rewire_cluster_scale_free(rewire_scale_free(nx.gnm_random_graph(1846,2203),100000),150)
    
def Batch_yeast(n,iter):
    Global_CC = []
    Local_CC = []
    Degree_Dists = []
    Modularitys = []
    for t in range(n):
        Gnew = Saccharomyces_metabolic()
        Global_CC.append(Global_cluster_coeff(Gnew))
        Local_CC.append(local_clust_cdf(Gnew))
        Degree_Dists.append([list(Gnew.degree().values()).count(k) for k in range(100)])
        Modularitys.append(community.modularity(community.best_partition(Gnew),Gnew))
        print(t)
    np.savetxt(".\Data\DegreeDist_Yeast_model_"+str(iter)+".csv",Degree_Dists,fmt='%d', delimiter=",")
    np.savetxt(".\Data\GlobalClustering_Yeast_model_"+str(iter)+".csv",Global_CC,fmt='%f', delimiter=",")
    np.savetxt(".\Data\LocalClusteringCDF_Yeast_model_"+str(iter)+".csv",Local_CC,fmt='%f', delimiter=",")
    np.savetxt(".\Data\Modunlarity_Yeast_model_"+str(iter)+".csv",Modularitys,fmt='%f', delimiter=",")

    
    nx.powerlaw_cluster_graph(453, 5, 0.6)    
    
def Batch_c_celegans(n,iter):
    Global_CC = []
    Local_CC = []
    Degree_Dists = []
    Modularitys = []
    for t in range(n):
        Gnew = C_elegans_metabolic()
        Global_CC.append(Global_cluster_coeff(Gnew))
        Local_CC.append(local_clust_cdf(Gnew))
        Degree_Dists.append([list(Gnew.degree().values()).count(k) for k in range(100)])
        Modularitys.append(community.modularity(community.best_partition(Gnew),Gnew))
        print(t)
    np.savetxt(".\Data\DegreeDist_worm_model_"+str(iter)+".csv",Degree_Dists,fmt='%d', delimiter=",")
    np.savetxt(".\Data\GlobalClustering_worm_model_"+str(iter)+".csv",Global_CC,fmt='%f', delimiter=",")
    np.savetxt(".\Data\LocalClusteringCDF_worm_model_"+str(iter)+".csv",Local_CC,fmt='%f', delimiter=",")
    np.savetxt(".\Data\Modunlarity_worm_model_"+str(iter)+".csv",Modularitys,fmt='%f', delimiter=",")