import numpy as np
from sklearn.cluster import AgglomerativeClustering

D = np.array([[0, 1, 2, 4, 4], [1, 0, 2, 4, 4], [2, 2, 0, 2, 2], [4, 4, 2, 0, 1], [4, 4, 2, 1, 0]])
model = AgglomerativeClustering(linkage='average', metric='precomputed', distance_threshold = None)
#cluster the objects based on the distances
cluster = model.fit(D)
n_objects = len(cluster.labels_) # this is just the number of sequences
# cluster.children_ gives the hierarchical clustering
# leaves are the original objects and labelled from 0 to 4 here
# the internal nodes in the cluster are labelled from 5 onwards
next_node = n_objects 
for i, merge in enumerate(cluster.children_):
    print("Align ", merge[0], " with ", merge[1], " to give ", next_node)
    next_node+=1
print("done")