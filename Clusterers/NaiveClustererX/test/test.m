validationMatrix = [1     1     0     0     0
                    1     0     1     0     0
                    0     0     0     0     0
                    0     0     0     1     0
                    0     0     0     0     1
                    0     0     0     0     0
                    0     0     0     0     1
                    0     0     0     0     0];
                
clusterer = NaiveClustererX();

clusters_test = clusterer.cluster(validationMatrix);

load('eval.mat');
assert(isequaln(clusters_test,clusters_eval));