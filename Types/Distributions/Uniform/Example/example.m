% Example/Test script for UniformDistributionX

limits = [0, 10; ...
          0, 10]; 

unif_dist = UniformDistributionX(limits);

samples = unif_dist.random(1000);

like = unif_dist.pdf(randi([-5,5],2,50));