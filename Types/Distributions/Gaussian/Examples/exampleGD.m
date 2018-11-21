% Example usage of GaussianDistributionX

mu = rand(4,1);
P = diag(ones(1,4));

gd = GaussianDistributionX(mu,P);

mean = gd.Mean;
covar = gd.Covar;

samples = gd.random(10000);

lik = gd.pdf(samples);