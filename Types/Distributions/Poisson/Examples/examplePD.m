% Example usage of PoissonDistributionX

lambda = 5;

pd = PoissonDistributionX(lambda);

mu = pd.Mean;
covar = pd.Covar;
samples = pd.random(5000000);

lik = pd.pdf(samples);