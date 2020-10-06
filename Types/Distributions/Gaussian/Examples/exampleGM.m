% Example usage of GaussianMixtureX

% Initialise some components
mu  = randi(500,2,8);
p   = [10^2 0; 0 10^2];
P = repmat(p,1,1,8);
w   = rand(1,8);

% Test object construction from means, covars and weights
gm = GaussianMixtureX(mu,P,w);

% Create component array structure 
components = cell(0,8);
for i = 1:8
    components{i}.Mean = mu(:,i);
    components{i}.Covar = P(:,:,i);
    components{i}.Weight = w(:,i);
end
% Test object construction from structure
gm = GaussianMixtureX(components);

% Extract information
n = gm.NumComponents;
m = gm.NumVariables;
means = gm.Means;
covars = gm.Covars;
weights = gm.Weights;

% Sample some data
samples = gm.random(10000);

% Evaluare likelihood
prob = gm.pdf(samples);

% Cluster the samples
[a,b,c,d,e] = gm.cluster(samples');


plot(samples(1,:), samples(2,:), '.')
[bandwidth,density,X,Y]=kde2d(samples');
h = surf(X,Y,density);
shading interp
colormap(jet(3000))
hold on
plotgaussellipse(gm.Mean,gm.Covar)