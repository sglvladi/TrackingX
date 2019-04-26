

meas_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',2,'Mapping', [1,2], 'MeasurementErrVariance', 1);
H = meas_model.matrix();
R = meas_model.covar();
measurements = [10 15;10 15];

dist = GaussianDistributionX([1;1], 10*eye(2));

numParticles = 1000;
particles = dist.random(numParticles);
weights = repmat(1/numParticles,1,numParticles);

xPred = dist.Mean; 
PPred=dist.Covar; 
[yPred, S, K] = KalmanFilterX.predictMeasurement_(xPred,PPred,H,R);

lik = meas_model.pdf(measurements,particles);
l = meas_model.pdf(measurements,xPred,PPred);


% Data assoc
lambda = 1/10000;
Pd = 0.9;
hypothesiser = EfficientHypothesisManagementX();

W = hypothesiser.hypothesise([lambda*(1-Pd); sum(Pd*lik,2)]');
% W = lambda*(1-Pd);
% W = [W Pd*sum(l)];
% W = W./sum(W);
% 
% Wp = lambda*(1-Pd);
% Wp = [Wp Pd*sum(lik)];
% Wp = Wp./sum(Wp);
%W = [0.009 1-0.009];

%[xPost, PPost] = KalmanFilterX.update_(xPred,dist.Covar,measurements,yPred, S, K);
[xPost, PPost] = KalmanFilterX.updatePDA_(xPred,PPred,measurements,W,yPred, S, K);


[newWeights] = ParticleFilterX_UpdatePDA(@(y,x)meas_model.pdf(y,x),measurements,particles,weights,W,lik);
XPost = particles*newWeights';
resampler = SystematicResamplerX();
new_weights = weights.*lik;
%new_weights =  W(1)*dist.Weights + sum(W(2:end).*lik.*dist.Weights,1);
new_weights = new_weights./sum(new_weights);
[new_particles, new_weights] = resampler.resample(particles,new_weights);

new_dist = ParticleDistributionX(new_particles, new_weights);

% [X,Y] = meshgrid(-3:0.1:5,0:0.1:8);
% Z = dist.pdf([X;Y]);
% for i=1:size(X,1)
%     Z(i,:) = dist.pdf([X;Y]);
% end

figure;
hold on;
[bandwidth,Z,X,Y]=kde2d([particles, measurements+meas_model.random(numParticles)]');
%contour3(X,Y,density,50);
h = surf(X,Y,Z);  
shading interp
colormap(jet(3000))
[bandwidth,Z,X,Y]=kde2d([new_dist.Particles, measurements+meas_model.random(numParticles)]');
%contour3(X,Y,density,50);
h = surf(X,Y,Z);  
shading interp
colormap(jet(3000))
plot(measurements(1,:),measurements(2,:),'r*');
plot_gaussian_ellipsoid(dist.Particles*dist.Weights',weightedcov(dist.Particles,dist.Weights));
plot_gaussian_ellipsoid(new_dist.Particles*new_dist.Weights',weightedcov(new_dist.Particles,new_dist.Weights));

