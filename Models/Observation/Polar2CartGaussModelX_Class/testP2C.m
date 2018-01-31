config.xDim = 4;
config.yDim = 2;
config.R = diag([ 2*(pi/180); 10 ]);

model = PositionalObsModelX(config);%PositionalObsModelX(config);%Polar2CartGaussianModelX(config);

xkm1 = [0.5; 0.3; sqrt(2); sqrt(2)];
pkm1 = [0 0.5; 0 0.3; sqrt(2) sqrt(2); sqrt(2) sqrt(2)];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
yk = model.obs(1,xkm1);
Yk = model.sample(1, yk, 500000);
pk = model.obs_cov();
Pk = model.eval(1,Yk,xkm1);
[bandwidth,density,X,Y]=kde2d(Yk');
figure;
%contour3(X,Y,density,50);
surf(X,Y,density);
hold on;