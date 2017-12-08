config.xDim = 4;
config.yDim = 2;
config.r = 0.1;

model = PositionalObsModelX(config);

xkm1 = [0.5; 0.3; sqrt(2); pi/4];
pkm1 = [0 1; 0 1; sqrt(2) sqrt(2); pi/4 pi/4];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
yk = model.obs(1,xkm1);
pk = model.obs_cov();
Pk = model.eval_likelihood(1,yk,xkm1);