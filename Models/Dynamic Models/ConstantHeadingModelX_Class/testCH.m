config.q_vel = 0.01;
config.q_head = 0.16;

CHmodel = ConstantHeadingModelX(config);

CHmodel2 = CHmodel;

xkm1 = [0; 0; sqrt(2); pi/4];
pkm1 = [0 1; 0 1; sqrt(2) sqrt(2); pi/4 pi/4];
Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
xk = CHmodel.sys(1, xkm1,2);
pk = CHmodel.sys(1, pkm1,2);
Pk = CHmodel.sys_cov(1, Pkm1,2);