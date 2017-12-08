Params_dm.nDims = 2;
DynModel = DynamicModelX(Params_dm);
 
DynModel.sys('wk',[0.5 1; 0.2 12],'xkm1',ans,'k',1);