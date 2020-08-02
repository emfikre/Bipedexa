% BipedWork MAIN
close all
auxdata.g = 1;
auxdata.lmax = 1;
auxdata.T = 2;
auxdata.D = 1;
auxdata.m = 1;
auxdata.d = auxdata.D/2;
% The following parameters roughly match Hasaneini et al. 2013, table 1
auxdata.Fmax = 4*auxdata.m*auxdata.g; 
auxdata.Taumax = 1*auxdata.m*auxdata.g*auxdata.lmax;
auxdata.r = 0.4*auxdata.lmax;
auxdata.I = 0.09*auxdata.m*auxdata.lmax^2;

dFpen   = [1 0.1 0.1]*1e-1;
dTaupen = [1 0.1 0.1]*1e-1;

auxdata.scaling = 'none';

guess='rand';
auxdata.setup.mesh.tolerance = 1e-3;
auxdata.snoptiter = 500;
auxdata.meshiter = 2;
auxdata.c = [1,1,dFpen(1),dTaupen(1),0.1,0.1,0.1,1000,0.1,0.1];
out = BipedWork(auxdata,guess);
plotStates(out)
%%
auxdata.scaling = 'automatic-hybrid';

guess2 = out(1);
auxdata.setup.mesh.tolerance = 1e-4;
auxdata.snoptiter = 1000;
auxdata.meshiter = 3;
auxdata.c = [1,1,dFpen(2),dTaupen(2),10,10,10,1000,10,10];
out(2) = BipedWork(auxdata,guess2);
plotStates(out(2))
%%
guess3 = out(2);
auxdata.setup.mesh.tolerance = 1e-4;
auxdata.setup.nlp.snoptoptions.tolerance = 1e-7;
auxdata.snoptiter = 1500;
auxdata.meshiter = 4;
auxdata.c = [1,1,dFpen(3),dTaupen(3),100,100,100,10000,100,100];
out(3) = BipedWork(auxdata,guess3);
plotStates(out(3))
%%
if out(3).result.maxerror > out(3).result.setup.mesh.tolerance || out(3).result.nlpinfo > 10
    % hasn't converged; try again
   guess4 = out(3);
   auxdata.meshiter = 5;
   auxdata.snoptiter = 2000;
   out(4) = BipedWork(auxdata,guess4); 
   plotStates(out(4))
end
fname = [date_prefix('yyyymmddHHMM'),'_sim'];
try
export_fig([fname,'.pdf']) % a function I like for plotting stuff. You probably don't have it.; requires some setup
catch
end
%% Animate the solution

Biped_animate(out(end),fname)
%% Save the workspace
save(fname,'out*','auxdata')
