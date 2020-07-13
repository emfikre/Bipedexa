% bipedexa MAIN

auxdata.g = 1;
auxdata.lmax = 1;
auxdata.T = 1.5;
auxdata.D = 0.5;
auxdata.m = 1;
auxdata.d = auxdata.D/2;
auxdata.Fmax = 4*auxdata.g;
auxdata.Taumax = 4*auxdata.g*auxdata.lmax;
auxdata.r = 0.5*auxdata.lmax;
auxdata.I = auxdata.m*auxdata.g*auxdata.r^2;
auxdata.c = [10,1000,0.1,0.1,1,1];
guess=[];
out = Bipedexa(auxdata,guess);

%%
close all
plotStates(out)

%% 
close all
Bipedexa_animate(out,'test')