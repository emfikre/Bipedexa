% bipedexa MAIN

auxdata.g = 10;
auxdata.m=65;
auxdata.lmax = 0.8;
auxdata.T = 0.5;
auxdata.D = 0.25;
auxdata.d = auxdata.D/2;
auxdata.Fmax = 4*auxdata.g;
auxdata.Taumax = 4*auxdata.g*auxdata.lmax;
auxdata.r = 0.2*auxdata.lmax;
auxdata.I = auxdata.g*auxdata.r^2;
auxdata.c1=10;
auxdata.c2=10;
guess=[];
out = Bipedexa(auxdata,guess);

%%
close all
plotStates(out)