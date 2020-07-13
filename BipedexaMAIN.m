% bipedexa MAIN
close all
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

guess=[];
auxdata.snoptiter = 500;
auxdata.meshiter = 2;
auxdata.c = [10,1000,0.01,0.01,0.01,0.01];
out1 = Bipedexa(auxdata,guess);
plotStates(out1)

guess2 = out1;
auxdata.snoptiter = 1000;
auxdata.meshiter = 2;
auxdata.c = [10,1000,1,1,0.1,0.1];
out2 = Bipedexa(auxdata,guess2);
plotStates(out2)

guess3 = out2;
auxdata.snoptiter = 1500;
auxdata.meshiter = 4;
auxdata.c = [10,1000,10,10,1,1];
out3 = Bipedexa(auxdata,guess3);
plotStates(out3)


%% 
close all
Bipedexa_animate(out3,'test')