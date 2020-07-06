function plotGRF(GPOPSoutput)

t = GPOPSoutput.result.interpsolution.phase.time;
F = GPOPSoutput.result.interpsolution.phase.state(:,7:9);

ti = GPOPSoutput.result.solution.phase.time; % values at collocation points
Fi = GPOPSoutput.result.solution.phase.state(:,7:9);


figure('color','w')

plot(t,F)
hold on
resetcolor
plot(ti,Fi,'o')