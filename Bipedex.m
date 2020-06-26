function output = Bipedex(g,l,r)
% Data for the problem

% Variable Bounds
% States include x,y,xdot,ydot,theta,omega,F1,F2,Torque 1
i=1
bounds.phase(i).initialtime.lower= 0
bounds.phase(i).initialtime.upper= 0
bounds.phase.finaltime.lower= 0.05
bounds.phase.finaltime.upper= T
bounds.phase.initialstate.lower = [0,0,0,0,0,0]; %% xdot,ydot,theta, and omega are unbounded
bounds.phase.initialstate.upper = [0,0,0,0,0,0]; 
bounds.phase.state.lower = [-inf,-inf,0,0,0,0]; 
bounds.phase.state.upper = [inf,inf,F1upp,F2upp,tor1upp,tor2upp]; 
bounds.phase.finalstate.lower = [D,l,0,0,0,0]; 
bounds.phase.finalstate.upper = [D,l,F1upp,F2upp,tor1upp,tor2upp]; 
bounds.phase.control.lower = Fdotmin; 
bounds.phase.control.upper = Fdotmax;
end

