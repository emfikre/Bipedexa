function output = Bipedexa(auxdata)

%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
g= auxdata.g; %gravity
lmax=auxdata.lmax; %max leg length
d=auxdata.d; %step length
D=auxdata.D; %sride length
Fmax=auxdata.Fmax;
Taumax=auxdata.Taumax;
I=auxdata.I;
r=auxdata.r;
T=auxdata.T;
% specify auxdata if not already done

%-------------------------------------------------------------------%
%------------------------- Variable Bounds -------------------------%
%-------------------------------------------------------------------%
% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).initialtime.lower = 0;              % scalar
bounds.phase(i).initialtime.upper = 0;              % scalar
bounds.phase(i).finaltime.lower = T ;                % scalar
bounds.phase(i).finaltime.upper = T ;                % scalar
%States
%6 kinematic states
%3 Force states
%3 Torque states
%1 Integrated Force
%1 Integrated Torque
 xlow = 0;
 xupp = D;
 ylow = 0;
 yupp = Inf;
 Flow= zeros(1,3);
 Fupp= [1 1 1]*Fmax;
 Taulow= [1 1 1]*(-Taumax);
 Tauupp= [1 1 1]*Taumax;
bounds.phase(i).initialstate.lower = [xlow,ylow,0,0,0,0,Flow,Taulow,0,0];           % row vector, length = numstates
bounds.phase(i).initialstate.upper = [xlow,yupp,Inf,Inf,Inf,Inf,Fupp,Tauupp,0,Inf];           % row vector, length = numstates
bounds.phase(i).state.lower = [xlow,ylow,0,0,0,0,Flow,Taulow,0,0];             % row vector, length = numstates
bounds.phase(i).state.upper = [xupp,yupp,Inf,Inf,Inf,Inf,Fupp,Tauupp,Inf,Inf];                 % row vector, length = numstates
bounds.phase(i).finalstate.lower = [xupp,ylow,0,0,0,0,Flow,Taulow,0,0];             % row vector, length = numstates
bounds.phase(i).finalstate.upper = [xupp,yupp,Inf,Inf,Inf,Inf,Fupp,Tauupp,Fmax*T,Inf];             % row vector, length = numstates
% 3 Time derivative of force controls
% 3 time derivative of torque controls
neg=[1 1 1]*(-Inf);
pos= [1 1 1]*(Inf);
bounds.phase(i).control.lower = [neg,neg];                % row vector, length = numstates
bounds.phase(i).control.upper = [pos,pos];                % row vector, length = numstates
% ???
bounds.phase(i).integral.lower = 0;                 % row vector, length = numintegrals
bounds.phase(i).integral.upper = Inf;                 % row vector, length = numintegrals
% no parameters introduced
%bounds.parameter.lower = ;                      % row vector, length = numintegrals
%bounds.parameter.upper = ;                      % row vector, length = numintegrals

% Endpoint constraints (if required)

bounds.eventgroup.lower = zeros(1,2); % row vector
bounds.eventgroup.upper = zeros(1,2); % row vector

% Path constraints (if required)
% 6 complimentary limb length constaints 
% 2 complimentary exclusion constaints 

% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).path.lower = zeros(1,8); % row vector, length = number of path constraints in phase


bounds.phase(i).path.upper =[inf inf inf inf inf inf 0 0]; % row vector, length = number of path constraints in phase
%-------------------------------------------------------------------------%
%---------------------------- Provide Guess ------------------------------%
%-------------------------------------------------------------------------%
% ----- PHASE 1 ----- %
i = 1;
guess.phase(i).time    = [0;T];                % column vector, min length = 2
guess.phase(i).state   = rand(2,14);                % array, min numrows = 2, numcols = numstates
guess.phase(i).control = rand(2,6);               % array, min numrows = 2, numcols = numcontrols
guess.phase(i).integral = rand;               % scalar

%guess.parameter = [];                    % row vector, numrows = numparams


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
setup.mesh.maxiterations=2;

% not required

%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name                        = 'bipedalprob';
setup.functions.continuous        = @Continuous;
setup.functions.endpoint          = @Endpoint;
setup.auxdata                     = auxdata; % not necessary
setup.bounds                      = bounds;
setup.nlp.solver= 'snopt';
setup.guess                       = guess;

setup.derivatives.derivativelevel = 'first';


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
end


function phaseout = Continuous(input)

% extract data
t = input.phase(1).time;
X = input.phase(1).state;
U = input.phase(1).control;

auxdata = input.auxdata;
g= auxdata.g; %gravity
lmax=auxdata.lmax; %max leg length
d=auxdata.d; %step length
D=auxdata.D; %sride length
Fmax=auxdata.Fmax;
Taumax=auxdata.Taumax;
I=auxdata.I;
r=auxdata.r;
T=auxdata.T;
c1=auxdata.c1;
c2=auxdata.c2;

%P = input.phase(1).parameter;

x=X(:,1);
y=X(:,2);
theta=X(:,5);
xdot = X(:,3); % provide derivative
ydot = X(:,4);
thetadot=X(:,6);
F=X(:,7:9); %Collecting Forces 
Tau=X(:,10:12); %Collecting Torques
P= X(:,13);
Q= X(:,14);

Ftr=F(:,1);
Flead=F(:,2);
Fref=F(:,3);
Tautr=Tau(:,1);
Taulead=Tau(:,2);
Tauref=Tau(:,3);
Tautrsqr=Tautr.^2;
Tauleadsqr=Taulead.^2;
Taurefsqr=Tauref.^2;

Fdot=U(:,1:3); 
Taudot=U(:,4:6);

Pdot= Flead;
Qdot= Taulead.^2;

%ntime=length(x);
zs=zeros(size(x));
os=ones(size(x));
temp=[os,zs,zs];
dvec=temp*d;
Dvec=temp*D;

rvec=-r.*[cos(theta),sin(theta),zs];
xc=[x,y,zs];
ltr = xc +rvec;
llead=(xc+rvec)-dvec;
lref=(xc+rvec)-Dvec;

magnitudeltr=sqrt(dot(ltr,ltr,2));
magnitudellead=sqrt(dot(llead,llead,2));
magnitudelref= sqrt(dot(lref,lref,2));
ultr=ltr./magnitudeltr;
ullead=llead./magnitudellead;
ulref =lref./magnitudelref;

Ftrvec= Ftr.*ultr;
Fleadvec = Flead.*ullead;
Frefvec = Fref.*ulref;

crossFtr=cross(rvec,Ftrvec);
crossFtrz=crossFtr(:,3);  %Extracting z column 
crossFlead=cross(rvec,Fleadvec);
crossFleadz=crossFlead(:,3); %Extracting z column 
crossFref=cross(rvec,Frefvec);
crossFrefz=crossFref(:,3); %Extracting z column

xddotmat= Ftr.*(x./ltr)+Fref.*((x-dvec)./lref)+Flead.*((x-Dvec)./llead);
xddot=xddotmat(:,1);
yddotmat= Ftr.*(y./ltr)+Fref.*(y./lref)+Flead.*(y./llead)-g;
yddot=yddotmat(:,2);
thetaddot=Tautr+Taulead+Tauref+crossFtrz+crossFleadz+crossFrefz;

phaseout.dynamics = [xdot,ydot,xddot,yddot,thetadot,thetaddot,Fdot,Taudot,Pdot,Qdot];


%what is c1 and c2 respectively?
phaseout.integrand = c1*(Ftr.^2+Fref.^2+Flead.^2)+c2*(Tautr.^2+Tauref.^2+Taulead.^2);
%Path constraint
Ftrllc= Ftr.*(lmax-ltr);
magFtrllc= sqrt(dot(Ftrllc,Ftrllc,2));
Fleadllc= Flead.*(lmax-llead);
magFleadllc= sqrt(dot(Fleadllc,Fleadllc,2));
Frefllc= Fref.*(lmax-lref);
magFrefllc=sqrt(dot(Frefllc,Frefllc,2));
Tautrllc= Tautrsqr.*(lmax-ltr);
magTautrllc=sqrt(dot(Tautrllc,Tautrllc,2));
Tauleadllc=Tauleadsqr.*(lmax-ltr);
magTauleadllc=sqrt(dot(Tauleadllc,Tauleadllc,2));
Taurefllc= Taurefsqr.*(lmax-llead);
magTaurefllc=sqrt(dot(Taurefllc,Taurefllc,2));
Fxc=P.*Ftr;
Tauxc=Q.*Tautr;
phaseout.path = [magFtrllc,magFleadllc,magFrefllc,magTautrllc,magTauleadllc,magTaurefllc,Fxc,Tauxc]; % path constraints, matrix of size num collocation points X num path constraints
end

function output = Endpoint(input)
%Endpoint Constraint for Forces and Torques
Finalstates =input.phase(1).finalstate;
Initialstates= input.phase(1).initialstate;
Ftr=Initialstates(7); 
Flead=Finalstates(8);
Ttr=Initialstates(10);
Tlead=Finalstates(11);
output.eventgroup.event = [(Ftr-Flead) (Ttr-Tlead)];% event constraints (row vector)

J = input.phase(1).integral(1);
output.objective = J; % objective function (scalar)

end