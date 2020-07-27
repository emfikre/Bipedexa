function output = Bipedexa(auxdata,guess)

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
m=auxdata.m;

% specify auxdata if not already done
if ~isfield(auxdata,'scaling') 
    auxdata.scaling = 'none';
end

%-------------------------------------------------------------------%
%------------------------- Variable Bounds -------------------------%
%-------------------------------------------------------------------%
% ----- PHASE 1 ----- %

if strcmpi(auxdata.scaling,'none')
    bounds = getBounds('open',auxdata);
else
    bounds = getBounds('closed',auxdata); 
end

%-------------------------------------------------------------------------%
%---------------------------- Provide Guess ------------------------------%
%-------------------------------------------------------------------------%
% ----- PHASE 1 ----- %
% i = 1;
% guess.phase(i).time    = [0;T];                % column vector, min length = 2
% guess.phase(i).state   = rand(2,14);                % array, min numrows = 2, numcols = numstates
% guess.phase(i).control = rand(2,6);               % array, min numrows = 2, numcols = numcontrols
% guess.phase(i).integral = rand;               % scalar
if isempty(guess)
    i = 1;
    guess.phase(i).time    = [0;T];                % column vector, min length = 2
    guess.phase(i).state   = [0,lmax+r,D/T,0,pi/2,0,m/T/2,0,m/T/2,0,0,0,0,0; ...
        0,lmax+r,D/T,0,pi/2,0,0,m/T/2,m/T/2,0,0,0,m/T/2*(T/2),0];                % array, min numrows = 2, numcols = numstates
    guess.phase(i).control = [-m/T^2/2,m/T^2/2,0,0,0,0; ...
                              -m/T^2/2,m/T^2/2,0,0,0,0];               % array, min numrows = 2, numcols = numcontrols
    guess.phase(i).integral = [1,0.1*T];
    guess.parameter = zeros(1,8);

elseif strcmpi(guess,'rand')
    i = 1;
    guess = struct;
    guess.phase(i).time    = [0;T];                % column vector, min length = 2
    guess.phase(i).state   = rand(2,14);                % array, min numrows = 2, numcols = numstates
    guess.phase(i).control = rand(2,6);               % array, min numrows = 2, numcols = numcontrols
    guess.phase(i).integral = rand(1,2);               % scalar
    guess.parameter = rand(1,8);
elseif isstruct(guess)
     % it's an output struct from a previous trial
    guess1 = guess;
    guess = [];
    % pull out the guess
    guess.phase.time = guess1.result.solution.phase.time;
    guess.phase.state = guess1.result.solution.phase.state;
    guess.parameter = guess1.result.solution.parameter;
    guess.phase.control = guess1.result.solution.phase.control;
    guess.phase.integral = guess1.result.solution.phase.integral;
end
%guess.parameter = [];                    % row vector, numrows = numparams


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
setup = auxdata.setup; % can input any setup parameters into this field. May be overwritten!
setup.mesh.maxiterations= auxdata.meshiter;
setup.method= 'RPM-integration';

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
setup.nlp.snoptoptions.maxiterations = auxdata.snoptiter;
setup.guess                       = guess;
setup.scales.method               = auxdata.scaling;

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
Pa = input.phase(1).parameter;

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
c = auxdata.c;
m = auxdata.m;


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

% Collect force and torque rates of change
dF = U(:,1:3);
dTau = U(:,4:6);

% Collect relaxation parameters
sLimbF   = Pa(:,1:3);
sLimbTau = Pa(:,4:6);
sExclF   = Pa(:,7);
sExclTau = Pa(:,8);

Ftr=F(:,1);
Flead=F(:,2);
Fref=F(:,3);
Tautr=Tau(:,1);
Taulead=Tau(:,2);
Tauref=Tau(:,3);
Tautrsqr=Tautr.^2;
Tauleadsqr=Taulead.^2;
Taurefsqr=Tauref.^2;

Fdot=dF; 
Taudot=dTau;

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

Ftrvec = Ftr.*ultr;
Fleadvec = Flead.*ullead;
Frefvec = Fref.*ulref;

crossFtr=cross(rvec,Ftrvec);
crossFtrz=crossFtr(:,3);  %Extracting z column 
crossFlead=cross(rvec,Fleadvec);
crossFleadz=crossFlead(:,3); %Extracting z column 
crossFref=cross(rvec,Frefvec);
crossFrefz=crossFref(:,3); %Extracting z column

xddot= (Ftrvec(:,1)+Fleadvec(:,1)+Frefvec(:,1))/m;
yddot= (Ftrvec(:,2)+Fleadvec(:,2)+Frefvec(:,2))/m-g;
thetaddot=(Tautr+Taulead+Tauref+crossFtrz+crossFleadz+crossFrefz)/I;

phaseout.dynamics = [xdot,ydot,xddot,yddot,thetadot,thetaddot,Fdot,Taudot,Pdot,Qdot];


% the vector c contains scaling parameters, specifying the relative
% amplification of various terms

% First column: Main cost
% Second column: force and torque rate penalty
phaseout.integrand = ...
    [c(1)*(Ftr.^2+Fref.^2+Flead.^2)+c(2)*(Tautr.^2+Tauref.^2+Taulead.^2), ...
     c(3)*sum(dF.^2,2) + c(4)*sum(dTau.^2,2)]; 

%%% Path constraints

% Hips above ground; simply ltr(:,2).

% Force activation limb length constraints
trllc= (Ftr+Tautrsqr).*(lmax-magnitudeltr) - sLimbF(:,1);
leadllc= (Flead+Tauleadsqr).*(lmax-magnitudellead) - sLimbF(:,2);
refllc= (Fref+Taurefsqr).*(lmax-magnitudelref) - sLimbF(:,3);

% Torque activation limb length constraints
% Tautrllc= Tautrsqr.*(lmax-magnitudeltr) - sLimbTau(:,1);
% Tauleadllc=Tauleadsqr.*(lmax-magnitudellead) - sLimbTau(:,2);
% Taurefllc= Taurefsqr.*(lmax-magnitudelref) - sLimbTau(:,3);

% Limb exclusion constraints
Fxc=P.*Ftr - sExclF;
Tauxc=Q.*Tautr - sExclTau;

phaseout.path = [ltr(:,2),trllc,leadllc,refllc,Fxc,Tauxc]; % path constraints, matrix of size num collocation points X num path constraints
end

function output = Endpoint(input)

c = input.auxdata.c;
Pa = input.parameter;
sLimbF = Pa(1:3);
sLimbTau = Pa(4:6);
sExclF   = Pa(7);
sExclTau = Pa(8);

Finalstates =input.phase(1).finalstate;
Initialstates= input.phase(1).initialstate;

Ftr=Initialstates(7); 
Flead=Finalstates(8);

Ttr=Initialstates(10);
Tlead=Finalstates(11);


ybeg=Initialstates(2); %equal
yend=Finalstates(2);

xdotbeg=Initialstates(3);
xdotend=Finalstates(3);

ydotbeg=Initialstates(4);
ydotend=Finalstates(4);

thetabeg=Initialstates(5);
thetaend=Finalstates(5);

omegabeg=Initialstates(6);
omegaend=Finalstates(6);

output.eventgroup.event = [(Ftr-Flead) (Ttr-Tlead) (ybeg-yend) (xdotbeg-xdotend) (ydotbeg-ydotend) (thetabeg-thetaend) (omegabeg-omegaend)];% event constraints (row vector)

J1 = input.phase.integral(1); % F^2+Tau^2 cost
J2 = input.phase.integral(2); % Force rate penalty
J3 = c(5)*sum(sLimbF,2) + c(6)*sum(sLimbTau,2) + c(7)*sExclF + c(8)*sExclTau; % relaxation penalties
output.objective = J1+J2+J3; % objective function (scalar)

end

function bounds = getBounds(type,auxdata)

lmax = auxdata.lmax;
D = auxdata.D;
T = auxdata.T;
Fmax = auxdata.Fmax;
Taumax = auxdata.Taumax;
c = auxdata.c;
r = auxdata.r;

i = 1;
bounds = struct;
bounds.phase(i).initialtime.lower = 0;              % scalar
bounds.phase(i).initialtime.upper = 0;              % scalar
bounds.phase(i).finaltime.lower = T ;                % scalar
bounds.phase(i).finaltime.upper = T ;                % scalar

% Endpoint constraints (if required)

bounds.eventgroup.lower = [0,0,0,0,0,0,0]; % row vector
bounds.eventgroup.upper = [0,0,0,0,0,0,0]; % row vector

% Path constraints (if required)
% 1 normal (hip above ground)
% 6 complementarity limb length constaints 
% 2 complementarity exclusion constaints 

% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).path.lower = zeros(1,6); % row vector, length = number of path constraints in phase
bounds.phase(i).path.upper =[Inf inf inf inf 0 0]; % row vector, length = number of path constraints in phase

switch lower(type)
    case 'open'
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
        [ulow,vlow,wlow] = deal(-Inf);
        [uupp,vupp,wupp] = deal(Inf);
        thetalow = 0;
        thetaupp = pi;
        
        Flow= zeros(1,3);
        Fuppini=[1 0 0]*Fmax;
        Fuppfin=[0 1 0]*Fmax;
        Fupp= [1 1 1]*Fmax;
        Taulow= [1 1 1]*(-Taumax);
        Taulowini= [1 0 0]*(-Taumax);
        Tauuppini = [1 0 0]*(Taumax);
        Taulowfin =[0 1 0]*(-Taumax);
        Tauuppfin = [0 1 0]*(Taumax);
        Tauupp= [1 1 1]*Taumax;
        
        Plow = 0;
        Pupp = T*Fmax;
        Qlow = 0;
        Qupp = Taumax^2*T;
        
        % 3 Time derivative of force controls
% 3 time derivative of torque controls
dFneg=[1 1 1]*(-Inf);
dFpos= [1 1 1]*(Inf);
dTauneg = [1 1 1]*(-Inf);
dTaupos = [1 1 1]*(Inf);

% Integrals (cost)
[J1upp,J2upp] = deal(Inf);

% Parameters
% 3 relaxation parameters (leg length, Force)
% 3 relaxation parameters (leg length, torque)
% 1 relaxation paramter (leg exclusion, force)
% 1 relaxation paramter (leg exclusion, torque)
supp = Inf(1,8);                      % row vector, length = numintegrals

        
    case 'closed'
        %States
        %6 kinematic states
        %3 Force states
        %3 Torque states
        %1 Integrated Force
        %1 Integrated Torque
        xlow = 0;
        xupp = D;
        ylow = 0;
        yupp = (lmax+r)*4;
        [ulow,vlow] = deal(-4*D/T);
        [uupp,vupp] = deal( 4*D/T);
        
        thetalow = 0;
        thetaupp = pi;
        wlow = -4*pi/T;
        wupp =  4*pi/T;
        
        Flow= zeros(1,3);
        Fuppini=[1 0 0]*Fmax;
        Fuppfin=[0 1 0]*Fmax;
        Fupp= [1 1 1]*Fmax;
        Taulow= [1 1 1]*(-Taumax);
        Taulowini= [1 0 0]*(-Taumax);
        Tauuppini = [1 0 0]*(Taumax);
        Taulowfin =[0 1 0]*(-Taumax);
        Tauuppfin = [0 1 0]*(Taumax);
        Tauupp= [1 1 1]*Taumax;
        
        Plow = 0;
        Pupp = T*Fmax;
        Qlow = 0;
        Qupp = Taumax^2*T;
        
        % 3 Time derivative of force controls
        % 3 time derivative of torque controls
        dFneg=[1 1 1]*(-100*Fmax);
        dFpos= [1 1 1]*(100*Fmax);
        dTauneg = [1 1 1]*(-100*Taumax);
        dTaupos = [1 1 1]*(100*Taumax);
        
        % Integrals (cost)
        J1upp = c(1)*sum(Fmax.^2,2)*T + c(2)*sum(Taumax.^2,2)*T;
        J2upp = c(3)*sum(dFpos.^2,2)*T + c(4)*sum(dTaupos.^2,2)*T;
        
        % Parameters
        % 3 relaxation parameters (leg length, Force)
        % 3 relaxation parameters (leg length, torque)
        % 1 relaxation paramter (leg exclusion, force)
        % 1 relaxation paramter (leg exclusion, torque)
        supp = 1*ones(1,8);                      % row vector, length = numintegrals

        
end

bounds.phase(i).initialstate.lower =    [xlow,ylow,ulow,vlow,thetalow,wlow,Flow   ,Taulowini,Plow,Qlow];           % row vector, length = numstates
bounds.phase(i).initialstate.upper =    [xlow,yupp,uupp,vupp,thetaupp,wupp,Fuppini,Tauuppini,Pupp,Qupp];           % row vector, length = numstates
bounds.phase(i).state.lower =           [xlow,ylow,ulow,vlow,thetalow,wlow,Flow   ,Taulow   ,Plow,Qlow];             % row vector, length = numstates
bounds.phase(i).state.upper =           [xupp,yupp,uupp,vupp,thetaupp,wupp,Fupp   ,Tauupp   ,Pupp,Qupp];                 % row vector, length = numstates
bounds.phase(i).finalstate.lower =      [xupp,ylow,ulow,vlow,thetalow,wlow,Flow   ,Taulowfin,Plow,Qlow];             % row vector, length = numstates
bounds.phase(i).finalstate.upper =      [xupp,yupp,uupp,vupp,thetaupp,wupp,Fuppfin,Tauuppfin,Pupp,Qupp];             % row vector, length = numstates

bounds.phase(i).control.lower = [dFneg,dTauneg];                % row vector, length = numstates
bounds.phase(i).control.upper = [dFpos,dTaupos];                % row vector, length = numstates

bounds.phase(i).integral.lower = [0,0];                 % row vector, length = numintegrals
bounds.phase(i).integral.upper = [J1upp,J2upp];

bounds.parameter.lower = zeros(1,8);                      % row vector, length = numintegrals
bounds.parameter.upper = supp;

end