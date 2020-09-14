function plotVPP(GPOPSoutput) %Force vectors would probs go here 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
aux = GPOPSoutput.result.setup.auxdata;
D = aux.D;
T = aux.T;
g = aux.g;
c = aux.c;
r= aux.r;
d= aux.d;

t = GPOPSoutput.result.interpsolution.phase.time;
X = GPOPSoutput.result.interpsolution.phase.state; 

x= X(:,1);
y=X(:,2);
theta=X(:,5);
F = X(:,7:9);
Ftr=F(:,1);
Flead=F(:,2);
Fref=F(:,3);
Tau= X(:,10:12);

zs=zeros(size(x));
os=ones(size(x));
temp=[os,zs,zs];
dvec=temp*d;
Dvec=temp*D;

rvec=-r.*[cos(theta),sin(theta),zs];
xc=[x,y,zs];
ltr = xc +rvec;
llead=(xc+rvec)-Dvec;
lref=(xc+rvec)-dvec;

magnitudeltr=sqrt(dot(ltr,ltr,2));
magnitudellead=sqrt(dot(llead,llead,2));
magnitudelref= sqrt(dot(lref,lref,2));
ultr=ltr./magnitudeltr;
ullead=llead./magnitudellead;
ulref =lref./magnitudelref;

Ftrvec = Ftr.*ultr;
Fleadvec = Flead.*ullead;
Frefvec = ulref*1.25;
Tautr=Tau(:,1);
Taulead=Tau(:,2);
Tauref=Tau(:,3);

Taurefrfmag= Tauref./magnitudelref; %calculating magnitude of ground reaction force produced by torque
uTauref=[pi-(cos((pi/2)+acos(x./magnitudelref))),pi-(sin((pi/2)+acos(x./magnitudelref))),zs]; %calculating direction at which vector will point

Taureff=uTauref;


eqforceref=[(Taureff(:,1)+Frefvec(:,1)),(Taureff(:,2)+Frefvec(:,2))];




% extracted do the dynamics here.

%quiver function THIS IS WHERE YOU CAN IMPLEMNT THE COM FRAMED SYSTEM, THE
%com is always 0,0 so subtect x,y and you will get initial posistion of the

i = 1:10:length(t);

Px=d-x;
Py=0-y;
figure 
close all
axis equal
quiver(Px(i),Py(i),Frefvec(i,1),Frefvec(i,2),0)

end

