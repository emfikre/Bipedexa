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
Frefvec =1.25* ulref;
Tautr=Tau(:,1);
Taulead=Tau(:,2);
Tauref=Tau(:,3);

Taurefrfmag= Tauref./magnitudelref; %calculating magnitude of ground reaction force produced by torque
legangle= acos(ultr(:,1));
torangle= legangle+(pi/2);
j=1:length(t);
utor=zeros(length(t),2);
utor(j,1)=cos(torangle);
utor(j,2)=sin(torangle);


Taureff=utor.*Taurefrfmag;


eqforceref=[(Taureff(:,1)+Frefvec(:,1)),(Taureff(:,2)+Frefvec(:,2))];




% extracted do the dynamics here.

%quiver function THIS IS WHERE YOU CAN IMPLEMNT THE COM FRAMED SYSTEM, THE
%com is always 0,0 so subtect x,y and you will getinitial posistion of the
i = 1:15:length(t);
Px=d-x;
Py=0-y;

close all
axis equal

figure('color','w')
quiver(Px(i,1),Py(i,1),eqforceref(i,1),eqforceref(i,2),0)
hold on
grid on
plot(0,0,'o')
hold on 
plot(0,0.145,'or')
% title('GRFs during step of one leg relative to the COM')
xlabel('x/l_b')
ylabel('y/l_b')
errorbar(0,0.145,0.035,0.035)
legend({'GRFs','COM',['Empirical VPP', newline, 'Maus 2010'] ,'Range'})
end

