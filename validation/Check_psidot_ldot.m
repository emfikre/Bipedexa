X = out(end).result.interpsolution.phase.state;
t = out(end).result.interpsolution.phase.time;
aux = out(end).result.setup.auxdata;
d = aux.d;
D = aux.D;
r = aux.r;
lmax = aux.lmax;
i = 1;

P = [0,D,d];

F = X(:,7:9);
Tau = X(:,10:12);

[psidot,ldot] = psidot_ldot_fun(X(:,1:6),r,P(i));
[psi,l] = psi_l_fun(X(:,1:6),r,P(i));

psidot_num = [diff(psi)./diff(t);NaN];
ldot_num = [diff(l)./diff(t);NaN];

close all

subplot(2,1,1)
plot(t,ldot)
hold on
resetcolor
plot(t,ldot_num,'o')

yyaxis right
plot(t,l)

subplot(2,1,2)
plot(t,psidot)
hold on
resetcolor
plot(t,psidot_num,'o')

yyaxis right
plot(t,psi)


figure
plot(t,lmax - l)
hold on
plot(t,F(:,i)./max(F(:,i)))
plot(t,Tau(:,i))
