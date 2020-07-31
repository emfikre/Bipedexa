clear

syms x y u v o w r_c P F

xp = x-P;
co = cos(o);
so = sin(o);
W = x - r_c.*co - P;
H = y - r_c.*so;
l = sqrt(W.^2 + H.^2);
dW = u + r_c.*w.*so;
dH = v - r_c.*w.*co;

phidot = w - (W.*dH - H.*dW)./l.^2;
phidot = simplify(phidot);

phidot_paper = w + ( y.*u - xp.*v - r_c.^2.*w + (xp.*w + v).*r_c.*co -(-y.*w +u).*r_c.*so )./(xp.^2 + y.^2 +r_c.^2 - 2.*r_c.*(xp.*co+y.*so));


simplify(( y.*u - xp.*v - r_c.^2.*w + (xp.*w + v).*r_c.*co -(-y.*w +u).*r_c.*so ) -  ( -(W.*dH - H.*dW) ) ) % should be zero
simplify((xp.^2 + y.^2 +r_c.^2 - 2.*r_c.*(xp.*co+y.*so)) - l.^2) % should be zero


%% compute ldot

ldot = (W.*dW + H.*dH)./l;

ldot_paper = (xp.*u + y.*v - (u+y.*w).*r_c.*co + (-v + xp.*w).*r_c.*so)./l;

simplify(ldot-ldot_paper) % should be zero


%% create function to calculate ldot and phidot at once
% using one function instead of two saves some computational time

matlabFunction([phidot,ldot],'File','phidot_ldot_fun','vars',{[x y u v o w],r_c,P});
