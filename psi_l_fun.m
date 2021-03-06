function [psi,l] = psi_l_fun(in1,r_c,P)
%PSI_L_FUN
%    [PSI,L] = PSI_L_FUN(IN1,R_C,P)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    31-Jul-2020 14:05:11

o = in1(:,5);
x = in1(:,1);
y = in1(:,2);
t2 = cos(o);
t3 = sin(o);
t6 = -x;
t4 = r_c.*t2;
t5 = r_c.*t3;
t7 = -t5;
t9 = P+t4+t6;
t10 = (t5-y).^2;
t8 = t7+y;
t11 = t9.^2;
t12 = t10+t11;
psi = o+acos(t9.*1.0./sqrt(t12));
if nargout > 1
    l = sqrt(t12);
end
