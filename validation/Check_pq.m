% check pq

t = out(end).result.solution.phase.time;
U = out(end).result.solution.phase.control;
p_ax = U(:,7:9);
q_ax = U(:,10:12);

p_rot = U(:,13:15);
q_rot = U(:,16:18);

subplot(2,1,1)
plot(t,p_ax)
resetcolor;hold on
plot(t,q_ax,'--')

subplot(2,1,2)
plot(t,p_rot)
resetcolor;hold on
plot(t,q_rot,'--')

clc
trapz(t,sum(p_ax + p_rot - q_ax - q_rot,2))

max(max(p_ax.*q_ax))

max(max(p_rot.*q_rot))