%% system
F = 1;
M = 1;
L = 0.842;
g = 9.8093;

A = [             0 1 0 0;
               0 -F/M 0 0;
                  0 0 0 1;
             0 F/(M*L) g/L 0];

B = [0 1/M 0 -1/(L*M)]';
C = [1 0 0 0; 0 0 1 0];
D = 0;
G = ss(A, B, C, 0);
P = [-1-sqrt(3)*j -1+sqrt(3)*j -5 -10]

%% Discretization
Ts = 0.1;
Gd = c2d(G,Ts, 'zoh')
[Ad,Bd,Cd,Dd,TS] = ssdata(Gd)
Pd = exp(P*Ts)
Kd = -1 * acker(Ad, Bd, Pd)
Ad1 = Ad + Bd * Kd
Cd1 = [1 0 0 0; 0 0 1 0]
Gd1 = ss(Ad1, Bd, Cd1, 0, Ts)
figure(1)
[Z P K] = ss2zp(Ad1, Bd, Cd1,[0;0]);
zplane(Z,P);
figure(2)
step(Gd1)
%% Plot 4 Cases
t = 0: 0.01: 5;
x01 = [0; 0.1; 0.1; 0]
[yd1, t, xd1]=initial(Gd1,x01,t);
x02 = [1; 0.1; 0.1; 0]
[yd2, t, xd2] = initial(Gd1,x02,t);
x03 = [0; 0.6; 0.1; 0]
[yd3, t, xd3] = initial(Gd1,x03,t);
x04 = [0; 0.1; 0.6; 1]
[yd4, t, xd4] = initial(Gd1,x04,t);
fig2 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
stairs(t,yd1(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd2(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd3(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd4(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
grid on
title('Cart Displacement in 4 Cases')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
xlabel('Time(s)')
ylabel('Cart Displacement(m)')
subplot(1,2,2)
stairs(t,yd1(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd2(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd3(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd4(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation in 4 Cases')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')


% [y, t, x] = initial(Gd1,x0);
% u = (Kd * x')';
% 
% subplot(2,1,1)
% plot(t, y(:,1))
% grid on
% title('s vs t')
% xlabel('t [s]')
% ylabel('s [m]')
% 
% subplot(2,1,2)
% plot(t, y(:,2))
% grid on
% title('\phi vs t')
% xlabel('t [s]')
% ylabel('phi [rad/s]')
% 
