%% Continuous Non-linear System
t0 = 0:0.01:5;
x01 = [0;0.1;0.1;0];
x02 = [1;0.1;0.1;0];
x03 = [0;0.6;0.1;0];
x04 = [0;0.1;0.6;1];
[t,x1] = ode45(@Nonlinear, t0, x01)
[t,x2] = ode45(@Nonlinear, t0, x02)
[t,x3] = ode45(@Nonlinear, t0, x03)
[t,x4] = ode45(@Nonlinear, t0, x04)
fig1 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
plot(t, x1(:,1),t, x2(:,1),t, x3(:,1),t, x4(:,1),'LineWidth',2)
grid on
title('Cart Displacement in 4 Cases')
xlabel('Time(s)')
ylabel('Cart Displacement(m)')
legend('Case 1','Case 2','Case 3','Case 4')
subplot(1,2,2)
plot(t, x1(:,3),t, x2(:,3),t, x3(:,3),t, x4(:,3),'LineWidth',2)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation in 4 Cases')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')
legend('Case 1','Case 2','Case 3','Case 4')

% Compare with linear system
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
K = [29.3864 26.8975 121.5071  34.5548];
G = ss(A + B*K, B, C, 0)
[yl1, t, xl1] = initial(G,x01,t0);
[yl4, t, xl4] = initial(G,x04,t0);
fig2 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
plot(t, x1(:,1),t, xl1(:,1),'--',t, x4(:,1),t, xl4(:,1),'--','LineWidth',2)
hold on 
grid on
title('Cart Displacement in Case 1 and 4')
xlabel('Time(s)')
ylabel('Cart Displacement(m)') 
legend('Nonlinear system(Case 1)','Linear system(Case 1)','Nonlinear system(Case 4)','Linear system(Case 4)')
subplot(1,2,2)
plot(t, x1(:,3),t, xl1(:,3),'--',t, x4(:,3),t, xl4(:,3),'--','LineWidth',2)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation in Case 1 and 4')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')
legend('Nonlinear system(Case 1)','Linear system(Case 1)','Nonlinear system(Case 4)','Linear system(Case 4)')

figure(3)
u1 = K*x1';
u4 = K*x4';
ul1 = K*xl1';
ul4 = K*xl4';
plot(t, u1,t, ul1,'--',t, u4,t, ul4,'--','LineWidth',2)
grid on
title('State Feedback in Case 1 and 4')
xlabel('Time(s)')
ylabel('State Feedback(N)')
legend('Nonlinear system(Case 1)','Linear system(Case 1)','Nonlinear system(Case 4)','Linear system(Case 4)')



% figure(5)
% subplot(4,1,1)
% plot(t, xc1(:,1))
% grid on
% title('s vs t')
% xlabel('t [s]')
% ylabel('s [m]')
% 
% subplot(4,1,2)
% plot(t, xc1(:,2))
% grid on
% title('v vs t')
% xlabel('t [s]')
% ylabel('v [m/s]')
% 
% subplot(4,1,3)
% plot(t, xc1(:,3))
% grid on
% title('\phi vs t')
% xlabel('t [s]')
% ylabel('\phi [rad]')

%% Discrete Non-linear System
% Fs = 1000,Ts = 0.001
% Ts = 0.001;
% Fs = 1/Ts;
% n = 0:10000;
% nTs = n*Ts;
% [t1,x1] = ode45(@Nonlinear, nTs, x0)
%  figure(2)
% subplot(3,1,1)
% plot(t1, x1(:,1))
% grid on
% hold on
% % stem(n*Ts,x1(:,1))
% title('s vs t')
% xlabel('t [s]')
% ylabel('s [m]')
% 
% subplot(3,1,2)
% plot(t1, x1(:,2))
% grid on
% title('v vs t')
% xlabel('t [s]')
% ylabel('v [m/s]')
% 
% 
% subplot(3,1,3)
% plot(t1, x1(:,3))
% grid on
% title('\phi vs t')
% xlabel('t [s]')
% ylabel('\phi [rad]')
% 


% stairs(ScopeData(:,1),ScopeData(:,2))
function dx = Nonlinear(t,x)
F = 1;
M = 1;
L = 0.842;
g = 9.8093;
K = [17.1674   14.7339   94.9922   25.8779];
u = K * x;
dx = zeros(4,1);
dx(1) = x(2);
dx(2) = -F*x(2)/M + u/M;
dx(3) = x(4);
dx(4) = F*cos(x(3))*x(2)/(L*M) + g*sin(x(3))/L - cos(x(3))*u/(L*M);
end
