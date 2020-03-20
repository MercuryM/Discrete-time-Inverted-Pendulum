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
G = ss(A, B, C, 0)

%% Discretization
Ts = 0.01;
Gd = c2d(G,Ts, 'zoh')
[Ad,Bd,Cd,Dd,TS] = ssdata(Gd)

%% LQR
Q1 = [100,0,0,0;
       0,1,0,0;
     0,0,100,0;
       0,0,0,1];
R = 1;
[Kdo1,S1,e1] = dlqr(Ad,Bd,Q1,R)
Ado1 = Ad - Bd * Kdo1;
Cdo1 = [1 0 0 0; 0 0 1 0; -Kdo1]
Gdo1 = ss(Ado1, Bd, Cdo1, 0, Ts)
step(Gdo1)
%% Given initial state & Plot the result
t = 0: 0.01: 5;
x01 = [0; 0.1; 0.1; 0]
[yd1, t, xd1]=initial(Gdo1,x01,t);
x02 = [1; 0.1; 0.1; 0]
[yd2, t, xd2] = initial(Gdo1,x02,t);
x03 = [0; 0.6; 0.1; 0]
[yd3, t, xd3] = initial(Gdo1,x03,t);
x04 = [0; 0.1; 0.6; 1]
[yd4, t, xd4] = initial(Gdo1,x04,t);
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
%% Compare different Q
Q2 = [1000,0,0,0;
       0,1,0,0;
     0,0,100,0;
       0,0,0,1];
R = 1;
[Kdo2,S2,e2] = dlqr(Ad,Bd,Q2,R)
Ado2 = Ad - Bd * Kdo2;
Cdo2 = [1 0 0 0; 0 0 1 0; -Kdo2]
Gdo2 = ss(Ado2, Bd, Cdo2, 0, Ts)
Q3 = [100,0,0,0;
       0,1,0,0;
     0,0,1000,0;
       0,0,0,1];
R = 1;
[Kdo3,S3,e3] = dlqr(Ad,Bd,Q3,R)
Ado3 = Ad - Bd * Kdo3;
Cdo3 = [1 0 0 0; 0 0 1 0; -Kdo3]
Gdo3 = ss(Ado3, Bd, Cdo3, 0, Ts)
Q4 = [1000,0,0,0;
       0,1,0,0;
     0,0,1000,0;
       0,0,0,1];
R = 1;
[Kdo4,S4,e4] = dlqr(Ad,Bd,Q4,R)
Ado4 = Ad - Bd * Kdo4;
Cdo4 = [1 0 0 0; 0 0 1 0; -Kdo4]
Gdo4 = ss(Ado4, Bd, Cdo4, 0, Ts)
t = 0: 0.01: 5;
x01 = [0; 0.1; 0.1; 0]
[yd1, t, xd1]=initial(Gdo1,x01,t);
[yd12, t, xd12]=initial(Gdo2,x01,t);
[yd13, t, xd13]=initial(Gdo3,x01,t);
[yd14, t, xd14]=initial(Gdo4,x01,t);
fig3 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
stairs(t,yd1(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd12(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd13(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(t,yd14(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
grid on
title('Cart Displacement for different Q')
legend('Q11=100, Q33=100', 'Q11=1000, Q33=100', 'Q11=100, Q33=1000', 'Q11=1000, Q33=1000')
xlabel('Time(s)')
ylabel('Cart Displacement(m)')
subplot(1,2,2)
stairs(t,yd1(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd12(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd13(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(t,yd14(:,2),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation for different Q')
legend('Q11=100, Q33=100', 'Q11=100, Q33=1000', 'Q11=1000, Q33=100', 'Q11=1000, Q33=1000')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')



