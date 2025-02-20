%% 4 Cases
x01 = [0; 0.1; 0.1; 0];
x02 = [1; 0.1; 0.1; 0];
x03 = [0; 0.6; 0.1; 0];
x04 = [0; 0.1; 0.6; 1];
T = 0.1;
Period = 5;
sample_number = fix(Period/T);
linstep = linspace(0, T*sample_number, sample_number + 1);
y1 = zeros(sample_number + 1, 4); y2 = zeros(sample_number + 1, 4); y3 = zeros(sample_number + 1, 4); y4 = zeros(sample_number + 1, 4);
y1(1,:) = x01; y2(1,:) = x02; y3(1,:) = x03; y4(1,:) = x04;
K = [17.1674   14.7339   94.9922   25.8779];
for i = 1: sample_number
    t = [0 T];
    u1 = K*x01; u2 = K*x02; u3 = K*x03; u4 = K*x04;
    [t1,x1] = ode45(@(t,x)Nonlinear(t,x,u1),t,x01); [t2,x2] = ode45(@(t,x)Nonlinear(t,x,u2),t,x02); [t3,x3] = ode45(@(t,x)Nonlinear(t,x,u3),t,x03); [t4,x4] = ode45(@(t,x)Nonlinear(t,x,u4),t,x04);
    x01 = x1(end,:); x02 = x2(end,:); x03 = x3(end,:); x04 = x4(end,:);
    y1(i+1,:) = x1(end,:); y2(i+1,:) = x2(end,:); y3(i+1,:) = x3(end,:); y4(i+1,:) = x4(end,:);
    x01 = x01'; x02 = x02'; x03 = x03'; x04 = x04';
end
fig1 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
stairs(linstep,y1(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y2(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y3(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y4(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
grid on
title('Cart Displacement in 4 Cases')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
xlabel('Time(s)')
ylabel('Cart Displacement(m)')
subplot(1,2,2)
stairs(linstep,y1(:,3),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y2(:,3),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y3(:,3),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep,y4(:,3),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation in 4 Cases')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')
%% Different T
t0 = 0:0.01:5;
x0 = [0; 0.1; 0.1; 0];
[t,x1] = ode45(@NonlinearCT, t0, x0);
T2 = 0.08; T3 = 0.05; T4 = 0.01;
x012 = [0; 0.1; 0.1; 0]; x013 = [0; 0.1; 0.1; 0]; x014 = [0; 0.1; 0.1; 0];
sample_number2 = fix(Period/T2); sample_number3 = fix(Period/T3); sample_number4 = fix(Period/T4);
linstep2 = linspace(0, T2*sample_number2, sample_number2 + 1);
linstep3 = linspace(0, T3*sample_number3, sample_number3 + 1);
linstep4 = linspace(0, T4*sample_number4, sample_number4 + 1);
y12 = zeros(sample_number2 + 1, 4); y13 = zeros(sample_number3 + 1, 4); y14 = zeros(sample_number4 + 1, 4);
y12(1,:) = x012; y13(1,:) = x013; y14(1,:) = x014;
K = [17.1674   14.7339   94.9922   25.8779];
for i = 1: sample_number2
    t2 = [0 T2];
    u12 = K*x012;
    [t12,x12] = ode45(@(t,x)Nonlinear(t,x,u12),t2,x012);
    x012 = x12(end,:);
    y12(i+1,:) = x12(end,:);
    x012 = x012'; 
end
for i = 1: sample_number3
    t3 = [0 T3];
    u13 = K*x013;
    [t13,x13] = ode45(@(t,x)Nonlinear(t,x,u13),t3,x013);
    x013 = x13(end,:);
    y13(i+1,:) = x13(end,:);
    x013 = x013'; 
end
for i = 1: sample_number4
    t4 = [0 T4];
    u14 = K*x014;
    [t14,x14] = ode45(@(t,x)Nonlinear(t,x,u14),t4,x014);
    x014 = x14(end,:);
    y14(i+1,:) = x14(end,:);
    x014 = x014'; 
end
fig2 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
stairs(linstep,y1(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep2,y12(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep3,y13(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
stairs(linstep4,y14(:,1),'Linewidth',2,'Marker','o','MarkerSize',2)
hold on
plot(t, x1(:,1),'Linewidth',2)
grid on
title('Cart Displacement in different Sampling Time')
legend('T = 0.1', 'T = 0.08', 'T = 0.05', 'T = 0.01','Continuous')
xlabel('Time(s)')
ylabel('Cart Displacement(m)')
subplot(1,2,2)
stairs(linstep,y1(:,3),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(linstep2,y12(:,3),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(linstep3,y13(:,3),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
stairs(linstep4,y14(:,3),'Linewidth',2,'Marker','o','MarkerSize',1)
hold on
plot(t, x1(:,3),'Linewidth',2)
hold on
yline(0.758,'-.b','Upper Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(-0.758,'-.b','Lower Threshold for \phi(t)','LabelHorizontalAlignment','center','LabelVerticalAlignment','middle');
yline(0,'-.b','Equilibrium \phi(t)');
grid on
title('Pendulum Angular Rotation in different Sampling Time')
legend('T = 0.1', 'T = 0.08', 'T = 0.05', 'T = 0.01','Continuous')
xlabel('Time(s)')
ylabel('Pendulum Angular Rotation(rad)')






lsiminfo(y1(:,1),linstep,0)
lsiminfo(y1(:,3),linstep,0)
lsiminfo(y12(:,1),linstep2,0)
lsiminfo(y12(:,3),linstep2,0)
lsiminfo(y13(:,1),linstep3,0)
lsiminfo(y13(:,3),linstep3,0)
lsiminfo(y14(:,1),linstep4,0)
lsiminfo(y14(:,3),linstep4,0)



function dx = Nonlinear(t,x,u)
F = 1;
M = 1;
L = 0.842;
g = 9.8093;

dx = zeros(4,1);
dx(1) = x(2);
dx(2) = -F*x(2)/M + u/M;
dx(3) = x(4);
dx(4) = F*cos(x(3))*x(2)/(L*M) + g*sin(x(3))/L - cos(x(3))*u/(L*M);
end


% 
function dx = NonlinearCT(t,x)
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
