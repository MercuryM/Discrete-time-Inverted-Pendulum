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
A*B
A*A*B
A*A*A*B
R = [B A*B (A^2)*B (A^3)*B]
Tc = ctrb(A,B)
rankr = rank(Tc)
rank (R)
C =[1 0 0 0; 0 0 1 0]
C*A
C*A*A
C*A*A*A
O = [C; C*A; C*A^2; C*A^3]
To = obsv(A,C)
ranko = rank(To)
rank(O)
T = 0.01

Ad = exp(A*T)


% # inverted-pendulum.m
% #
% # Control system simulation for inverted pendulum. Based on Example 12-5 from
% # Modern Control Engineering, Fourth Edition by K. Ogata.

% pkg load control

% clear all
% close all
% 
% % # Inverted pendulum parameters
% g = 9.81;   %Gravity constant [m/s]
% M = 1.1;    %Cart mass [kg]
% m = 0.062;  %Pendulum mass [kg]
% l = 0.1425; %Pendulum length to center of mass [m]
% 
% % # Desired closed loop poles
% J = [-1+j*sqrt(3) -1-j*sqrt(3) -5 -5 -5];
% 
% % # State space representation of system
% A = [             0 1 0 0;
%       (M+m)/(M*l)*g 0 0 0;
%                   0 0 0 1;
%              -m/M*g 0 0 0];
% 
% B = [0 -1/(M*l) 0 1/M]';
% 
% C = [0 0 1 0];
% 
% D = [0];
% 
% Ahat = [ A zeros(4,1);
%         -C          0];
% 
% Bhat = [B; 0];
% 
% % Determine the state feedback gain matrix
% Khat = acker(Ahat, Bhat, J)
% K = Khat(1:4);
% kI = -Khat(5);
% 
% % # Calculate the step response
% ts = 0.02;
% t = 0:ts:6;
% 
% AA = [(A - B*K) B*kI;
%              -C    0];
% 
% BB = [0 0 0 0 1]';
% 
% CC = [C 0];
% 
% DD = [0];
% 
% sys = ss(AA, BB, CC, DD);
% [y, t, x] = step(sys, t);
% 
% u = (-Khat * x')';
% 
% % Plot the step response
% subplot(3,2,1)
% plot(t, x(:,1))
% grid
% title('\theta vs t')
% xlabel('t [s]')
% ylabel('\theta [rad]')
% 
% subplot(3,2,2)
% plot(t, x(:,2))
% grid
% title('\theta dot vs t')
% xlabel('t [s]')
% ylabel('\theta dot [rad/s]')
% 
% subplot(3,2,3)
% plot(t, x(:,3))
% grid
% title('x vs t')
% xlabel('t [s]')
% ylabel('x [m]')
% 
% subplot(3,2,4)
% plot(t, x(:,4))
% grid
% title('x dot vs t')
% xlabel('t [s]')
% ylabel('x dot [m/s]')
% 
% subplot(3,2,5)
% plot(t, x(:,5))
% grid
% title('\zeta vs t')
% xlabel('t [s]')
% ylabel('\zeta [m]')
% 
% subplot(3,2,6)
% plot(t, u)
% grid
% title('u vs t')
% xlabel('t [s]')
% ylabel('u [N]')
% 
% zeros(4,2)