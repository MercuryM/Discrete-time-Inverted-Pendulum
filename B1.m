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

figure(1)
step(G)
grid on

%% state feedback control
P1 = [-1-sqrt(3)*j -1+sqrt(3)*j -5 -10]
P2 = [-7-7*j -7+7*j -35 -50]

K1 = -1 * place(A,B,P1)
A1 = A + B * K1
G1 = ss(A1, B, C, 0)
figure(2)
step (G1)
grid on

K2 = -1 * place(A,B,P2)
A2 = A + B *K2
G2 = ss(A2, B, C, 0)
figure(3)
step(G2)
grid on

ts = 0.01;
t = 0 : ts : 5;
x0 = [0; 0.1; 0.1; 0]
[y1,t1,x1] = initial(G1,x0,t)
[y2,t2,x2] = initial(G2,x0,t)
figure(4)
yyaxis left;
plot(t1, x1(:,2),'LineWidth',2)
hold on
plot(t2, x2(:,2),'LineWidth',2)
ylabel('Cart Velocity(m/s)')
yyaxis right;
plot(t1, x1(:,4),'LineWidth',2)
hold on
plot(t2, x2(:,4),'LineWidth',2)
ylabel('Pendulum Angular Velocity(rad/s)')
xlabel('Time(s)');
legend('Trial 1 - x2','Trial 2 - x2','Trial 1 - x4','Trial 2 - x4')
grid on
title('Cart Velocity and Pendulum Angular Velocity in 2 Cases');


%% initial conditions
ts = 0.01;
t = 0:ts:5;
%Case1
x01 = [0; 0.1; 0.1; 0]
C1 = [1 0 0 0; 0 0 1 0; K1];
G3 = ss(A1, B, C1, 0);
[yc1, t, xc1] = initial(G3,x01,t);
u = (K2 * xc1')';
figure(5)
subplot(4,1,1)
plot(t, xc1(:,1))
grid on
title('s vs t')
xlabel('t [s]')
ylabel('s [m]')

subplot(4,1,2)
plot(t, xc1(:,2))
grid on
title('v vs t')
xlabel('t [s]')
ylabel('v [m/s]')

subplot(4,1,3)
plot(t, xc1(:,3))
grid on
title('\phi vs t')
xlabel('t [s]')
ylabel('\phi [rad]')

subplot(4,1,4)
plot(t, xc1(:,4))
grid on
title('\phi vs t')
xlabel('t [s]')
ylabel('\phi [rad]')


% subplot(3,1,3)
% plot(t, u)
% grid on
% title('u vs t')
% xlabel('t [s]')
% ylabel('u [N]')



%Case2
x02 = [1; 0.1; 0.1; 0]
[yc2, t, xc2] = initial(G3,x02,t);
x03 = [0; 0.6; 0.1; 0]
[yc3, t, xc3] = initial(G3,x03,t);
x04 = [0; 0.1; 0.6; 1]
[yc4, t, xc4] = initial(G3,x04,t);


fig6 = figure('Renderer', 'painters', 'Position', [10 10 1200 500]);
subplot(1,2,1)
plot(t, xc1(:,2),t, xc2(:,2),t, xc3(:,2),t, xc4(:,2),'LineWidth',2)
hold on
yline(0,'-.b','Equilibrium v(t)');
ylabel('Cart Velocity (m/s)')
xlabel('Time(s)');
legend('Case 1','Case 2','Case 3','Case 4')
grid on
title('Cart Velocity in 4 Cases');
subplot(1,2,2)
plot(t, xc1(:,4),t, xc2(:,4),t, xc3(:,4),t, xc4(:,4),'LineWidth',2)
hold on
yline(0,'-.b','Equilibrium  \partial{\phi(t)}');
ylabel('Pendulum Angular Velocity (rad/s)')
xlabel('Time(s)');
legend('Case 1','Case 2','Case 3','Case 4')
grid on
title('Pendulum Angular Velocity in 4 Cases');
% 
% 
% figure(8)
% pzmap(G1)

% u=ones(size(t));
% y=lsim(G1,u,t);
% % y1=lsim(G,u,t);
% plot(t,y,'-')
% ylim([-0.0002 0.0002]);

% [y,t,x]=initial(G1,x0,t);
% plot(t,y)
% K = [k1 k2 k3 k4]







% X = s*eye(4)-A-B*K
% M = det(X)
% collect(M)
% 
% syms EPS
% ra=routh([1 (297*k4)/250 - k2 + 1 (297*k3)/250 - k1 + (547*k4)/250 - 233/20 (233*k2)/20 + (547*k3)/250 - 233/20 (233*k1)/20],EPS)
% 
% cond1 = (297*k4)/250 - k2 + 1 > 0;
% cond2 = -(k1 + k3 + (58261*k4)/5000 - k1*k2 + (297*k1*k4)/250 + (297*k2*k3)/250 + (547*k2*k4)/250 - (88209*k3*k4)/62500 - (162459*k4^2)/62500)/((297*k4)/250 - k2 + 1) > 0;
% cond3 = (((297*k4)/250 - k2 + 1)*(136750000*k1*k3 - 8484258125*k4 - 728125000*k3 + 865012500*k1*k4 - 136887500*k2*k3 + 6891120625*k2*k4 + 2621073200*k3*k4 + 1027634850*k1*k4^2 + 162459000*k2*k3^2 + 865012500*k2^2*k3 - 1892647350*k2*k4^2 + 1593137500*k2^2*k4 - 355460292*k3*k4^2 - 193001292*k3^2*k4 + 136750000*k3^2 + 1892647350*k4^2 - 136750000*k1*k2*k3 - 865012500*k1*k2*k4 + 162459000*k1*k3*k4 - 728425850*k2*k3*k4))/(250000*(297*k4 - 250*k2 + 250)*(k1 + k3 + (58261*k4)/5000 - k1*k2 + (297*k1*k4)/250 + (297*k2*k3)/250 + (547*k2*k4)/250 - (88209*k3*k4)/62500 - (162459*k4^2)/62500)) > 0;
% cond4 = (233*k1)/20 > 0;
% 
% conds = [cond1 cond2 cond3 cond4];
% 
% sol = solve(conds, [k1 k2 k3 k4], 'ReturnConditions', true);

% sol.conditions

function RA=routh(poli,epsilon);
%ROUTH   Routh array.
%   RA=ROUTH(R,EPSILON) returns the symbolic Routh array RA for 
%   polynomial R(s). The following special cases are considered:
%   1) zero first elements and 2) rows of zeros. All zero first 
%   elements are replaced with the symbolic variable EPSILON
%   which can be later substituted with positive and negative 
%   small numbers using SUBS(RA,EPSILON,...). When a row of 
%   zeros is found, the auxiliary polynomial is used.
%
%	Examples:
%
%	1) Routh array for s^3+2*s^2+3*s+1
%
%		>>syms EPS
%		>>ra=routh([1 2 3 1],EPS)
%		ra =
%
% 		   1.0000    3.0000
% 		   2.0000    1.0000
% 		   2.5000         0
% 		   1.0000         0
%
%	2) Routh array for s^3+a*s^2+b*s+c
%	
%		>>syms a b c EPS;
%		>>ra=routh([1 a b c],EPS);
%		ra =
% 
%		[          1,          b]
%		[          a,          c]
%		[ (-c+b*a)/a,          0]
%		[          c,          0]
%
% 
%   Author:Rivera-Santos, Edmundo J.
%   E-mail:edmundo@alum.mit.edu
%
if(nargin<2),
	fprintf('\nError: Not enough input arguments given.');
	return
end
dim=size(poli);		%get size of poli		
coeff=dim(2);				%get number of coefficients
RA=sym(zeros(coeff,ceil(coeff/2)));	%initialize symbolic Routh array 
for i=1:coeff,
	RA(2-rem(i,2),ceil(i/2))=poli(i); %assemble 1st and 2nd rows
end
rows=coeff-2;		%number of rows that need determinants
index=zeros(rows,1);	%initialize columns-per-row index vector
for i=1:rows,
	index(rows-i+1)=ceil(i/2); %form index vector from bottom to top
end
for i=3:coeff,				%go from 3rd row to last
	if(all(RA(i-1,:)==0)),		%row of zeros
			fprintf('\nSpecial Case: Row of zeros detected.');
			a=coeff-i+2;		%order of auxiliary equation
			b=ceil(a/2)-rem(a,2)+1; %number of auxiliary coefficients
			temp1=RA(i-2,1:b);	%get auxiliary polynomial
			temp2=a:-2:0;		%auxiliry polynomial powers
			RA(i-1,1:b)=temp1.*temp2;	%derivative of auxiliary
	elseif(RA(i-1,1)==0),		%first element in row is zero
			fprintf('\nSpecial Case: First element is zero.');
			RA(i-1,1)=epsilon;	%replace by epsilon
	end
				%compute the Routh array elements
	for j=1:index(i-2),	
		RA(i,j)=-det([RA(i-2,1) RA(i-2,j+1);RA(i-1,1) RA(i-1,j+1)])/RA(i-1,1);
	end
end
end