% Data Entry is a code that solves the displacement of regular and irregular
% structures taking into account the effect of live load  as a percentaje
% of the mass of the structure.

%% Clearing Command Window, Workspace and Plots
clc
clear
close all

%% Building Properties
% General Properties
nFloors = 10;                % Number of Floors
masspf = 0.3;


% Structure Geometry parameters
Lc  = 0.025;                 % [m] Length of the column in X direction
Bc  = 0.002;                 % [m] Breadth of the column in Y direction
hFloor = 0.165;                 % [m] Height of the floors/Columns

% Material Properties
E = 65;                     % GPa

Xi = 0.05;                 % [%] Damping
tri_lin = [1,1];        % tri-linear state of the material. first yield/elastic modulus, second yield/elastic modulus
SY1 = 240;                   % [MPa] First Yiled point
SY2 = 260;                   % [MPa] Second Yiled point

System.masspf= masspf;
System.n     = nFloors;
System.Lc    = Lc;
System.Bc    = Bc;
System.H     = hFloor;
System.E     = E;
System.TL    = tri_lin;
System.Xi    = Xi;
System.SY1   = SY1;
System.SY2   = SY2;

%Typical Square shaped building plan with side 3m
x_col = [0 0.45 0.45 0];          % x coordinate location
y_col = [0 0 0.315 0.315];          % y coordinate location
Geometry.x = x_col;
Geometry.y = y_col;

% Loading Ground Acceleration from the file
% R = 0.5                   % Energy dissipation coefficient NSR-10
% ga = load('GAcc.mat','ag');
% ga = struct2array(ga);
ga = xlsread('Ansys.xlsx');
GAx = ga;
GAy = ga;
Loading.signalX = GAx;
Loading.signalY = GAy;

% B_THARK Calculates dispalcement doing Time History Analisys 
[ux,uy,ur,t,Track]=B_THARK(System,Geometry,Loading);
% 

ux = ux';
uy = uy';

%% Export results to compare
% U = [t', ux'*100, uy'*100, ur'*100];
% writematrix(U,'Displacement.xlsm');
% xlswrite('Displacement.xlsm',U);

% %% PLOTS
% % 
track_b = Track.beta;
bC1L1(1,:) = track_b(1,1,:);
figure();
plot(bC1L1(1,:));
% 
Stress = Track.Stress;
sc1l1 = Stress(1,1,:);
for i = 1:length(sc1l1)
    stressplot(i,1) = Stress(1,1,i);
end

sc1l1 = Stress(1,1,:);
figure();
plot(sc1l1(1,:)/10^6);
xlabel('Time');ylabel('Stress in MPa');title('Stress');


%% Earthquake
% Imput Signal in X
% figure();
% subplot(1,2,1)
% plot(t,GAx);
% title('Acceleration Records');
% xlabel('Time [s]')
% ylabel('Acceleration in "g" [m/s^2]')
% legend('Ground Accl in "x" [m/s^2]','Ground Accl in "y" [m/s^2]');
% % 
% % Imput Signal in Y
% subplot(1,2,2)
% plot(t,GAy);
% title('Acceleration Records');
% xlabel('Time [s]')
% ylabel('Acceleration in "g" [m/s^2]')
% legend('Ground Accl in "y" [m/s^2]','Ground Accl in "y" [m/s^2]');
% % 
%% Displacements of the slab
% X Displacement
figure();
subplot(3,1,1)
plot(t,ux*100);
title('Slab displacement in "x" v. Time');
xlabel('Time [s]')
ylabel('Displacement in "x" [cm]')
legend; %('1','2','3','4','5','6');

% Y Displacement
subplot(3,1,2)
plot(t,uy*100);
title('Slab displacement in "y" v. Time');
xlabel('Time [s]')
ylabel('Displacement in "y" [cm]')
legend; %('1','2','3','4','5','6');

% Rotational Displacement
subplot(3,1,3)
plot(t,ur);
title('Slab rotation "{\theta}" v. Time');
xlabel('Time [s]')
ylabel('Rotation "{\theta}" [rad]');
legend; %('1','2','3','4','5','6');
% 

