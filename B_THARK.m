function [ux,uy,ur,t,Track] = B_THARK(System,Geometry,Loading)
% Time History Análisis using 4th order Runge Kutta Method to solve
% [mDL+mLL]*[u'']+[C]*[u']+[K]*[u]=[mDL+mLL]*[ga]

%% Clearing Command Window, Workspace and Plots
clc
clear
close all

%% Building Properties
% General Properties
nFloors = 10;                % Number of Floors
masspf = 0.300000000;


% Structure Geometry parameters
Lc  = 0.025;                 % [m] Length of the column in X direction
Bc  = 0.003;                 % [m] Breadth of the column in Y direction
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


%%
tic
% Retrieve Building Properties
% General Properties
xcol    = Geometry.x;             % Column coordinates in a Counter-ClockWise sense
ycol    = Geometry.y;
nFloors = System.n;               % Number of Floors
masspf  = System.masspf;          % [kg] mass per floor

% Geometry of the elements that makes the structure
Lc     = System.Lc;               % [m] Length of the column in X direction
Bc     = System.Bc;               % [m] Breadth of the column in Y direction
H      = System.H;                % [m] Height of floors/Columns

% Material Properties
tri_lin = System.TL;              % 1 X 2 array of tri-linear state of the material. first yield/elastic modulus, second yield/elastic modulus
Xi      = System.Xi;              % [%] Damping
SY1     = System.SY1;             % [MPa] First Yiled point
SY2     = System.SY2;             % [MPa] Second Yiled point
E       = System.E;               % [GPa] Elastic Modulus

% Signal Properties
GAx  = Loading.signalX;
GAy  = Loading.signalY;
Nx   = length(GAx);
Ny   = length(GAy);
L_ga = max(Nx,Ny);

% Define Time array for analysis
Ti = 0;
Dt = 0.0001;
% Dt = 0.01;
Tf = 0.1;
t  = Ti:Dt:Tf;

%% Structural Data
% ---------------------
% Geometry Properties
% ---------------------

E   = E*10^9;                    % [N/m^2] Elastic modulus of concrete
G   =  E/2;                      % [N/m^2] Shear modulus of concrete
% CA  = Lc*Bc;                   % [m^2] Cross sectional area of columns
J   = 0.18*Bc^3*Lc;              % [m^4] Torsional Constant
Ixx = Bc*Lc^3/12;                % [m^4] Inertia along X axis
Iyy = Lc*Bc^3/12;                % [m^4] Inertia along Y axis

%% System Data
% ---------------------
%         MASS
% ---------------------

% Polygeom function calculates Area, Centroid, Geometry and 
% Inertia of the polygon formed by index vectors xcol and ycol
[Area,COMass,P_inertia] = f_polygeom(xcol,ycol);

Mslab = masspf;         % [kg] Mass per m^2 on each floor

COM0(1,1) = COMass(2);            % X-COM of the building plan
COM0(2,1) = COMass(3);            % Y-COM of the building plan

mI = 0.0011;       % [kg/m^2] Mass Inertia of the selected geometry MUST BE 1 X nFloors dimenssional array

% ---------------------
%       STIFFNESS
% ---------------------
nCol = length(xcol);              % Number of columns per floor
beta = zeros(nCol,nFloors) + 1;   % Defines the coefficient of shear modulus in each computation default set to 1.
beta_y1 = zeros(nCol,nFloors) + 1;% Initiates the column stiffness multiplier to 1 (to be changed later according to 1st yielding)
beta_y2 = zeros(nCol,nFloors) + 1;% Initiates the column stiffness multiplier to 1 (to be changed later according to 2nd yielding)
                                  % If a particular columns is found to yield, the coefficient is reduced to account for the reduction of stiffness
Kx = beta*12*E*Iyy/H^3;           % X directional Stiffness
Ky = beta*12*E*Ixx/H^3;           % Y directional Stiffness
Kr = beta.*(Kx.*abs((ycol'-COM0(2,1)))+Ky...
    .*abs((xcol'-COM0(1,1))));    % Torsional Stiffness

% ---------------------
%   FORCING FUNCTION
% ---------------------
F = zeros(3*nFloors,length(t)); 
 
% Loop to set forcing function in x and y directions
for fm = 1:nFloors
%     F(1 + 3*(fm-1),:) = 175000*(1 - exp(-10.*t));
    F(1 + 3*(fm-1),:) = 0;
    F(2 + 3*(fm-1),:) = 0;
end

%% System Variables

% Initializing Vectors

Stiffx = zeros(1,nFloors);        % Stiffness in x
Stiffy = zeros(1,nFloors);        % Stiffness in y
Stiffr = zeros(1,nFloors);        % Rotational stiffness

Track_stress = zeros (nCol,nFloors,length(t)-1);
Track_beta   = zeros (nCol,nFloors,length(t)-1);

CSx = zeros(nFloors,nCol);
CSy = zeros(nFloors,nCol);
CSX = zeros(1,nFloors);
CSY = zeros(1,nFloors);

Mark_i=zeros(nCol,nFloors);

% Location of Center Of Mass (COM) at the start of the calculation
COM(1,:) = (zeros(1,nFloors) + 1).*COM0(1,1);
COM(2,:) = (zeros(1,nFloors) + 1).*COM0(2,1);

X     = zeros(3*nFloors,3);       % [[Displacemets],[Velocities],[Accelerations]] ie: [[Ux;Uy;Ur],[Vx;Vy;Vr],[Ax;Ay;Ar]]
U     = zeros(3*nFloors,1);       % Vector of displacements
Ux    = zeros(nCol,nFloors);      % Vector of displacements in x
Uy    = zeros(nCol,nFloors);      % Vector of displacements in y
Uy1   = zeros(3*nFloors,1);       % First Yield displacements
Uy2   = zeros(3*nFloors,1);       % 2nd Yield displacements

for i = 1:nFloors
    U(1+3*(i-1),1) = 0;
    U(2+3*(i-1),1) = 0.001;
    U(3+3*(i-1),1) = 0;
end
v   = X(:,2);                     % Vector of velocities
a   = X(:,3);                     % Vector of accelerations
Ax0 = zeros(1,length(t));         % Vector of accelerations in x
Ay0 = zeros(1,length(t));         % Vector of accelerations in y



%% TIME HISTORY
for i=1:(length(t)-1)
    % See that the block is confined withing the boundary
    % (Taken from: www.mathworks.com/help/matlab/ref/inpolygon.html)
    % Before that, Building Plan vertices and Blocks coordinates at any
    % given time 't' need to be determined with respect to moving COM
    % (since the rotation is defined about the COM)
    
    if i~=1
        for j=1:nFloors
            
            % co-ordinates of polygon vertices with respect to moving COM
            xcol_modified = (xcol + U(1 + 3*(j-1)) - COM(1,j));
            ycol_modified = (ycol + U(2 + 3*(j-1)) - COM(2,j));
            % taking rotation of the plan into account
            xcol_modified = xcol_modified*cos(U(3 + 3*(j-1)));
            ycol_modified = ycol_modified + xcol_modified*tan(U(3 + 3*(j-1)));
                        
                      
            % inout==0 indicates that the block has crossed the boundary, when the
            % previous displacement ((i-1)th iteration) coordinates are re-assigned to the current iteration "i'th" dispalcement value
            
        end
    end

    
%% Analyse the forces in each columns
% f_displacement calculates displacement of the columns in deformed configuration

    for j = 1:nFloors
        [Ux(:,j),Uy(:,j)] = f_displacement(U(:,i),[xcol' ycol'],COM,j-1);
    end
    
% f_stress calculates column Stresses and Strains
% Stress array has number of columns as rows and number of floors as columns

VStress = zeros(nCol,nFloors);    % Von-Mises Stress
    
    r = zeros(nCol,1);
    for j = 1:nFloors
        
        % r is rotation
        r(1:nCol,j) = U(3*(j),i).*(zeros(nCol,1) + 1);
        if j==1
            VStress(:,j) = f_stress(Ux(:,j),Uy(:,j),r(:,j),0,0,0,H,E,G,J,Ixx,Iyy,Lc,Bc);
            
        else
            VStress(:,j) = f_stress(Ux(:,j),Uy(:,j),r(:,j),Ux(:,j-1),Uy(:,j-1),r(:,j-1),H,E,G,J,Ixx,Iyy,Lc,Bc);
           
        end
    end
    Track_stress(:,:,i) = VStress(:,:);     % just for ploting purpuses

% Check if any column yields according to Von-Mises Stress criteria
    for j = 1:nFloors
        for k = 1:nCol
            if Mark_i(k,j)==0
                if SY2*10^6 >= VStress(k,j) && VStress(k,j)> SY1*10^6
                    Mark_i(k,j) = 100;                                                                           
                    beta_y1(k,j) = tri_lin(1); 
                    % 1st Yield displacements (storey wise)
                    Uy1(1+3*(j-1)) = U(1+3*(j-1),i);
                    Uy1(2+3*(j-1)) = U(2+3*(j-1),i);
                    Uy1(3+3*(j-1)) = U(3+3*(j-1),i);                 
                    
                elseif VStress(k,j)>SY2*10^6
                    Mark_i(k,j) = 10000;                    
                    beta_y2(k,j) = tri_lin(2);
                    % 2nd Yield displacements (storey wise)
                    Uy2(1+3*(j-1)) = U(1+3*(j-1),i);
                    Uy2(2+3*(j-1)) = U(2+3*(j-1),i);
                    Uy2(3+3*(j-1)) = U(3+3*(j-1),i);
                else                    
                    Mark_i(k,j) = 0;
                end
            elseif Mark_i(k,j)==100
                if VStress(k,j)>SY2*10^6                                                                            
                    beta_y2(k,j) = tri_lin(2);           
                    Mark_i(k,j) = 10000;
                    % 2nd Yield displacements (storey wise)
                    Uy2(1+3*(j-1)) = U(1+3*(j-1),i);
                    Uy2(2+3*(j-1)) = U(2+3*(j-1),i);
                    Uy2(3+3*(j-1)) = U(3+3*(j-1),i);
                    
                else
                    beta_y1(k,j) = tri_lin(1);             
                    Mark_i(k,j) = 100;
                end
            elseif Mark_i(k,j)==10000
                beta_y2(k,j) = tri_lin(2);                
                Mark_i(k,j) = 10000;
            end
        end
    end
    
    % Setup a counter to judge if the yielding has occured
    
    % 'Mark' sums up the counters initiated in the above if-for loops. if
    % the sum of all Marks =0, no yielding anywhere has taken place. If the
    % sum of all the Marks is anywhere between 1 to 9999, first yielding
    % has occured in one or more locations. If the sum of Marks crosses
    % 10000 to any value more than that, that would indicate that 2nd
    % yielding has occured at one or more places. 
    
    for j = 1:nFloors
    Mark_j(j) = sum(Mark_i(:,j));
    end
    Mark = sum(Mark_j(:));
    
    Track_beta(:,:,i) = beta_y1;

%% Track the Centre of Stiffness (CS) wrt to the coordinate system at the center 
% of the slab (in deformed configuration)
    % Centre of Stiffness for rectangular plan area

    for j = 1:nFloors
        for k = 1:nCol
            CSx(j,k) = (xcol(k) - xcol(1))*Ky(k,j);
            CSy(j,k) = (ycol(k) - ycol(1))*Kx(k,j);
        end
        CSX(j) = sum(CSx(j,:))/sum(Ky(:,j)) + xcol(1) - COM0(1);
        CSY(j) = sum(CSy(j,:))/sum(Kx(:,j)) + ycol(1) - COM0(2);
    end
    CS(1,:) = CSX ;
    CS(2,:) = CSY ;
    
%% Define e, e_x and e_y  eccentricites
% Distance of the COM of the system from the CS
    
    ex = zeros(1,nFloors);
    ey = zeros(1,nFloors);
    e = zeros(1,nFloors);
    
    % Track_ex(i,:) = ex;
    
%% Sytem Parameters
% ---------------------
%      MASS MATRIX
% ---------------------
    
    Mass = (zeros(3,nFloors)) + [Mslab*(zeros(1,nFloors)+1); Mslab*(zeros(1,nFloors)+1); mI*(zeros(1,nFloors)+1)];
    % f_Mmatrix calculates the Mass matrix of the structure 
    M = f_Mmatrix(nFloors,Mass);
    
    % Track_Mmat(i) = det(M);
    
% ---------------------
%    STIFFNESS MATRIX
% ---------------------   
% 1st - In-Elastic Stiffness Matrix components --------------------------
    Kx_ie1 = beta_y1*12*E*Iyy/H^3;       % X directional Stiffness
    Ky_ie1 = beta_y1*12*E*Ixx/H^3;       % Y directional Stiffness
    Kr_ie1 = beta_y1.*(Kx.*abs((ycol'-COM(2,:)))+Ky...
        .*abs((xcol' - COM(1,:))));  % Rotational Stiffness
  
    for j = 1:nFloors
        Stiffx_ie1(1,j) = sum(Kx_ie1(:,j));
        Stiffy_ie1(1,j) = sum(Ky_ie1(:,j));
        Stiffr_ie1(1,j) = sum(Kr_ie1(:,j));
    end    
% 1st In-Elastic Stiffness Matrix Assembly
   K_ie1 = f_Kmatrix(nFloors,Stiffx_ie1,Stiffy_ie1,Stiffr_ie1,ex,ey); 
   
% 2nd - In-Elastic Stiffness Matrix components --------------------------
    Kx_ie2 = beta_y2*12*E*Iyy/H^3;       % X directional Stiffness
    Ky_ie2 = beta_y2*12*E*Ixx/H^3;       % Y directional Stiffness
    Kr_ie2 = beta_y2.*(Kx.*abs((ycol'-COM(2,:)))+Ky...
        .*abs((xcol' - COM(1,:))));  % Rotational Stiffness
  
    for j = 1:nFloors
        Stiffx_ie2(1,j) = sum(Kx_ie2(:,j));
        Stiffy_ie2(1,j) = sum(Ky_ie2(:,j));
        Stiffr_ie2(1,j) = sum(Kr_ie2(:,j));
    end    
% 2nd In-Elastic Stiffness Matrix Assembly
   K_ie2 = f_Kmatrix(nFloors,Stiffx_ie2,Stiffy_ie2,Stiffr_ie2,ex,ey);
   
% Elastic Stiffness Matrix components -----------------------------------
      for j = 1:nFloors
        Stiffx(1,j) = sum(Kx(:,j));
        Stiffy(1,j) = sum(Ky(:,j));
        Stiffr(1,j) = sum(Kr(:,j));
      end
% Elastic Stiffness Matrix Assembly 
    K = f_Kmatrix(nFloors,Stiffx,Stiffy,Stiffr,ex,ey);

% ---------------------
%    DAMPING MATRIX
% ---------------------
    % f_Cmatrix calculates the damping matrix
    nModes = 3*nFloors;
    C = f_Cmatrix(M,K,Xi,nModes);

    
    
%% Runge Kutta 4th Order Method
     h = Dt;
    % Check if the structure is yielded or not
 
      if Mark ==0
        % Time Step Size

        uk1 = v(:,i);
        vl1 = a(:,i);
        uk2 = uk1 + h/2*vl1;
        vl2 = M\((F(:,i)+ F(:,i+1))/2 - K*(U(:,i) + uk1*h/2) - C*(v(:,i) + vl1*h/2));
        uk3 = uk2 + h/2*vl2;
        vl3 = M\((F(:,i)+F(:,i+1))/2   - K*(U(:,i) + uk2*h/2) - C*(v(:,i) + vl2*h/2));
        uk4 = uk3 + h*vl3;
        vl4 = M\((F(:,i) + F(:,i+1))/2  - K*(U(:,i) + uk3*h) - C*(v(:,i) + vl3*h));

        U(:,i+1) = U(:,i) + 1/6*h*(uk1 + 2*uk2 + 2*uk3 + uk4);
        v(:,i+1) = v(:,i) + 1/6*h*(vl1 + 2*vl2 + 2*vl3 + vl4);
        a(:,i+1) = M\(F(:,i+1) - K*U(:,i+1) - C*v(:,i+1));
        
      elseif Mark < 9999
          
        uk1 = v(:,i);
        vl1 = a(:,i);
        uk2 = uk1 + h/2*vl1;
        vl2 = M\((F(:,i)+ F(:,i+1))/2  - K_ie1*(U(:,i) + uk1*h/2 - Uy1) - K*(Uy1) - C*(v(:,i) + vl1*h/2));
        uk3 = uk2 + h/2*vl2;
        vl3 = M\((F(:,i)+F(:,i+1))/2  - K_ie1*(U(:,i) + uk2*h/2 - Uy1) - K*(Uy1) - C*(v(:,i) + vl2*h/2));
        uk4 = uk3 + h*vl3;
        vl4 = M\((F(:,i) + F(:,i+1))/2  - K*(U(:,i) + uk3*h - Uy1) - K*(Uy1) - C*(v(:,i) + vl3*h));

        U(:,i+1) = U(:,i) + 1/6*h*(uk1 + 2*uk2 + 2*uk3 + uk4);
        v(:,i+1) = v(:,i) + 1/6*h*(vl1 + 2*vl2 + 2*vl3 + vl4);
        a(:,i+1) = M\(F(:,i+1)  - K_ie1*(U(:,i+1) - Uy1) - K*(Uy1) - C*v(:,i+1));
        
      else
        uk1 = v(:,i);
        vl1 = a(:,i);
        uk2 = uk1 + h/2*vl1;
        vl2 = M\((F(:,i)+ F(:,i+1))/2  - K_ie2*(U(:,i) + uk1*h/2 - Uy1 - Uy2) - K*(Uy1) - K_ie1*(Uy2 - Uy1) - C*(v(:,i) + vl1*h/2));
        uk3 = uk2 + h/2*vl2;
        vl3 = M\((F(:,i)+F(:,i+1))/2   - K_ie2*(U(:,i) + uk2*h/2 - Uy1 - Uy2) - K*(Uy1) - K_ie1*(Uy2 - Uy1) - C*(v(:,i) + vl2*h/2));
        uk4 = uk3 + h*vl3;
        vl4 = M\((F(:,i) + F(:,i+1))/2  - K_ie2*(U(:,i) + uk3*h - Uy1 - Uy2) - K*(Uy1) - K_ie1*(Uy2 - Uy1) - C*(v(:,i) + vl3*h));

        U(:,i+1) = U(:,i) + 1/6*h*(uk1 + 2*uk2 + 2*uk3 + uk4);
        v(:,i+1) = v(:,i) + 1/6*h*(vl1 + 2*vl2 + 2*vl3 + vl4);
        a(:,i+1) = M\(F(:,i+1)  - K_ie2*(U(:,i+1) - Uy1 - Uy2) - K*(Uy1) - K_ie1*(Uy2 - Uy1) - C*v(:,i+1));
      end
          
   
    
%% 
end

%% Eigenvalues to Compare with Modal Analysis in SOFTWARE
[~,eigVs] = eig(K,M);
Frequencies = sqrt(diag(eigVs))/(2*pi);
Period = 1/Frequencies;
display(Period);

%% Slab displacement
ux  = zeros(nFloors,length(t));
uy  = zeros(nFloors,length(t));
ur  = zeros(nFloors,length(t));


for j = 1:nFloors
    ux(j,:)  = U(1+3*(j-1),:);
    uy(j,:)  = U(2+3*(j-1),:);
    ur(j,:)  = U(3+3*(j-1),:); 
   
end

%% Tracking Pointers
Track.beta = Track_beta;
Track.Stress = Track_stress;

toc

% end

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
end
