function [A,COMass,J] = f_polygeom(x_col,y_col) 
% f_polygeom calculates Area, Center of Mass, and the Torcional Constant
%  and then, plots the system's plan view with its Columns, Center of mass and
%  Center of Stiffness
tic
x_plot = x_col;
y_plot = y_col;
 
% check if inputs are same size
if ~isequal( size(x_col), size(y_col) )
  error( 'X and Y must be the same size');
end
%% 
% temporarily shift data to mean of vertices for improved accuracy
xm = mean(x_col);
ym = mean(y_col);

x_col = x_col - xm;
y_col = y_col - ym;
%%  
% summations for Counter Clock Wise boundary
x_shift = x_col([2:end 1]);
y_shift = y_col([2:end 1]);
%%
a = x_col.*y_shift - x_shift.*y_col;
%% 
A   = sum(a)/2;

xc  = sum((x_col+x_shift).*a)/6/A;
yc  = sum((y_col+y_shift).*a) /6/A;
Ixx = sum((y_col.*y_col +y_col.*y_shift + y_shift.*y_shift).*a)/12;
Iyy = sum((x_col.*x_col +x_col.*x_shift + x_shift.*x_shift).*a)/12;
Ixy = sum((x_col.*y_shift +2*x_col.*y_col +2*x_shift.*y_shift + x_shift.*y_col).*a )/24;
%%
dx = x_shift - x_col;
dy = y_shift - y_col;
P = sum(sqrt(dx.*dx+dy.*dy));
%%
% check for CCW versus CW boundary
if A < 0
   A = -A;
   Ixx = -Ixx;
   Iyy = -Iyy;
   Ixy = -Ixy;
end
%% 
% centroidal moments
Iuu = Ixx - A*yc*yc;
Ivv = Iyy - A*xc*xc;
Iuv = Ixy - A*xc*yc;
J = Iuu + Ivv;
%%
% replace mean of vertices
x_com = xc + xm;
y_com = yc + ym;
% Ixx = Iuu + A*y_com*y_com;
% Iyy = Ivv + A*x_com*x_com;
% Ixy = Iuv + A*x_com*y_com;
%%
% principal moments and orientation
I = [ Iuu  -Iuv ;
     -Iuv   Ivv ];
% [ eig_vec, eig_val ] = eig(I);
% I1 = eig_val(1,1);
% I2 = eig_val(2,2);
[eig_vec, ~] = eig(I);
ang1 = atan2( eig_vec(2,1), eig_vec(1,1) );
ang2 = atan2( eig_vec(2,2), eig_vec(1,2) );

%% Initial Center of Stiffness

n   = length(x_col);
csx = zeros(1,n);
csy = zeros(1,n);
% Assume all stiffnesses to be unity 
for col = 1:n
    csx(col) = (x_col(col) - x_col(1))*1;
    csy(col) = (y_col(col) - y_col(1))*1;
end

CSx = sum(csx)/length(x_col);
CSy = sum(csy)/length(y_col);

%%
% return values
COMass =  [A  x_com  y_com P];
% iner = [ Ixx  Iyy  Ixy  Iuu  Ivv  Iuv ];
% P_Inertia  =[ I1  ang1  I2  ang2  J ];

% constants
% d2r = pi / 180;
 
% % show results
% area = COMass(1);
% x_com = COMass(2);
% y_com = COMass(3);
% perimeter = COMass(4);
% disp( [ ' ' ] );
% disp( [ ' ' ] );
% disp( [ '      area     x_cen     y_cen     perim' ] );
% disp( [ area  x_com  y_com  perimeter ] );
%  
% I1 = P_Inertia(1);
% angle1 = P_Inertia(2);
% I2 = P_Inertia(3);
% angle2 = P_Inertia(4);
% disp( [ ' ' ] )
% disp( [ '        I1        I2' ] );
% disp( [ I1 I2 ] );
% disp( [ '    angle1    angle2' ] );
% disp( [ angle1/d2r angle2/d2r ] );
 
%% plot outline plan
x_col = x_plot;
y_col = y_plot;
xplot = x_col([1:end 1]);
yplot = y_col([1:end 1]);
rad = 15;         % Limits of the function

x1 = [x_com-rad*cos(ang1)  x_com+rad*cos(ang1)];
y1 = [y_com-rad*sin(ang1)  y_com+rad*sin(ang1)];
x2 = [x_com-rad*cos(ang2)  x_com+rad*cos(ang2)];
y2 = [y_com-rad*sin(ang2)  y_com+rad*sin(ang2)];

plot(xplot,yplot,'b',x1,y1,'g:', x2,y2,'g:');
hold on
p1=plot(x_plot,y_plot,'s','MarkerSize',10,'MarkerFaceColor',[1 1 0]);
p2=plot(x_com, y_com,'o');
p3=plot(CSx,CSy,'*','MarkerSize',10);  
hold off
axis( [ -rad/2  rad/2  -rad/2  rad/2 ] )
axis square
title('Building Plan');
legend([p1 p2 p3],{'Columns','Center of Mass','Center of Stiffness'})

end