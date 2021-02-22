% f_stress calculates Von-Mises stresses in each columns per floor
% Input as X displacement, Y displacement, rotation and Height of the Column

function[StressFloor] = f_stress(Dx, Dy, r0, Dxb, Dyb, Drb, H, E, G, J,Ixx,Iyy,Lc,Bc)


Dx = (Dx - Dxb);
Dy = (Dy - Dyb);
r0 = (r0 - Drb);

X_Stress_Bending = Dx.*(6*E/H^2)*Lc/2;
Y_Stress_Bending = Dy.*(6*E/H^2)*Bc/2;
Rot_Stress = r0/H*G.*J;
X_Stress_Shear = Dx.*(12*E*Iyy/H^3)/(Lc*Bc);
Y_Stress_Shear = Dy.*(12*E*Ixx/H^3)/(Lc*Bc);


Sigma_xx = 0;
Sigma_yy = 0;
Sigma_zz = X_Stress_Bending + Y_Stress_Bending;
Tau_zx = X_Stress_Shear;
Tau_zy = Y_Stress_Shear;
Tau_xy = Rot_Stress;

% Stress = [Sigma_x, Sigma_y, Tor_z];
StressFloor = sqrt(1/2*(Sigma_xx.^2 + Sigma_yy.^2 + Sigma_zz.^2) + 3*(Tau_zx.^2 + Tau_zy.^2 + Tau_xy.^2));

end