
% Following functions calculate Displacement, Velocity and Accelerations of
% any point in the plane
% Inputs are V_any(X-Accl of COM, Y-Accl of COM, Angular Accl
% of COM, x-coordinate of the POI, y-coordinate of POI)
function[Ax,Ay] = f_acceleration(a,nFloors)

Ax = a(1+3*nFloors) ;
Ay = a(2+3*nFloors) ;

end
