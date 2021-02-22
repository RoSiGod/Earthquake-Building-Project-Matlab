% f_displacement calculates the displacement of any point in the plane in
% a deformed configuration

function[Dx,Dy] = f_displacement(u,x,COM,nFloors)

Dx = u(1+3*nFloors) - u(3+3*nFloors).*(COM(2,nFloors+1) - x(:,2));
Dy = u(2+3*nFloors) + u(3+3*nFloors).*(COM(1,nFloors+1) - x(:,1));

end