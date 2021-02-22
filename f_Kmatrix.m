function[K] = f_Kmatrix(nFloors,Stiffx,Stiffy,Stiffr,ex,ey)
% f_Kmatrix calculates the Stiffness matrix based on number of floors.
% Exploiting the band matrix properties of the stiffness the code
% attempts to find indices of a matrix which has non zero entries. 
% Each floor has 3 DOFs and hence, n stories will have 3n DOFs. 

% Each floor has X and Y lateral sitffness and a rotational stiffness about Z.
% Rotational stiffess has the term 'e' eccentricity embedded into it.

% Kx, Ky, Kr are lateral X, Y and rotation stiffess
% ex, ey, e are eccentricities in x, y, and net exentricity
% n is the number of floors
% let Kx = [Kx0, Kx1, Kx2,... Kxn];
Kx = Stiffx;
Ky = Stiffy;
Kr = Stiffr;
 
for j = 1:nFloors-1
    
    Ke(:,:,j) = [  Kx(j)+Kx(j+1),           0,                   -ey(j)*(Kx(j)+Kx(j+1));
                        0,             Ky(j)+Ky(j+1),             ex(j)*(Ky(j)+Ky(j+1));
           -ey(j)*(Kx(j)+Kx(j+1)), ex(j)*(Ky(j)+Ky(j+1)), (ex(j)^2*(Ky(j)+Ky(j+1))+ey(j)^2*(Kx(j)+Kx(j+1)))+Kr(j)+Kr(j+1)];
end
% Elemental Matrix
Ke(:,:,nFloors) = [Kx(nFloors),                   0,                   -ey(nFloors)*(Kx(nFloors));
                       0,                   Ky(nFloors),                ex(nFloors)*(Ky(nFloors));
           -ey(nFloors)*(Kx(nFloors)), ex(nFloors)*(Ky(nFloors)), (ex(nFloors)^2*Ky(nFloors)+ey(nFloors)^2*Kx(nFloors))+Kr(nFloors)];

% Initialize K as a zero matrix
K = zeros(3*nFloors,3*nFloors);

% Global Assembly
for j = 0:nFloors-1
    K(1+3*j:3+3*j, 1+3*j:3+3*j) = Ke(:,:,j+1);
end

% Connectivity between adjacent elements
for i=1:nFloors-1
    index_x = 3*(i-1);
    index_y = 3+3*(i-1);
    
    % Lateral to Lateral
    K(index_x+1, index_y+1) = -Kx(i+1);
    K(index_x+2, index_y+2) = -Ky(i+1);
    K(index_x+3, index_y+3) = -Kr(i+1);
    
    K(index_y+1, index_x+1) = -Kx(i+1);
    K(index_y+2, index_x+2) = -Ky(i+1);
    K(index_y+3, index_x+3) = -Kr(i+1);
    
    % Torsional to lateral
    K(index_y+1, index_y)   =  Kx(i+1)*(ey(i+1));
    K(index_y+2, index_y)   = -Ky(i+1)*(ex(i+1));
    
    K(index_y, index_y+1)   =  Kx(i+1)*(ey(i+1));
    K(index_y, index_y+2)   = -Ky(i+1)*(ex(i+1));
    
    % Lateral to Torsional
    K(index_x+1, index_y+3) =  Kx(i+1)*ey(i+1);
    K(index_x+2, index_y+3) = -Ky(i+1)*ex(i+1);
    
    K(index_y+3, index_x+1) =  Kx(i+1)*ey(i+1);
    K(index_y+3, index_x+2) = -Ky(i+1)*ex(i+1);
    
end

end
