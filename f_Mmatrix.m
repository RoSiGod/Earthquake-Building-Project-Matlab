function[M] = f_Mmatrix(nFloors,Mass)
% f_Mmatrix calculates the Mass matrix based on the number of floors
% Exploiting the band matrix properties of the Mass, the code
% attempts to find indices of a matrix which has non zero entries.
% each floor has 3 DOFs and hence, n stories will have 3n DOFs.

% Each floor has a Mass contrubution in X, Y and a rotational Inertia

% let M be [M0; M0; I0, M1; M1; I1,... Mn; Mn; In]
% where n is the number of floors


for j = 1:nFloors
    Me(:,:,j) = [Mass(1,j),  0,     0;
                   0,   Mass(2,j),  0;
                   0,     0,   Mass(3,j)];
end


for j = 1:nFloors
    M(1+3*(j-1):3+3*(j-1), 1+3*(j-1):3+3*(j-1)) = Me(:,:,j);
end

end
