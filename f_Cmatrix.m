function[C] = f_Cmatrix(M,K,Xi,nModes)
% f_Cmatrix calculates Damping matrix according to Rayleigh damping criteria


if nModes>length(M)
    error('Error = nModes MUST be smaller or equal to Maximum DoF');
end

[phi,eval] = eig(K,M);

Iphi  = inv(phi);
Tphi  = transpose(phi);
ITphi = inv(Tphi);

AngFreq  = sqrt(diag(eval));
M = Tphi*M*phi;
K = Tphi*K*phi;

for i=1:length(M)
    
    if i<= nModes
        Epsilon(i) = (0.2 - Xi)/(AngFreq(nModes) - AngFreq(1))*(AngFreq(i) - AngFreq(1)) + Xi;
    else
        Epsilon(i) = (0.2 - Xi)/(AngFreq(nModes) - AngFreq(1))*(AngFreq(i) - AngFreq(nModes)) + Xi;
    end
    
end

B1 = (2*Epsilon(1)*AngFreq(1) - 2*Epsilon(nModes)*AngFreq(nModes))/(AngFreq(1)^2 - AngFreq(nModes)^2);
A1 = 2*Epsilon(1)*AngFreq(1) - B1*AngFreq(1)^2;

C = A1*M + B1*K;

C = C*Iphi;
C = ITphi*C;

end
