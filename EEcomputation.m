function EEvalue = EEcomputation(SNR,lambda,m,k,gammaval,alpha,omegaSigma2,eta,epsilon,tau,A,C0,C1,D0,D1)
%Calculates the energy efficiency (EE) according to the theoretical
%formulas provided in Section IV of the article:
%
%Emil Bjornson, Luca Sanguinetti, Marios Kountouris, "Deploying Dense
%Networks for Maximal Energy Efficiency: Small Cells Meet Massive MIMO,"
%IEEE Journal on Selected Areas in Communications, to appear.
%
%Download article: http://arxiv.org/pdf/1505.01181.pdf
%
%This is version 1.0 (Last edited: 2016-01-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%This is version 1.0. 
% 
%INPUT: 
%SNR         = SNR value (rho/sigma^2)
%lambda      = BS density (BS/m)
%m           = Number of BS antennas
%k           = Number of users
%gammaval    = SINR value
%alpha       = Pathloss exponent
%omegaSigma2 = Propagation loss multiplied with noise variance
%eta         = Power amplifier efficiency
%epsilon     = Level of hardware impairments
%tau         = Length of coherence block (in symbols)
%A           = Power consumed by coding, decoding, and backhaul (in J/bit)
%C0          = Static energy consumption
%C1          = Circuit energy per active UE
%D0          = Circuit energy per active BS antenna
%D1          = Signal processing coefficient
% 
%OUTPUT: 
%EEvalue     = Energy Efficiency

%Compute the B1 and B2 from Eq. (18) and Eq. (19)
B1 = (4*k/(alpha-2)^2 + (k+m*(1-epsilon^2))/(alpha-1) + 2*(k+1/SNR)/(alpha-2));
B2 = (k+1/SNR + 2*k/(alpha-2))*(1+1/SNR) + (1-epsilon^2)*epsilon^2*m;

%Compute beta using Eq. (17)
beta = B1*gammaval / (m*(1-epsilon^2)^2-B2*gammaval);

%Check if the problem is feasible (satisfied the two constraints)
if beta>=1 && k*beta<=tau
    
    %Compute the ASE and AEC (both normalized by lambda) based on Eq. (10)
    %and Eq. (12).
    ASE = k *(1-beta*k/tau) * log2(1+ (1-epsilon^2)^2*m/ ( (k+1/SNR)*(1+2/(beta*(alpha-2))+1/SNR)+2*k/(alpha-2)*(1+1/SNR)+k/beta*(4/(alpha-2)^2 + 1/(alpha-1)) + (1-epsilon^2)*m/(alpha-1)/beta + (1-epsilon^2)*epsilon^2*m ));
    AEC = ( (1-(beta*k-1)/tau)*(SNR*omegaSigma2/eta)*gamma(alpha/2+1)/(pi*lambda)^(alpha/2) * k + C0 + C1*k + D0*m + D1*m*k) + A*ASE;
    
    EEvalue = (ASE/AEC)/1e6; %Compute EE normalized to 1e6 to get a better scale for optimization
    
else
    
    EEvalue = 0; %If the current problem was not feasible
    
end
