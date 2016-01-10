%This Matlab script can be used to generate Figure 5 in the article:
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


%Initialization
close all;
clear all;


%%Simulation parameters


%Select the values of log2(1+gamma) that should be considered
rateValues = [1 2 3];

%Propagation parameters
alpha = 3.76; %Pathloss exponent
tau = 400; %Length of coherence block (in symbols)

%Hardware characterization
eta = 0.39; %Power amplifier efficiency
epsilonValues = 0:0.01:0.2; %Range of levels of hardware impairments

%Spectral resources
T = 1/(2e7); %Symbol time (based on 20 MHz)

%Energy parameters
A = 1.15e-9; %Power consumed by coding, decoding, and backhaul (in J/bit)
C0 = 10 * T; %Static energy consumption (10 W divided over the symbols)
C1 = 0.1 * T; %Circuit energy per active UE
D0 = 0.2 * T; %Circuit energy per active BS antenna
D1 = 1.56e-10; %Signal processing coefficient


%Maximal number of antennas and users considered in the simulation. These
%numbers need to selected carefully so that the maximal value is not at the
%edge of the considered region
Mmax = 200;
Kmax = 30;

%Placeholders for storing of simulation results
EE_theory = zeros(length(epsilonValues),length(rateValues)); %EE for different epsilon and gamma values (using theoretical formulas)


%Go through all gamma values
for j = 1:length(rateValues)
    
    %Extract the current gamma value
    gammaval = 2^rateValues(j) - 1;
    
    %Go through all epsilon values
    for n = 1:length(epsilonValues)
        
        %Extract the current epsilon value
        epsilon = epsilonValues(n);
        
        %Prepare to store the best parameters for each (M,K)-value
        EEtmp = zeros(Mmax,Kmax);
        
        %Go through range of K values
        for k = 1:Kmax
            
            %Go through range of M values
            for m = 1:Mmax
                
                %Compute B1bar and B2bar from Eq. (18) and Eq. (19), respectively.
                B1bar = k*(4/(alpha-2)^2 + 1/(alpha-1) + 2/(alpha-2)) + m*(1-epsilon^2)/(alpha-1);
                B2bar = k*(1 + 2/(alpha-2)) + m*(1-epsilon^2)*epsilon^2;
                
                
                %Check if the two constraints in Eq. (23) are satisfied
                if B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval)>=1 && (k*B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval))<=tau %gammaval>=m*(1-epsilon^2)^2/(B1+B2) &&
                    
                    %Compute the objective function in Eq. (23) and store the EE
                    ASE = k *(1-(B1bar*gammaval / (m*(1-epsilon^2)^2-B2bar*gammaval))*k/tau) * log2(1+ gammaval);
                    AEC = C0 + C1*k + D0*m + D1*m*k + A*ASE;
                    EEtmp(m,k) = ASE/AEC;
                    
                end
                
                
            end
            
        end
        
        %Find the M and K values that maximize the EE
        [EEmaxM,optM] = max(EEtmp,[],1);
        [EEmax,optK] = max(EEmaxM);
        
        %Store the maximal EE
        EE_theory(n,j) = EEmax;
        
    end
    
end


%Plot Figure 5 from the paper
figure; hold on; box on;

plot(epsilonValues,EE_theory(:,1)/1e6,'r--','LineWidth',1);
plot(epsilonValues,EE_theory(:,2)/1e6,'b-','LineWidth',1);
plot(epsilonValues,EE_theory(:,3)/1e6,'k-.','LineWidth',1);

xlabel('Level of hardware impairments (\epsilon)');
ylabel('Energy efficiency [Mbit/Joule]');
ylim([0 12]);

legend('\gamma = 1','\gamma = 3','\gamma = 7','Location','SouthWest');
