%This Matlab script can be used to generate Figure 3 and Figure 4 in the article:
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

%Select the value of log2(1+gamma) that should be considered
rateValues = 2;
gammaval = 2^rateValues - 1;

%Propagation parameters
alpha = 3.76; %Pathloss exponent
tau = 400; %Length of coherence block (in symbols)

%Hardware characterization
eta = 0.39; %Power amplifier efficiency
epsilon = 0.05; %Level of hardware impairments

%Spectral resources
T = 1/(2e7); %Symbol time (based on 20 MHz)

%Energy parameters
A = 1.15e-9; %Power consumed by coding, decoding, and backhaul (in J/bit)
C0 = 10 * T; %Static energy consumption (10 W divided over the symbols)
C1 = 0.1 * T; %Circuit energy per active UE
D0 = 0.2 * T; %Circuit energy per active BS antenna
D1 = 1.56e-10; %Signal processing coefficient


%Maximal number of antennas and users considered when plotting the EE.
Mmax = 100;
Kmax = 15;


%Placeholders for storing of simulation results
EEtheory = zeros(Mmax,Kmax); %EE for different lambda and gamma values (using theoretical formulas)
betaOptimal = zeros(Mmax,Kmax); %Store optimal beta for each point using the theoretical formulas


%Go through all number of users
for k = 1:Kmax
    
    %Go through all number of antennas
    for m = 1:Mmax
        
        %Compute B1bar and B2bar from Eq. (18) and Eq. (19), respectively.
        B1bar = k*(4/(alpha-2)^2 + 1/(alpha-1) + 2/(alpha-2)) + m*(1-epsilon^2)/(alpha-1);
        B2bar = k*(1 + 2/(alpha-2)) + m*(1-epsilon^2)*epsilon^2;
        
        %Check if the two constraints in Eq. (23) are satisfied
        if B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval)>=1 && (k*B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval))<=tau

            %Compute the objective function in Eq. (23) and store the EE
            ASE = k *(1-(B1bar*gammaval / (m*(1-epsilon^2)^2-B2bar*gammaval))*k/tau) * log2(1+ gammaval);
            AEC = C0 + C1*k + D0*m + D1*m*k + A*ASE;
            EEtheory(m,k) = ASE/AEC;
            
            %Compute and store the optimal beta according to Eq. (17) using
            %B1bar and B2bar since we consider the asymptotic regime.
            betaOptimal(m,k) = B1bar*gammaval / (m*(1-epsilon^2)^2-B2bar*gammaval);
            
        end
        
    end
    
end





%%Run the alternating optimization algorithm from Section IV.D

%Number of iterations
iterations = 20;

%Prepare to store the progress (notice that the first element in the
%starting point)
Kiterations = zeros(iterations+1,1);
Miterations = zeros(iterations+1,1);
EEiterations = zeros(iterations+1,1);

%Initiate the algorithm
Miterations(1) = 20;
Kiterations(1) = 1;
EEiterations(1) = EEtheory(Miterations(1),Kiterations(1));

%Go through all iterations
for iter = 2:iterations+1

    
    %Step 2: Update K using Theorem 3
    
    %Compute current cbar
    cbar = Miterations(iter-1)/Kiterations(iter-1);
    
    %Compute G in Eq. (27)
    G = (( 4*gammaval/(alpha-2)^2 + gammaval/(alpha-1) + 2*gammaval/(alpha-2) ) + cbar*gammaval*(1-epsilon^2)/(alpha-1) ) /tau /( (1-epsilon^2)*(1-(1+gammaval)*epsilon^2)*cbar - (1+2/(alpha-2))*gammaval);
    
    %Compute K* in Eq. (26)
    Kstar = (sqrt( (G*C0)^2 + C0*D1*cbar + C0*G*(C1+D0*cbar)) - G*C0)/(D1*cbar + G*(C1+D0*cbar));
    
    %Store the new optimized K
    Kiterations(iter) = Kstar;
    
    
    %Step 3: Update M using Theorem 4
    
    %Compute the parameters in Eqs. (31)-(36)
    a0 = gammaval*Kstar*(1-epsilon^2)/(tau*(alpha-1));
    a1 = (Kstar/tau)*(4*gammaval/(alpha-2)^2+gammaval/(alpha-1) + 2*gammaval/(alpha-2));
    a2 = (1-epsilon^2)*(1-(1+gammaval)*epsilon^2);
    a3 = (1+2/(alpha-2))*gammaval;
    a4 = C0+C1*Kstar;
    a5 = D0*Kstar+D1*Kstar^2;
    
    %Compute the solution in Eq. (29)
    cbar1 = (a1 + a3 + sqrt(a1*a3 + a1^2 + a1*a2*a4/a5 + a0*a3*a4/a5 - a0*a1*a4/a5 -a0^2*a3*a4/(a2*a5) + a0*a1*a3/a2 + a0*a3^2/a2 ) ) / (a2-a0);

    %Compute the solution in Eq. (30)
    cbar2 = gammaval*(1+4/(alpha-2)^2 + 1/(alpha-1) + 4/(alpha-2)) / ( (1-epsilon^2)*(1-(1+gammaval)*epsilon^2) - gammaval*(1-epsilon^2)/(alpha-1) );

    %Check if the solution in Eq. (29) satisfies the first constraint
    if ((4/(alpha-2)^2 + 1/(alpha-1) + 2/(alpha-2)) + cbar1*(1-epsilon^2)/(alpha-1))*gammaval/(cbar1*(1-epsilon^2)^2-((1 + 2/(alpha-2)) + cbar1*(1-epsilon^2)*epsilon^2)*gammaval) >= 1
        cbarStar = cbar1;
    else
        cbarStar = cbar2;
    end

    %Compute the resulting M value
    Miterations(iter) = cbarStar*Kstar;

    
    %Extract the current K and M values, which are used to compute the
    %current EE value
    k = Kiterations(iter);
    m = Miterations(iter);
    
    %Compute B1bar and B2bar from Eq. (18) and Eq. (19), respectively.
    B1bar = k*(4/(alpha-2)^2 + 1/(alpha-1) + 2/(alpha-2)) + m*(1-epsilon^2)/(alpha-1);
    B2bar = k*(1 + 2/(alpha-2)) + m*(1-epsilon^2)*epsilon^2;
    
    %Check if the two constraints in Eq. (23) are satisfied
    if B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval)>=1 && (k*B1bar*gammaval/(m*(1-epsilon^2)^2-B2bar*gammaval))<=tau %gammaval>=m*(1-epsilon^2)^2/(B1+B2) &&
        
        %Compute the objective function in Eq. (23) and store the EE
        ASE = k *(1-(B1bar*gammaval / (m*(1-epsilon^2)^2-B2bar*gammaval))*k/tau) * log2(1+ gammaval);
        AEC = C0 + C1*k + D0*m + D1*m*k + A*ASE;
        EEiterations(iter) = ASE/AEC;
        
        %Create a vector with the energy used by each part of the AEC
        powerParts = [C0 C1*k D0*m D1*m*k A*ASE];
        
    end

end



%Density of the lines that are used in the 3d plot to make it easier to
%see the shape
gridDensity = 10;



%Plot Figure 3 from the paper
figure(3);
hold on; box on; grid on;

surface(1:Kmax,1:Mmax,EEtheory/1e6,'EdgeColor','none');
colormap(autumn);

view([-17 32]);

xlabel('Number of UEs (K)')
ylabel('Number of BS antennas (M)');
zlabel('Energy efficiency [Mbit/Joule]');

%Plot lines on top of the 3d surface, to make it easier to see the shape
for m = [1 gridDensity:gridDensity:Mmax]
    plot3(1:Kmax,m*ones(1,Kmax),EEtheory(m,:)/1e6,'k-');
end

for k = [1 gridDensity:gridDensity:Kmax]
    plot3(k*ones(1,Mmax),1:Mmax,EEtheory(:,k)/1e6,'k-');
end

%Plot the progress of the alternating optimization algorithm
convergence = 4;
plot3(Kiterations(1:convergence-1),Miterations(1:convergence-1),EEiterations(1:convergence-1)/1e6,'ko-','LineWidth',1);

%Plot the optimal solution to the EE maximization problem
[EEmaxM,optM] = max(EEtheory,[],1);
[EEmax,optK] = max(EEmaxM);

plot3(optK,optM(optK),EEmax/1e6,'k-*','LineWidth',1);



%Plot Figure 4 from the paper
figure(4);
powerPartsNorm = powerParts/sum(powerParts);
pie(powerPartsNorm,{'C_0 (31%)','C_1 K (3%)','D_0 M (56%)','D_1 M K (9%)','A ASE (1%)'});
