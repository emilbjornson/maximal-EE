%This Matlab script can be used to generate Monte-Carlo realizations that
%are used to plot Figure 2 in the article:
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

%Pathloss exponent
alpha = 3.76;

%Define the coverage area (as a square with wrap-around). This dimension
%has no impact on  the Monte-Carlo results since we compute metrics that
%are scale invariant.
squareLength = 1000;

%Maximal number of Monte-Carlo setups
monteCarloSetups = 100;

%Average number of base stations per setup
averageNumberofBSs = 1000;

%Maximal number of users per cell
K = 100;

%Go through all setups
for n = 1:monteCarloSetups
    
    %Display simulation progress
    disp(['Realization ' num2str(n) ' out of ' num2str(monteCarloSetups)]);
    
    %Generate the number of BSs in the area
    nbrBSs = poissrnd(averageNumberofBSs);
    
    %If the number is zero, then make a new realization
    while nbrBSs == 0
        nbrBSs = poissrnd(averageNumberofBSs);
    end
    
    %Random BS locations with uniform distribution
    BSpositions = (rand(nbrBSs,1) + 1i*rand(nbrBSs,1)) * squareLength;
    
    %Compute a size of the distance from a BS where its users need to be
    %contained.
    distances = abs(BSpositions);
    sortedBSdistances = sort(distances,'ascend');
    maxDistance = 2*sortedBSdistances(min([10 nbrBSs]));
    
    %Compute alternative BS locations by using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrBSs 1]);
    
    
    %Prepare to put out users in the cells
    UEpositions = zeros(K,nbrBSs);
    perBS = zeros(nbrBSs,1);
    
    %Initiate matrices where first and order interference are computed
    interference1 = zeros(K,nbrBSs);
    interference2 = zeros(K,nbrBSs);
    
    %Go through all the cells
    for j = 1:nbrBSs
        
        
        %Put out K users in each cell
        while min(perBS(j))<K
            
            %One random location in the vicinity of BS j
            posX = rand(1,1)*maxDistance + real(BSpositions(j)) - maxDistance/2;
            posY = rand(1,1)*maxDistance + imag(BSpositions(j)) - maxDistance/2;
            position = mod(posX,squareLength) + 1i*mod(posY,squareLength);
            
            %Find closest BS (with wrap around)
            distancesWrapped = min(abs(BSpositionsWrapped - repmat(position,size(BSpositionsWrapped))),[],2);
            [~,index] = min(distancesWrapped);
            
            %Add the user if BS j is the closest one.
            if index == j
                perBS(j) = perBS(j) + 1;
                UEpositions(perBS(j),j) = position;
            end
            
        end
        
        %Compute the distance from the users in cell j to BS j
        distancesSquaredBSj = min(abs( repmat(UEpositions(:,j),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        for l = 1:nbrBSs
            
            %Compute the distance from the users in cell j to BS l
            distancesSquaredBSl = min(abs( repmat(UEpositions(:,j),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(l,:),[K 1]) ),[],2);
            
            %Compute inteference terms of the types that show up in Eq. (7)
            interference1(:,l) = interference1(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(alpha);
            interference2(:,l) = interference2(:,l) + (distancesSquaredBSj./distancesSquaredBSl).^(2*alpha);
            
        end
        
    end
    
    %Store the results so that the Monte-Carlo realizations can be reused
    results{n}.interference1 = interference1;
    results{n}.interference2 = interference2;
    results{n}.nbrBSs = nbrBSs;
    
end

%Save the results in file that can be used for performance analysis
save resultsMC results;
