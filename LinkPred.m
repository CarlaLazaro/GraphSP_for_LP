%% Link prediction methods
clear all
close all
format long

%% Select Data base
addpath(append(pwd,'/LoadData'))

% st='Lawyers'; G=LawyersData();

% db='USAir'; G=USAirData();

% db='NS';  G=NSData();

% db='Yeast'; G=YeastData();

% db='CEle'; G=CEleData();

% db='Power'; G=PowerData();

% db='Router'; G=RouterData(); 

% db='EColi';  G=EColiData();

% db='PB';  G=PBData();

%% Simulation parameters
numMeth = 9; %Number of methods to try
nR = 10; %Number of repetitions
unk_p = 10; %Unkown percentage of edges

%% Graph parameters
n_nodes = size(G.Nodes,1);
X = triu(ones(size(G.Nodes,1))); X = X-diag(diag(X));
[row,col] = find(X); st = [row col]; %All possible edges in a graph
n_edges = size(row,1); %Total number of possible edges
nM = fix(n_edges*unk_p/100); nO = n_edges-nM; %Number of Missing and Observed

A = adjacency(G);
A2 = A^2; A2 = full(A2);
%% Define Max Threshold and Threshold Step for each method

%-------(1)Common Neighbours-CN--------------------------------------------
thMaxCN = max(A2,[],'all')+1;

%------Random Set for threshold definition---------------------------------
thS_pos = randperm(n_edges,fix(n_edges*20/100))'; %Size=20% of original Gph
thS_st = st(thS_pos,:); %Nodes od the Set for Threshold Definition
thS_s = thS_st(:,1); thS_t = thS_st(:,2); 

%-------(2)Jaccard-J & (3)Adamic-Adar-AA-----------------------------------
thMaxJ = 0; thMaxAA = 0;
for j = 1:size(thS_s)
    
    score_CN = A2(thS_s(j),thS_t(j)); 
    s_neighbors = neighbors(G,thS_s(j));
    t_neighbors = neighbors(G,thS_t(j));
    
    U = union(s_neighbors,t_neighbors);
    size_U = size(U,1);
    score_J = score_CN/size_U;
    if score_J>thMaxJ
        thMaxJ = score_J;
    end
    
    isec = intersect(s_neighbors,t_neighbors);
    d = degree(G, isec);
    score_AA = ones(1,size(isec,1))*(1./log(d));
    if score_AA>thMaxAA
        thMaxAA = score_AA;
    end
end
%--------------(4)Katz-K---------------------------------------------------
b = 0.001;
I = eye(size(A));
scMAT_K = inv(I-b*A)-I;
thMaxK = max(scMAT_K,[],'all');
%-------------(5)Page-Rank-PR----------------------------------------------
alfa = 0.85; 
d = degree(G); d(d==0) = 1;
D_inv = diag(1./d);
P = A*D_inv;
P_aux = inv(eye(n_nodes)-alfa.*P);
scMAT_PR =(1-alfa).*P_aux;
thMaxPR = max(scMAT_PR,[],'all');
%----------(6)(7)(8)(9)FIR-Graph-Filters-Based-Methods---------------------
thMaxFIR = 1;  

%------------thMaxM & thStpM-----------------------------------------------
thMaxM = zeros(1,numMeth); thStpM = zeros(1,numMeth);

thMaxM(1) = thMaxCN; thMaxM(2) = thMaxJ;
thMaxM(3) = thMaxAA; thMaxM(4) = thMaxK;
thMaxM(5) = thMaxPR; thMaxM(6) = thMaxFIR;
thMaxM(7) = thMaxFIR; thMaxM(8) = thMaxFIR;
thMaxM(9) = thMaxFIR;

thStpM(1) = 1; thStpM(2) = 0.00001;
thStpM(3) = 0.001; thStpM(4) = 0.00000001;
thStpM(5) = 0.00001; thStpM(6) = 0.0001;
thStpM(7) = 0.0001; thStpM(8) = 0.0001; 
thStpM(9) = 0.0001; 

%% Run all methods
%----------------Vector-initializations------------------------------------
dimMeth = max(round(thMaxM./thStpM))+1;
Pd = zeros(dimMeth,numMeth); Pfa = Pd;
Pd_mean = Pd; Pfa_mean = Pd;

sqsumM_auc = zeros(1,numMeth);
aucM = zeros(1,numMeth);
aucM_mean = zeros(1,numMeth);

%---------------Start-nR-realizations--------------------------------------
rng(1) %Control the random set
for r=1:nR
    %% Missing and Observed Set
    M_set_pos = randperm(n_edges, fix(n_edges*unk_p/100))'; %Posit.MissSet
    st_M = st(M_set_pos,:); %Nodes Missing Set
    sM = st_M(:,1); tM = st_M(:,2); %Source and target nodes of MissSet
    ind_m = sub2ind(size(A),sM,tM);
    st_O = st; st_O(M_set_pos,:) = []; %Nodes Observed Set
    
    aM = full(A(ind_m)); %Value of the Missing Edges
    n_true = sum(aM); %Number of links to be detected
    n_false = nM-n_true;
    
    ind_o = sub2ind(size(A),st_O(:,1),st_O(:,2));
    ind_M = sub2ind(size(A),[sM;tM],[tM; sM]);
    Ao = A; Ao(ind_M) = 0; %Observed adjacency matrix
    Ao2 = Ao^2; Ao2 = full(Ao2);
    
    w = 0.35; %Weight for unknown links
    Aw = A; Aw(ind_M) = w; 
    
    unk_nod = [sM(aM==1);tM(aM==1)];
    dO = d; dO(unk_nod) = dO(unk_nod)-1; %Degree of the observed graph
    
    dO(dO==0) = 1; DinvO = diag(1./dO); %Observed degree matrix
    Po = Ao*DinvO; %Observed Random Walk 
    
   %% Calculate Scores from the observed graph 
    scoreCN = zeros(nM,1); scoreJ = zeros(nM,1); scoreAA = zeros(nM,1);
   for i = 1:nM %for each unknown edge
    %----(1)Common neighbours-CN-----------------
    scoreCN(i) = Ao2(sM(i),tM(i));
    %--------------------------------------------   

    s_neighbors = neighbors(G,sM(i));
    t_neighbors = neighbors(G,tM(i)); 
    %-----(2)Jaccard-J---------------------------      
    U = union(s_neighbors,t_neighbors); 
    scoreJ(i) = scoreCN(i)/size(U,1);
    %--------------------------------------------

    %-----(3)Adamic-Adar-AA----------------------
    isec = intersect(s_neighbors,t_neighbors);
    di = degree(G, isec);
    scAA = ones(1,size(isec,1))*(1./log(di));
    if isempty(scAA)==1
        scAA = 0;
    end
    scoreAA(i) = scAA;
    %-------------------------------------------  
   
   end
   
    %----(4)Katz-K-------------------------------
    sc_auxK = inv(I-b*Ao)-I;
    scoreK=sc_auxK(ind_m);
    %--------------------------------------------

    %----(5)Page Rank-PR-------------------------
    P_auxO = inv(eye(n_nodes)-alfa.*Po);
    sc_auxPR = (1-alfa).*P_auxO;
    scorePR = sc_auxPR(ind_m);
    %--------------------------------------------

   %----(6)MLE_A-LA------------------------------
    L = 8; a = A(ind_o);
    beta0 = zeros(L,1);
    it_max = 200;

    AL = zeros(n_nodes,n_nodes,L);
    AL(:,:,1) = Aw; mem = Aw;
    for j = 2:L
        mem = Aw*mem;
        norm_mem = norm(mem,1);
        AL(:,:,j) = mem./norm_mem;
    end
    
    [betaLA,Z_A] = MLE (beta0,a,AL,nO,L,ind_o,it_max);

    sM = betaLA(1,1).*ones(n_nodes);
    for i=2:L
        sM = sM+betaLA(i).*AL(:,:,i);
    end
    scLA = exp(sM)./(1+exp(sM));
    scoreLA = scLA(ind_m);
    
    %-----(7)MLE_P-LP----------------------------------
    PL = zeros(n_nodes,n_nodes,L);
    PL(:,:,1) = Po; mem = Po;
    for j = 2:L
        mem = Po*mem;
        PL(:,:,j) = mem;
    end
    
    [betaLP,Z_P] = MLE (beta0,a,PL,nO,L,ind_o,it_max); 
    
    sM = betaLP(1,1).*ones(n_nodes);
    for i = 2:L
        sM = sM+betaLP(i).*PL(:,:,i);
    end
    scLP = exp(sM)./(1+exp(sM));
    scoreLP = scLP(ind_m);
    
    %-----(8)MMSE_A-MSA--------------------------------
    betaMSA = (Z_A'*Z_A)\Z_A'*a;
    
    sM = betaMSA(1,1).*ones(n_nodes);
    for i = 2:L
        sM = sM+betaMSA(i).*AL(:,:,i);
    end
    scMSA = exp(sM)./(1+exp(sM));
    scoreMSA = scMSA(ind_m);
    
    %------(9)MMSE_P-MSP-------------------------------
    betaMSP = (Z_P'*Z_P)\Z_P'*a;
    
    sM = betaMSP(1,1).*ones(n_nodes);
    for i = 2:L
        sM = sM+betaMSP(i).*PL(:,:,i);
    end
    scMSP = exp(sM)./(1+exp(sM));
    scoreMSP = scMSP(ind_m);
        
%------Score vectors for each method-----------------------------
    scoreMeth = zeros (nM, numMeth);
    scoreMeth(:,1) = scoreCN; scoreMeth(:,2) = scoreJ;
    scoreMeth(:,3) = scoreAA; scoreMeth(:,4) = scoreK;
    scoreMeth(:,5) = scorePR; scoreMeth(:,6) = scoreLA;
    scoreMeth(:,7) = scoreLP; scoreMeth(:,8) = scoreMSA;
    scoreMeth(:,9) = scoreMSP;

    %% Obtain PD and PFA for each Method
    for k = 1:numMeth
        t = 1;
        for th = 0:thStpM(k):thMaxM(k)
            n_pd = 0; n_pfa = 0;
            for i = 1:nM
                if scoreMeth(i,k) >= th
                    if aM(i)==1
                       n_pd = n_pd+1; 
                    else
                       n_pfa = n_pfa+1;
                    end
                end
            end
          Pd(t,k) = n_pd/n_true;
          Pfa(t,k) = n_pfa/n_false;
          t = t+1;
        end
    end

    %% Mean values & AUC
    for k = 1:numMeth
      aucM(k) = trapz(fliplr(Pfa(:,k)'),fliplr(Pd(:,k)'));
      sqsumM_auc(k) = sqsumM_auc(k)+aucM(k)^2;
      aucM_mean(k) = 1/r*aucM(k)+(r-1)/r*aucM_mean(k); 
      Pd_mean(:,k) = 1/r*Pd(:,k)+(r-1)/r*Pd_mean(:,k);
      Pfa_mean(:,k) = 1/r*Pfa(:,k)+(r-1)/r*Pfa_mean(:,k);
    end

end
%% Standard deviation 
stdM = sqrt((sqsumM_auc/(nR-1))-((nR/(nR-1))*aucM_mean.^2));

%% Save results
filename = "Results_multMethod_"+db+".mat";
save(filename)

%% Plot
figure(1) 
for k = 1:numMeth
   plot(Pfa_mean(:,k),Pd_mean(:,k),'-','LineWidth',1)
   hold on
end
xlabel('False Positive Rate (1 - Specificity)') 
ylabel('True Positive Rate (Sensitivity)')
title(db) 
auc_lg=aucM_mean*100;
std_lg=stdM*100;
legend("Common Neighbours: " + round(auc_lg(1),2) + "\pm" + round(std_lg(1),2),...
    "Jaccard: " + round(auc_lg(2),2) + "\pm" + round(std_lg(2),2),...
    "Adamic-Adar: " +round(auc_lg(3),2) + "\pm" + round(std_lg(3),2),...
    "Katz: " + round(auc_lg(4),2) + "\pm" + round(std_lg(4),2) ,...
    "Page-Rank: " + round(auc_lg(5),2) + "\pm" + round(std_lg(5),2),...
    "MLE-A : "+ round(auc_lg(6),2) + "\pm" + round(std_lg(6),2),...
    "MLE-P: " + round(auc_lg(7),2) + "\pm" + round(std_lg(7),2),...
    "MMSE-A : " + round(auc_lg(8),2) + "\pm" + round(std_lg(8),2),...
    "MMSE-P : " + round(auc_lg(9),2)+ "\pm" + round(std_lg(9),2),...
    'location', 'southeast');