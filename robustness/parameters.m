%   These are the parameters of the Anthitetic circuit Sigma-Anti-sigma and
%   LuxR protein.
%   Update 05/08/2019 by Yadira Boada

function [p] = parameters(Ncell,sd, ODmax, Vext,States)
%%%%%%%%%%%%%%%%%%%%%%%%  General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    p.Size =  States-1;
    p.Ncell = Ncell;
    p.doubling = 80;          % doubling time [min]
    p.mu = 0.88/60;%log(2)/p.doubling;        % growth rate [1/min]
    p.nA = 6.023e23;                          % Avogadro's number: # particles/mol
    p.Vcell =  1.1e-15;                   % typical volume of E. coli (liters). Source: Bionumbers
    p.Vext = Vext;                   %culture medium volume.  microfluidic device = 1e-9 liters
                                     % From Solvej Siedler, Novel biosensors based on flavonoid
    p.OD_to_cells = 8e11;               % cells/liter for OD=1, Agilent, E. coli Cell Culture
    p.cellmax = ODmax*p.Vext*p.OD_to_cells;     % maximum number of cells
    p.molec_to_nM = 1/(p.Vcell*p.nA)*1e9; % uM Conversion factor from number of particles to concentration (1 nMolar=nanomols/liter)
    p.ahlmolec_to_nM = 1/(p.Vext*p.nA)*1e9;        %nM in Molarity
    p.ahl_to_molec = p.Vext*p.nA*1e-9; %%in nM       %nM in Molarity
    %p.muVector=[];
    
    for k = 1:Ncell
    p.w(k) = 0;   % PERTURBATION
    %p.pN_luxI = 17;                   % plasmid number  pBR322 (15-20 copies/cell)
    p.CN(k) = 10 + sd*randn(1);   % plasmid number  pACYC184 (10 copies/cell)
    
    %Sigma
    p.dms(k) = log(2)/3 + sd*randn(1);             % degradation rate mRNA [1/min]
    p.dm = p.dms;
    p.ds(k) = 0.0003 + sd*randn(1);                 % degradation rate  [1/min]. Rapid degradation: RR Burgess, doi: 10.1006/rwgn.2001.1192
    p.ps(k) =  2*0.45 + sd*randn(1); %b*dms;    % translation rate 0.09-0.13 [1/min] from our rates calculator
    p.ks(k) = 1.1*1.5*3.2 + sd*randn(1);           % transcription rate [0.18653 - 1.8653] from our rates calculator    
    p.kd20(k) = 2*600*100 + sd*randn(1);        % dissociation cte to promoter [molecules], Annunziata 2017
    
    %Anti-sigma
    p.alpha(k) = 0.05 + sd*randn(1);             % basal expresion pLux
    p.beta(k) = 0.01 + sd*randn(1);             % basal expresion pUxuRp
    p.dma(k) = log(2)/3 + sd*randn(1);           % mRNA degradation rate [1/min]
    p.da(k) = 0.0003 + sd*randn(1);               % protein degradation rate [1/min]
    p.pa(k) =  0.18 + sd*randn(1);               % translation rate 0.15-0.2  [1/min] [3.9648 - 7.9295]
    p.ka(k) = 1.5*2.06 + sd*randn(1);              % transcription rate 0.18653 - 1.8653 [1/min] from our rates calculator
        
    %SigmaComplex
    p.kdc = 0.01 + sd*randn(1);                  % dissociation constant [molecules] Annunziata 2017 an orthogonal multi-input
    p.k_c = 1.8e-3 + sd*randn(1);                % [1/min] Annunziata 2017 an orthogonal multi-input
    p.kd9 = 10 + sd*randn(1); %les than kd2
    %p.kc = p.k_c/p.kdc + sd*randn(1);            % binding rate sigma to anti-sigma [min^-1 molecules^-1]
    p.dc = 0.001 + sd*randn(1);                   % degradation rate [1/min] Annunziata 2017 an orthogonal multi-input
    
    %LuxR
    p.beta1(k) = 0.01 + sd*randn(1);                % basal expresion p20
    p.dmR(k) = log(2)/3 + sd*randn(1);              % mRNA degradation rate  [1/min]
    p.dR(k) = 0.02 + sd*randn(1);                   % degradation rate [1/min]
    p.pR(k) = 0.05 + sd*randn(1);                  % translation rate [1/min] from our rates calculator
    p.kR(k) = 1 + sd*randn(1);                   % transcription rate [0.39164 - 3.9164] from our rates calculator
    
    %LuxI
    p.dI(k) = 0.02 + sd*randn(1);                   % degradation rate [1/min]
    p.pI(k) = 0.075 + sd*randn(1);                  % translation rate  [1/min] [0.05228 - 0.073192] from our rates calculator
    p.kI(k) = 2.5 + sd*randn(1);                   % transcription rate [1/min] [0.22689 - 2.2689] from our rates calculator
    p.kA = 0.04;         %Synthesis rate of AHL
    
    % Monomer LuxR.AHL
    p.kd1(k) = 100 + sd*randn(1);                   % dissociation constant of R to A [nM], Urbanowski etal. 2004
    p.k_1(k) = 10 + sd*randn(1);                   % unbinding rate LuxR to AHL [1/min]
    p.k1(k) = p.k_1(k)/p.kd1(k) + sd*randn(1);            % binding rate LuxR to AHL [1/min]
    p.dRA(k) = log(2)/5 + sd*randn(1);              % degradation rate of (LuxR.A) [1/min]. Buchler et al. 2004 Monomer half-life is just few minutes.
    
    % Dimer (R.A)2
    p.kd2(k) = 0.7*20 + sd*randn(1);                   %dissociation cte (LuxR.A) to (LuxR.A) [nM], Buchler et al. 2003
    %Koren, R. & Hammes, G. G. (1976) Biochemistry 15, 1165�1171.
    %Northrup, S. H. & Erickson, H. P. (1992) Proc. Natl. Acad. Sci. USA 89, 3338�3342.
    p.k_2(k) = 1 + sd*randn(1);                     % dissociation rate [1/min]
    p.k2(k) = p.k_2(k)/p.kd2(k) + sd*randn(1);               % binding rate LuxR to AHL [1/min]
    p.kdlux(k) = 300 + sd*randn(1);           % 600 nM dissociation cte (LuxR.A)2 to promoter [nM], Bucler et al [1 1000]nM
    
    %Sigma dimer
    p.kds(k) = 0.05*200 + sd*randn(1);
    
    %GFP
    %p.pg(k) =  0.08 + sd*randn(1); %0.06-0.09 Re-written in COST.m   CONSTITUTIVE translation rate   
    %p.kg(k) = 1.38 + sd*randn(1);   %0.25-2.5 1/min
    %p.dmg(k)= log(2)/3 + sd*randn(1); 
    %p.dg(k) = log(2)/420 + sd*randn(1);   %7 hours, or 40 min with LVA. New Unstable Variants of Green Fluorescent Protein for Studies of Transient Gene Expression in Bacteria. 
   
    %mRFP
    %p.pf(k) =  0.08 + sd*randn(1); % 0.06-0.09 40min of maturation time BIONUMBERS
    %p.kf(k) = 0.82 + sd*randn(1);
    %p.dmf(k)= p.dmg(k); 
    %p.df(k) = p.dg(k); 

    %AHL and AHLe
    p.D = 2;        %kinetic rate of AHL external transport [1/min] across the cell membrane, calculated
    
    p.dA(k) = 0.0004 + sd*randn(1);    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    p.dAe(k) = 0.0000481 + sd*randn(1);          % Horswill et al., 2007  %0.0164, Degradation rate for AHL. From Kaufmann etal. 2005. Similar to You etal. Nature, 2004
    %[0.05 0.03 min^-1]Degradation from Bionumbers online
    %0.000282 Degradation rate for external AHL. Fitted using half-life of 180 minutes, from Englmann etal. 2007
    % In Kauffmann & Sartorio, 2005 they use 0.0018

    %PmHS2
    p.d3(k) = 0.002 + sd*randn(1);
    p.p3(k) = 0.02 + sd*randn(1);    %[0.010589 - 0.014824]
    p.k3(k)=0.4;     %0.22 1/min [0.039173 - 0.39173] from rates calculator
    
    %GlmM
    p.d3(k) = 0.002 + sd*randn(1);
    p.p4(k) =  0.02+ sd*randn(1);    %[0.021416 - 0.03]
    p.d4(k) =  p.d3(k)+ sd*randn(1);    %[0.021416 - 0.03]
    
    %GalU
    p.d5(k) = 0.002 + sd*randn(1);
    p.p5(k) = 0.05 + sd*randn(1);    %[0.041872 - 0.05862]
    
    %KfiD
    p.d6(k) = 0.002 + sd*randn(1);
    p.p6(k) = 0.04 + sd*randn(1);    %[0.027455 - 0.038436]
    
    %UxUR 
    p.d7(k) = 0.002 + sd*randn(1);
    p.p7(k) = 0.06 + sd*randn(1);    %[0.056667 - 0.079333]
    p.k7(k) = 0.6;     %Ref  1/min [0.07, 0.10127 - 0.68,1.0127] from rates calculator
    p.kd9(k) = 80*30 + sd*randn(1);
    
    %NagCm
    p.d8(k) = 0.002 + sd*randn(1);
    p.p8(k) = 0.03 + sd*randn(1);    %[0.025297 - 0.035416]
    p.k8(k) = p.k7(k);     % 1/min [0.10127 - 1.0127] from rates calculator
    p.kd10(k) = 60*20 + sd*randn(1);
    
    %Silencing & RNA of MurA
    p.kmT = 1.33e-4*1000;            %1/ (molec.min) M A C Groenenboom, Athanasius F M Mar´ee, et al. The rna silencing pathway: the bits and pieces that matter. PLoS Computational Biology, 1(2):e21, 2005.
    p.k12 = 2*2* 1.24;    %mean=6.82, 1.2414 - 12.4138   from calculator
    %p.k13(k) =  0.05;     %mean=0.85 [0.15418 - 1.5418] 1/min
    p.d14(k) = 0.001 + sd*randn(1);
    %p.p14(k) = 0.05 + sd*randn(1);    %[0.023894 - 0.033451]
    
    
    %Nuevas constantes:
    p.kdUxuR = 15*20*0.4*71101; %molecules  71101.686;
    p.kdNagC = 0.0005*10*150.937; %molecules
    p.k13 = 4.286; %1/min
    p.p14 = 0.09*1.718; %1/min
    p.catPm = 240; %1/min
    p.catMurA = 228; %1/min
    p.KmUA = 6037462; %molecules
    p.KmNAc = 150936; %molecules
    p.wUA =194.139; %g/mol
    p.wNAc = 221.210; %g/mol
    
    %For N-acetyl 6P, 
    p.k19 = 0.1; %1/min
    p.p19 = 0.08;       % [0.034229 - 0.047921] from our rates 1/min  
    
    %Enzyme NagK
    p.catNagK = 12120/15000; %1/min
    p.KmNAc6p = 147917.830; %molecules
    p.wNAc6p = 301.188; %g/mol
    %p.catNagk = 250; %g/mol
            
%     p.catKfiD =    %1 - 6.6 mM  EC1.1.1.22


end