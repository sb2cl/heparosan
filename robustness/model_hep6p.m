% Naringerin Metabolic, Anthitetic controller and QdoR biosensor model.
% Updated 2022 Alejandro Vignoni, Yadira Boada

function [dxdt] = model_hep6p(t,x,p,pert,pertA)
MurAFlux = 1;
    if pert==1 %GlcUA pert
        FluxUA = 1*pertA;
        PertmurA = 1;
        NoContr = 0;
    elseif pert==2 %MurA pert
        FluxUA = 1;
        PertmurA = 1.5; 
        NoContr = 0;
    elseif pert==3 %MurA pert / without control
        FluxUA = 0.8;
        PertmurA = 1
        NoContr = 1;
        MurAFlux = 2.3384;
    else 
        FluxUA = 1;
        PertmurA = 1;
        NoContr = 0;
    end
Ncell = p.Ncell;
Size = p.Size; 

for k = 1:Ncell

%% A:QUORUM SENSING
%x1 = LuxR
dxdt(1,1)= p.pR(k)*p.CN(k)*p.kR(k)./(p.dm(k)+p.mu) - (p.dR(k)+p.mu)*x(1); 

%x2 = LuxI
dxdt(2,1)= p.pI(k)*p.CN(k)*p.kR(k)./(p.dm(k)+p.mu)-(p.dI(k)+p.mu)*x(2); 

%% B:OVEREXPRESSION

%x3 = PmHS2 AHL=x(m+19)
c3 =p.p3(k)*p.CN(k)*p.k3(k)./(p.dm(k)+p.mu);
Plux = p.alpha(k) + (1-p.alpha(k))*x(22)^2./( p.kdlux(k)*(p.kd2(k)*p.CN(k)./x(1))^2 + x(22)^2);
dxdt(3,1)=c3*Plux - (p.d3(k)+p.mu)*x(3);

%x4 = GlmM
c4 =p.p4(k)*p.CN(k)*p.k3(k)./(p.dm(k)+p.mu);
dxdt(4,1)=c4*Plux - (p.d4(k)+p.mu)*x(4);

%x5 = GalU
c5 =p.p5(k)*p.CN(k)*p.k3(k)./(p.dm(k)+p.mu);
dxdt(5,1)=c5*Plux - (p.d5(k)+p.mu)*x(5);

%x6 = KfiD
c6 =p.p6(k)*p.CN(k)*p.k3(k)./(p.dm(k)+p.mu);
dxdt(6,1)=c6*Plux - (p.d6(k)+p.mu)*x(6);

%% C: BIOSENSORS
%x7 = UxuR
c7 =p.p7(k)*p.CN(k)*p.k7(k)./(p.dm(k)+p.mu);
dxdt(7,1)=c7*Plux - (p.d7(k)+p.mu)*x(7);

%x8 = NagC
c8 =p.p8(k)*p.CN(k)*p.k7(k)./(p.dm(k)+p.mu);
dxdt(8,1)=c8*Plux - (p.d8(k)+p.mu)*x(8);

%% D: BIOCONTROLLER
%x9 = Sigma drive by UxuR promoter (GlcUA biosensor)
c9 = p.ps(k)*p.CN(k)*p.ks(k)./( p.dm(k)+p.mu );
PuxuRp = p.beta(k) +(1-p.beta(k))*x(15)^2./( 0.05*p.kdUxuR(k)*(p.kd9(k)*p.CN(k)./x(7))^2 + x(15)^2);
dxdt(9,1)=10*c9*PuxuRp - p.k_c(k)./p.kdc(k)*x(9)*x(10) + p.k_c(k)*x(11) - (p.ds(k)+p.mu)*x(9);
      
%x10 = AntiSigma drive by NagC promoter (GlcNAc6P biosensor)
c10 =p.pa(k)*p.CN(k)*p.ka(k)./(p.dm(k)+p.mu);
dxdt(10,1)=c10*(p.beta(k) +(1-p.beta(k))*x(17)^2./( p.kdNagC(k)*(p.kd10(k)*p.CN(k)./x(8))^2 + x(17)^2))-...
             p.k_c(k)./p.kdc(k)*x(9)*x(10) + p.k_c(k)*x(11) - (p.da(k)+p.mu)*x(10);
         
%x11 = Complex Sigma-Asigma
dxdt(11,1) = p.k_c(k)/p.kdc(k)*x(9)*x(10)-p.k_c(k)*x(11)-(p.dc(k)+p.mu)*x(11);

%% E: SILENCING
%x12 = siRNA - promoter P20
dxdt(12,1) = p.CN(k)*p.k12(k)*(p.alpha(k) + (1-p.alpha(k))*x(9)^2./( p.kd20*(p.kds(k)*p.CN(k))^2 + x(9)^2 ))-...
              p.kmT*x(12)*x(13)-(p.dm(k)+p.mu)*x(12);
          
%x13 = mRNAMurA
dxdt(13,1) = p.CN(k)*p.k13(k) - p.kmT*x(12)*x(13)- (p.dm(k)+p.mu)*x(13);

%x14 = MurA
dxdt(14,1) = PertmurA*(p.p14(k)*x(13)*(~NoContr)+MurAFlux*(NoContr))  - (p.d14(k)+p.mu)*x(14);
%p.p14(k)*x(13)
%% METABOLIC PATHWAY
V_H = p.catPm(k)*x(3)*x(15)*x(16)./( 1/(p.KmUA(k)*p.KmNAc(k))+x(15)/p.KmNAc(k)+x(16)/p.KmUA(k)+x(15)*x(16));
V_EP = p.catMurA(k)*x(14)*x(16)./(p.KmNAc(k)+x(16));
V_6P = p.catNagK(k)*x(20)*x(16)./(p.KmNAc6p(k)+x(16));
%x15=GlcUA   Concentrations in wild type [1.2e-4  2.67e-3] M
dxdt(15,1) =  p.GlcUA_factor.*1.1*1.45*FluxUA*0.0134*421610 - V_H - p.mu*x(15);   %Flux0 units: #molecules/min

%x16: GlcNAc  Concentrations in wild type [6.79e-3  1.26e-2] M
dxdt(16,1) = p.GlcNAc_factor.*1.2*1.01 * 0.08*421610-V_H - V_EP -V_6P- p.mu*x(16); %8-fold,then increase parameters 

%x17: GlcNAc6p
dxdt(17,1) = V_6P - p.mu*x(17);
          
%x18: HEPAROSAN
dxdt(18,1) = V_H - p.mu*x(18);
           
%x19: Peptidoglycans EP     
dxdt(19,1) = V_EP  - p.mu*x(19);

%x20 = NagK
c19 =p.p19(k)*p.CN(k)*p.k19(k)./(p.dm(k)+p.mu);
dxdt(20,1)=c19*Plux - (p.d8(k)+p.mu)*x(20);

%x21 = Number of cells
dxdt(21,1) = p.mu*x(21)*(1-x(21)/p.cellmax);

%x22 = AHLint
dxdt(22,1) = p.D*p.Vcell/p.Vext*x(Size*Ncell+1)-p.D*x(22)+ p.kA*x(2)-(p.dA(k)+ p.mu)*x(22);

end
%AHLext, AHL=x(20)
    dxdt((Size*Ncell+1),1) = -p.D*x(21)*p.Vcell/p.Vext*x(Size*Ncell+1) +...
                              p.D*x(21)*sum(x(22:Size:Size*Ncell)) - p.dAe*x(Size*Ncell+1);
end