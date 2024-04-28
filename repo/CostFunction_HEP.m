function [J,g]=CostFunction_HEP(param,plot_out)
%   Cost function to optimize the PDI
%   J is the cost
%   g are the constrains in the GlcNAc and GlcUA concentrations
%   Updated 28/04/2024 by Alejandro Vignoni and Yadira Boada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load heparosan_GP.mat
Variance = 0;
Stdeviation =0;
Ncell = 1;

%Initial Optical density
ODinitial = 0.01;
ODmax = 140;
total_time=48;    %Hours

%ODE
options = odeset('AbsTol',1e-8,'RelTol',1e-6);      % for ode function
step = 5;
growth_time=24;

%System size
species={'LuxR','LuxI','PmHS2','GlmM','GalU','KfiD','UxuR','NagCm',...
    'Sigma','Asigma','Complex','RNAi','mRNA_MurA','MurA',...
    'GlcUA','GlcNAc','GlcNAc6p','Heparosan','Peptidoglycans','NagK',...
    'OD','AHL','AHLe'};
States = length(species);

%Experiment
labels = {'0'};    %Induction in the lab [nM]
induction = [0];  %Induction in the lab [nM]

%Fermentation & initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vext_0 = 4e-3; %L
Vext_reactor = 7; %L
Cellinitial = 1;

%% 0 Null initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = parameters(Ncell,Stdeviation,ODmax,Vext_0,States);
p.k12 = param(3);
p.kmT = 0.133;
p.GlcNAc_factor = param(1);
p.GlcUA_factor = param(2);

Initial = zeros(1,States);
Initial(length(Initial)-2)=1; %ini conditions cells
tfin = 60*16;     %simulation time
tspan = 0:step:tfin-step;
pert =0;
[t0,x0] = ode23t(@(t,x) model_hep6p(t,x,p,pert),tspan, Initial, options);
%% 1 GROWTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = parameters(Ncell,Stdeviation,ODmax,Vext_reactor,States);

p.k12 = param(3); %4.96;
p.kmT = 0.133;
p.GlcNAc_factor = param(1);
p.GlcUA_factor = param(2);

Initial(end-1:end)=0; %AHL & AHLext
ahle0 = induction.*p.ahl_to_molec; %AHLe vector
tfin = 60*growth_time;  % Time (min)
tspan = 0:step:tfin-step;

for j=1:length(induction)
    Initial = x0(end,:);
    Initial(length(Initial)-2) = ODinitial*p.Vext*p.OD_to_cells;
    Initial(end)= ahle0(j); %AHLext

    [t1,x1] = ode23t(@(t,x) model_hep6p(t,x,p,pert),tspan, Initial, options);
    for k=1:length(species)
        if k==21 %OD
            data = x1(:,k)./(p.OD_to_cells*p.Vext); %Cells to OD
            data2=data;
        elseif k==length(species)
            data = x1(:,k);
            data2 = x1(:,k).*p.ahlmolec_to_nM;  %molecules to nM
        else
            data = x1(:,k);
            data2 = data.*p.molec_to_nM;
        end
        T1molec(:,k) =table(data);
        T1conc(:,k) =table(data2);
    end
end
out.time=t1./60;

%% 2 PERTURBATION GlcUA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ahle0 = induction.*p.ahl_to_molec; %AHLe vector
tfin = t1(end)+60*(total_time-growth_time); % Tiempo de simulacion (min)
tspan = t1(end):step:tfin-step;
pert =1;

for j=1:length(induction)
    Initial=x1(end,:);
    Initial(15)=Initial(15)*1; %GlcUA decreases 40%

    [t1,x1] = ode23t(@(t,x) model_hep6p(t,x,p,pert),tspan, Initial, options);
    for k=1:length(species)
        if k==21 %OD
            data = x1(:,k)./(p.OD_to_cells*p.Vext); %Cells to OD
            data2=data;
        elseif k==length(species)
            data = x1(:,k);
            data2 = x1(:,k).*p.ahlmolec_to_nM;  %molecules to nM
        else
            data = x1(:,k);
            data2 = data.*p.molec_to_nM;
        end
        T2molec(:,k) =table(data);
        T2conc(:,k) =table(data2);
    end
    Tmolec = [T1molec;T2molec];
    Tconc = [T1conc;T2conc];
    Tmolec.Properties.VariableNames = species;
    Tconc.Properties.VariableNames = species;
    T1molec.Properties.VariableNames = species;
    T1conc.Properties.VariableNames = species;
    field=['AHL_' num2str(induction(j))];
    out.(field).molec=table2struct(Tmolec,"ToScalar",true);
    out.(field).conc=table2struct(Tconc,"ToScalar",true);
end
out.time=[out.time; t1./60];
%% Heparosan polydispersion index
%% Caluclation of the PDI index and molecular weight of Heparosan

GlcNAc_actual_concentration = [T1conc.GlcNAc(end); Tconc.GlcNAc(end)]*1e-6; %mM
GlcUA_actual_concentration = [T1conc.GlcUA(end); Tconc.GlcUA(end)]*1e-6; %mM
%Tconc.GlcNAc(end)*1e-6;
if plot_out == 1
    figure
    c = categorical({'1. Before perturbation' '2. After perturbation'});
    bar1 = bar(c,[GlcUA_actual_concentration GlcNAc_actual_concentration]);
    set(bar1(1),'DisplayName','GlcUA');
    set(bar1(2),'DisplayName','GlcNAc');
    legend;
    ylabel('Substrate concentration [mM]');
end
PDI_estimated = predict(gprMdl_pdi, [GlcNAc_actual_concentration GlcUA_actual_concentration]);
Hep_Mw_estimated = predict(gprMdl_Mw, [GlcNAc_actual_concentration GlcUA_actual_concentration]);

Restriction_Mw  = Hep_Mw_estimated(1);
Restriction_Mwp  = Hep_Mw_estimated(2); 
Restriction_NAc = GlcNAc_actual_concentration(1);
Restriction_UA  = GlcUA_actual_concentration(1);
Restriction_NAcp = GlcNAc_actual_concentration(2); 
Restriction_UAp  = GlcUA_actual_concentration(2); 

J = (PDI_estimated(1) + PDI_estimated(2))/2; 

g(1) = Restriction_Mw;
g(2) = Restriction_Mwp;
g(3) = Restriction_NAc;
g(4) = Restriction_UA;
g(5) = Restriction_NAcp;
g(6) = Restriction_UAp;
g(7) = Restriction_UA- Restriction_NAc;
g(8) = Restriction_UAp -Restriction_NAcp;

%% Plots

if plot_out == 1
    J = [PDI_estimated(1)  PDI_estimated(2)]; %
    total_rows=3;
    total_cols=2;
    total_fig= round( length(species)/(total_rows*total_cols) );
    j=0;
    counter=0;
    counter=j + counter;

    for f = 1:total_fig
        if f==total_fig
            plot_number= length(species)-counter;
        else
            plot_number=total_rows*total_cols;
        end
        figure;   %figure('Color',[1 1 1]);
        for j=1:plot_number
            hs=subplot(total_rows,total_cols,j);
            for k=1:length(induction)
                field=['AHL_' num2str(induction(k))];
                field2=char(species(j+counter));
                plot(out.time,out.(field).conc.(field2),'LineWidth',3); %molecules or concentration
                hold on;
            end
            title(species(j+counter)); xlabel('Time (h)'); ylabel('nM');
            %legend(labels);
            grid on;
            hold off;
        end
        counter=j + counter;
    end

    %% Promoters figures
    Plux = p.alpha+(1-p.alpha)*Tmolec.AHL.^2./( p.kdlux*(p.kd2*p.CN./Tmolec.LuxR).^2 + Tmolec.AHL.^2);
    PUxuRp=p.beta +(1-p.beta)*Tmolec.GlcUA.^2./( p.kdUxuR*(p.kd9*p.CN./Tmolec.UxuR).^2 +Tmolec.GlcUA.^2);
    PNagBp=p.beta+(1-p.beta)*Tmolec.GlcNAc6p.^2./( p.kdNagC*(p.kd10*p.CN./Tmolec.NagCm).^2 +Tmolec.GlcNAc6p.^2);
    P20=p.alpha+ (1-p.alpha)*Tmolec.Sigma.^2./(p.kd20*(p.kds*p.CN).^2 +Tmolec.Sigma.^2);

    figure;
    subplot(221)
    plot(Tconc.AHL,Plux,'LineWidth',3);
    xlabel('AHL (nM)'); ylabel('Plux promoter');

    subplot(222)
    plot(Tconc.GlcUA,PUxuRp,'LineWidth',3);
    xlabel('GlcUA (nM)'); ylabel('UxuRp promoter');

    subplot(223)
    plot(Tconc.GlcNAc6p,PNagBp,'LineWidth',3);
    xlabel('GlcNAc6p (nM)'); ylabel('NagBp promoter');

    subplot(224)
    plot(Tconc.Sigma,P20,'LineWidth',3);
    xlabel('Sigma (nM)'); ylabel('P20 promoter');


    [GlcNAc_I, GlcUA_I] = meshgrid(0:0.5:6, 0:0.5:6);

    % Predict PDI values using the Gaussian Process Regression Model
    [PDI_pred,~,PDI_int] = predict(gprMdl_pdi, [GlcNAc_I(:), GlcUA_I(:)]);

    %Mw_pred = predict(gprMdl_Mw, [GlcNAc_I(:), GlcUA_I(:)]);

    % Reshape the predicted PDI values into the same size as the grid
    PDI_predr = reshape(PDI_pred, size(GlcNAc_I));
    PDI_int_u = reshape(PDI_int(:,1), size(GlcNAc_I));
    PDI_int_l = reshape(PDI_int(:,2), size(GlcNAc_I));

    % Plot the surface of GlcNAc, GlcUA, and PDI
    figure
    colormap(winter);
    surf(GlcNAc_I, GlcUA_I, PDI_predr)
    hold on
    % plot3(X(:,1),X(:,2), Z_pdi,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    %     'MarkerSize',30,...
    %     'Marker','.',...
    %     'LineStyle','none',...
    %     'Color',[0 0 0]);
    plot3(GlcNAc_actual_concentration,GlcUA_actual_concentration,PDI_estimated,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',30,...
        'Marker','.',...
        'LineStyle','none',...
        'Color',[1 0 0]);
    surf(GlcNAc_I, GlcUA_I, PDI_int_l,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    surf(GlcNAc_I, GlcUA_I, PDI_int_u,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    xlabel('N-Acetyl Glucosamine [mM]')
    ylabel('Glucoronic Acid [mM]')
    zlabel('Polydispersion Index')
      
    % Predict PDI values using the Gaussian Process Regression Model
    % Int is the 95 percent confidence interval
    [Mw_pred,Mw_sd,Mw_int] = predict(gprMdl_Mw, [GlcNAc_I(:), GlcUA_I(:)]);

    % Reshape the predicted PDI values into the same size as the grid
    Mw_predr = reshape(Mw_pred, size(GlcNAc_I));
    Mw_int_u = reshape(Mw_int(:,1), size(GlcNAc_I));
    Mw_int_l = reshape(Mw_int(:,2), size(GlcNAc_I));


    %Mw_predr = reshape(Mw_pred, size(GlcNAc_I));
    figure
    colormap(parula);
    % Plot the surface of GlcNAc, GlcUA, and PDI
    surfc(GlcNAc_I, GlcUA_I, Mw_predr);
    hold on
    % plot3(X(:,1),X(:,2), Z_Mw,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    %     'MarkerSize',30,...
    %     'Marker','.',...
    %     'LineStyle','none',...
    %     'Color',[0 0 0]);
    plot3(GlcNAc_actual_concentration,GlcUA_actual_concentration,Hep_Mw_estimated,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
        'MarkerSize',30,...
        'Marker','.',...
        'LineStyle','none',...
        'Color',[1 0 0]);
    %surf(GlcNAc_I, GlcUA_I, Mw_int_l,'FaceAlpha',0.2,'EdgeAlpha',0.2);
    %surf(GlcNAc_I, GlcUA_I, Mw_int_u,'FaceAlpha',0.2,'EdgeAlpha',0.2);
    xlabel('N-Acetyl Glucosamine [mM]')
    ylabel('Glucoronic Acid [mM]')
    zlabel('Heparosan Molecular weight (Mw)')
    view([74 34])


end

end








