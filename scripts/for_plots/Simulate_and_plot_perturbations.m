load('results_cL1_10_cL2_10_cU1_20_cU2_20_date_20240331_154057.mat')
p = Results.xbest
[a, b]=CostFunction_HEP(p,0.6,3,0)
fig_handle = gcf;
pert = [ 0.8 1 1.2 1.4]
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),3,fig_handle);
end
clim([30 110]);

%%
[a, b]=CostFunction_HEP(p,0.6,2,0);
fig_handle = gcf;
for i=1:numel(pert)
    [ax, bx]=CostFunction_HEP(p,pert(i),2,fig_handle);
    b = [b; bx];
    a = [a; ax];
end

tableOutput = table(a,b(:,1:2),b(:,[3 5]),b(:,[4 6]),b(:,[7 8]),'VariableNames',{'PDI','MW','UDP-GlcNAc','UDP-GlcUA','Titer'});


%%
[a, b]=CostFunction_HEP(p,1,41,0)
fig_handle = gcf;
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),41,fig_handle);
end

%%
[a, b]=CostFunction_HEP(p,1,43,0)
fig_handle = gcf;
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),43,fig_handle);
end

%%
[a, b]=CostFunction_HEP(p,1,42,0)
fig_handle = gcf;
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),42,fig_handle);
end

%%
[a, b]=CostFunction_HEP(p,1,44,0)
fig_handle = gcf;
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),44,fig_handle);
end
%%


[a, b]=CostFunction_HEP(p,0.6,5,0)
fig_handle = gcf;
for i=1:numel(pert)
    [a, b]=CostFunction_HEP(p,pert(i),5,fig_handle);
end
hlg=legend('-40%','-20%','Typical conditions','20%','40%','Orientation','horizontal')
hlg.Title.String='UDP-GlcUA availability'

