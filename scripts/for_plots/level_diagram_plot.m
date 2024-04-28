clear
load tablaOutput.mat
load CM50.mat
%% 

figure
colormap(Colormap_MW50)
subplot(1,3,1)
plot(tableRandom.MW,tableRandom.PDI,'.','MarkerFaceColor',[0.1 0.1 0.1],'MarkerSize',10);
hold on
for i=1:5
    plot(tableOutput.MW(i,2),tableOutput.PDI(i,2),'.','MarkerSize',20);
end
    xlabel('Heparosan Mw (kDA)');
ylabel('Heparosan polydispersity index (PDI)');
ylim([0.95 1.65])
hold off
grid on
subplot(1,3,2)
plot(tableRandom.("UDP-GlcNAc"),tableRandom.PDI,'.','MarkerFaceColor',[0.1 0.1 0.1],'MarkerSize',10);
hold on
for i=1:5
plot(tableOutput.("UDP-GlcNAc")(i,2),tableOutput.PDI(i,2),'.','MarkerFaceColor',Colormap_MW50(i,:),'MarkerSize',20);
end
xlabel('UDP-GlcNAc Concentration (mM)');
ylim([0.95 1.65])
hold off
grid on
subplot(1,3,3)
plot(tableRandom.("UDP-GlcUA"),tableRandom.PDI,'.','MarkerFaceColor',[0.1 0.1 0.1],'MarkerSize',10);
hold on
for i=1:5
plot(tableOutput.("UDP-GlcUA")(i,2),tableOutput.PDI(i,2),'.','MarkerFaceColor',Colormap_MW50(i,:),'MarkerSize',20);
end
grid on
xlabel('UDP-GlcUA Concentration (mM)');
ylim([0.95 1.65])
