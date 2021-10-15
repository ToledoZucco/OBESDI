%% Modified 27 May 2021

function Save_Figures_2plot(z3D,zd3D,t3D,w,wh,k_vid,figureDirectory,filename,saveFigures,LegendSim1,LegendSim2,Color1,Color2)

if isfolder(figureDirectory) == 1 % If the file figureDirectory exists, delete it with all the files insisde
    rmdir(figureDirectory,'s');
    pause(0.1);
end

mkdir(figureDirectory); %Othercase, create the file figureDirectory
pause(0.1);


font=50; lw=3; ms = 20;
x0screen=100;y0screen=100;WidthScreen=1000;HeightScreen=750;

figure;
set(gcf,'units','points','position',[x0screen,y0screen,WidthScreen,HeightScreen])
hold on
p_w = plot(z3D(:,1),w(:,1),'color',Color1,'LineWidth',lw,'MarkerSize',ms);
p_wh = plot(zd3D(:,1),wh(:,1),'.-','color',Color2,'LineWidth',lw,'MarkerSize',ms);
%p_wh = plot(zd3D(:,1),wh(:,1),'.-','LineWidth',lw,'MarkerSize',ms);
%title({'$t=$'},'Interpreter','latex','FontSize',font)
t_w = title(sprintf('$t=%.2f$ $[s]$', t3D(1,1)),'Interpreter','latex','FontSize',font);
legend({['$w(\zeta,t)$ ',LegendSim1],['$\hat{w}(\zeta,t)$ ',LegendSim2]},'Location','northoutside','Interpreter','latex','FontSize',font)
ylabel({'$[m]$'},'Interpreter','latex','FontSize',font)
xlabel({'$\zeta$ $[m]$'},'Interpreter','latex','FontSize',font)
%grid on
xlim([0,1])
ylim([-1,1])
set(gca,'FontSize',font);

if saveFigures == true
    
    for i = 1:length(k_vid)

        k = k_vid(i);

%         set(p_wOL,'XData',z3D(:,k),'YData',wOL(:,k))
        set(p_w,'XData',z3D(:,k),'YData',w(:,k))
        set(p_wh,'XData',zd3D(:,k),'YData',wh(:,k))
        set(t_w,'String',sprintf('$t=%.2f$ $[s]$', t3D(1,k)))


%         print(gcf,'-dpng','-r100','-loose',[figureDirectory,'/',filename,num2str(i),'.png']);
        print(gcf,'-dpng','-r100',[figureDirectory,'/',filename,num2str(i),'.png']);


        pause(0.1);
    end
else
    
    for i = 1:length(k_vid)

        k = k_vid(i);

%         set(p_wOL,'XData',z3D(:,k),'YData',wOL(:,k))
        set(p_w,'XData',z3D(:,k),'YData',w(:,k))
        set(p_wh,'XData',zd3D(:,k),'YData',wh(:,k))
        set(t_w,'String',sprintf('$t=%.2f$ $[s]$', t3D(1,k)))

        pause(0.2);

    end

end


end