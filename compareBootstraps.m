function compareBootstraps(twoSigma,twoSigma_v2,twoSigma_v3,conf95,conf95_v2,conf95_v3, grouped)

if grouped==false
    len = size(twoSigma);
    for i=1:len(1)
        h = figure;
        if i < 3
            newTitle = sprintf('Linear - Coefficient %i',i);
            title(newTitle);
        elseif i < 6
            newTitle = sprintf('2nd Order Polynomial - Coefficient %i',i-2);
            title(newTitle);
        else
            newTitle = sprintf('3rd Order Polynomial - Coefficient %i',i-5);
            title(newTitle);
        end

        maxY = ylim;
        hold on
        plot(twoSigma(1,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
        hold on
        plot(conf95(1,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

        hold on
        plot(twoSigma_v2(1,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
        hold on
        plot(conf95_v2(1,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

        hold on
        plot(twoSigma_v3(1,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
        hold on
        plot(conf95_v3(1,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

        set(gca,'YTickLabel',[])
        ylim([0.95,1.3]);
        legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
               'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
               'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
               'location','northoutside');
        path = sprintf('/plots/bootstrapID%i.png',i);
        a = [pwd path];
        saveas(h,a);
        drawnow;
    end

    
    
    
    
    
    
    
else
    % linear 
    figure;
    subplot(1,2,1)
    title('Linear - Coefficient 1');
    maxY = ylim;
    hold on
    plot(twoSigma(1,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(1,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(1,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(1,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(1,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(1,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95,1.3]);
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');


    subplot(1,2,2)
    title('Linear - Coefficient 2');
    maxY = ylim;

    hold on
    plot(twoSigma(2,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(2,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(2,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(2,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(2,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(2,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');




    % 2nd Order
    figure;
    subplot(2,2,1)
    title('2nd Order Polynomial - Coefficient 1');
    maxY = ylim;
    hold on
    plot(twoSigma(3,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(3,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(3,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(3,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(3,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(3,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95,1.3]);
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');

    subplot(2,2,2)
    title('2nd Order Polynomial - Coefficient 2');
    maxY = ylim;

    hold on
    plot(twoSigma(4,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(4,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(4,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(4,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(4,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(4,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');


    subplot(2,2,3)
    title('2nd Order Polynomial - Coefficient 3');
    maxY = ylim;

    hold on
    plot(twoSigma(5,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(5,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(5,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(5,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(5,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(5,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');






    % 3rd Order
    figure;
    subplot(2,2,1)
    title('3rd Order Polynomial - Coefficient 1');
    maxY = ylim;
    hold on
    plot(twoSigma(6,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(6,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(6,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(6,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(6,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(6,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95,1.3]);
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');

    subplot(2,2,2)
    title('3rd Order Polynomial - Coefficient 2');
    maxY = ylim;

    hold on
    plot(twoSigma(7,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(7,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(7,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(7,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(7,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(7,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');


    subplot(2,2,3)
    title('3rd Order Polynomial - Coefficient 3');
    maxY = ylim;

    hold on
    plot(twoSigma(8,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(8,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(8,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(8,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(8,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(8,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');
       
       
    subplot(2,2,4)
    title('3rd Order Polynomial - Coefficient 4');
    maxY = ylim;

    hold on
    plot(twoSigma(9,:), 1.25*[maxY(2), maxY(2)], 'r--x','LineWidth',1.5);
    hold on
    plot(conf95(9,:), 1.2*[maxY(2), maxY(2)],'r-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v2(9,:), 1.15*[maxY(2), maxY(2)], 'g--x','LineWidth',1.5);
    hold on
    plot(conf95_v2(9,:), 1.1*[maxY(2), maxY(2)],'g-o','LineWidth',1.5);

    hold on
    plot(twoSigma_v3(9,:), 1.05*[maxY(2), maxY(2)], 'b--x','LineWidth',1.5);
    hold on
    plot(conf95_v3(9,:), 1*[maxY(2), maxY(2)],'b-o','LineWidth',1.5);

    set(gca,'YTickLabel',[])
    ylim([0.95 1.3])
    legend('Parametric Bootstrap - 95% range','Parametric Bootstrap - 2 Sigma',...
           'Residual Bootstrap - 95% range','Residual Bootstrap - 2 Sigma',...
           'Wild Bootstrap - 95% range','Wild Bootstrap - 2 Sigma',...
           'location','northoutside');
end

end



