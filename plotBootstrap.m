function plotBootstrap(twoSigma,conf95)

%linear
figure;
subplot(1,2,1)
title('Param1');
maxY = ylim;
hold on
plot(conf95(1,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(1,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(1,2,2)
title('Param2');
maxY = ylim;
hold on
plot(conf95(2,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(2,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');


%2nd order
figure;
subplot(2,2,1)
title('Param1');
maxY = ylim;
hold on
plot(conf95(3,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(3,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(2,2,2)
title('Param2');
maxY = ylim;
hold on
plot(conf95(4,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(4,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(2,2,3)
title('Param3');
maxY = ylim;
hold on
plot(conf95(5,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(5,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');


%3rd order
figure;
subplot(2,2,1)
title('Param1');
maxY = ylim;
hold on
plot(conf95(6,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(6,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(2,2,2)
title('Param2');
maxY = ylim;
hold on
plot(conf95(7,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(7,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(2,2,3)
title('Param3');
maxY = ylim;
hold on
plot(conf95(8,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(8,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');

subplot(2,2,4)
title('Param4');
maxY = ylim;
hold on
plot(conf95(9,:), 1.05*[maxY(2), maxY(2)],'--o');
hold on
plot(twoSigma(9,:), 1.1*[maxY(2), maxY(2)], '--o');
legend('coef','95% range','2 Sigma','location','northoutside');
end

