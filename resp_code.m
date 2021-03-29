%% 6.1
% load the data
n = length(images);

landmark = load("landmark_pos_phys.mat");
landmark_pos = landmark.landmark_pos_phys;

target_nii = load_untouch_nii('0500.nii.gz');
cpg1_nii = load_untouch_nii('0100_cpp_region1.nii.gz');
cpg2_nii = load_untouch_nii('0100_cpp_region2.nii.gz');
source_nii = load_untouch_nii('0007.nii.gz');
dist_nii = load_untouch_nii('0007_sdt.nii.gz');

% % display an image 
figure;
dispNiiSlice(target_nii,"z",1)
figure;
dispNiiSlice(source_nii,"z",1)
figure;
dispNiiSlice(cpg1_nii,"z",1);
figure;
dispNiiSlice(cpg2_nii,"z",1)

% % deform an image with a B-spline transformation 
% [def_vol_nii, def_field_nii, dis_field_nii] = ...
%     deformNiiWithCPGsSliding(cpg1_nii, cpg2_nii, dist_nii, source_nii, target_nii);
% figure;
% dispNiiSlice(def_vol_nii,"z",1);

%%
% Just display all images one after another as the tasks asks for: 
% "display all the MR images to get a feel..."

% for k = 1:n
%   disp(k);
%   dispNiiSlice(images(k),"z",1);
%   drawnow;
% end

%%
% Side by side comparison for all deformations. Runs slow or might be my
% laptop. Set to only first 100 images.
% Plot titles dont work for some reason. Might be dispNiiSlice funtion.

% I also suggest to edit dispNiiSlice to have 'show_axes' as folse
% directly in the function.

% for k = 1:100
%   disp(k);
%   [def_vol_nii, ~, ~] = ...
%         deformNiiWithCPGsSliding(cpg1(k), cpg2(k), dist_nii, source_nii, images(k));
%     
%   t = tiledlayout(1,3);
%   nexttile
%   dispNiiSlice(source_nii,"z",1);
%   title('Source Image')
%   
%   nexttile
%   dispNiiSlice(def_vol_nii,"z",1);
%   title('Deformed Image')
%   
%   nexttile
%   dispNiiSlice(images(k),"z",1);
%   title('Target Image')
% 
%   drawnow;
% end

%% 6.2 Calculate the surrogate signal

x_20 = nan(n,1);
for k = 1:n
    
    target_img = images(k).img;
    target_img_y150 = target_img(:,150);

    for i = 1:(length(target_img_y150)-1)
        if (target_img_y150(i) <=20) && (target_img_y150(i+1) >20)
            x0 = i;
            y0 = target_img_y150(i);
            x1 = i+1;
            y1 = target_img_y150(i+1);
            break
        end
    end
    x_20(k) = (20-y0)*(x1-x0)/(y1-y0) +x0;
end

figure;
plot(x_20(1:end))
ylabel("skin position[index]")
xlabel("image number")

%% 6.3 
% using the matrix format in the Appendix is recommended 

% fit CP image: (44,38,1,1,2)
SI_deform_CP_44_38 = nan(n,1);
for k = 1:n   
    cpg1_img = cpg1(k).img;
    SI_deform_CP_44_38(k) = cpg1_img(44,38,1,1,2);
end

figure;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),':')
hold on;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),'bx')
ylim([30,55])

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];
Y_linear = X_linear*inv(X_linear'*X_linear)*X_linear'*Y_train;

figure;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),':')
hold on;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),'bx');
hold on;
plot(x_20(1:100),Y_linear,'r-');
ylim([30,55])
ylabel("surrogate value")
ylabel("control-point value")
hold off

%% Compute the coefficients for all CP
% Do it in matrix form:
% Y_train_mat 100xM contains the all the CP of the 100 first images, M is
% of length 67*67*2*2 = 17956
% It is divided into four sections horizontal concatenated:
% CP Region 1 AP - CP Region 1 SI - CP Region 2 AP - CP Region 2 SI

Y_train_mat = zeros(100,67*67*4);

for i=1:100
    Y_train_mat(i,:) = [reshape(cpg1(i).img(:,:,1,1,1), [1,67*67]), ...
                        reshape(cpg1(i).img(:,:,1,1,2), [1,67*67]), ...
                        reshape(cpg2(i).img(:,:,1,1,1), [1,67*67]), ...
                        reshape(cpg2(i).img(:,:,1,1,2), [1,67*67])];
end

% Linear model
S_lin = [x_20(1:100),ones(length(x_20(1:100)),1)];
C = S_lin\Y_train_mat;

% Polynomial of 2nd order
S_p2 = [x_20(1:100).^2,x_20(1:100),ones(length(x_20(1:100)),1)];
C_p2 = S_p2\Y_train_mat;

% Polynomial of 3rd order
S_p3 = [x_20(1:100).^3, x_20(1:100).^2,x_20(1:100),ones(length(x_20(1:100)),1)];
C_p3 = S_p3\Y_train_mat;

%% Get CP 44-38 and plot it

% retrieve the coefficient of CP Region 1 AP and reshape them in a 67x67 grid
C_reg1_AP = reshape(C(:,(4489*0+1):4489*1),[],67,67);
C_p2_reg1_AP = reshape(C_p2(:,(4489*0+1):4489*1),[],67,67);
C_p3_reg1_AP = reshape(C_p3(:,(4489*0+1):4489*1),[],67,67);

figure;
imshow(squeeze(C_reg1_AP(1,:,:))','DisplayRange',[])
colormap 
colorbar

% retrieve the coefficient of CP Region 1 SI and reshape them in a 67x67 grid
C_reg1_SI = reshape(C(:,(4489*1+1):4489*2),[],67,67);
C_p2_reg1_SI = reshape(C_p2(:,(4489*1+1):4489*2),[],67,67);
C_p3_reg1_SI = reshape(C_p3(:,(4489*1+1):4489*2),[],67,67);

figure;
imshow(squeeze(C_reg1_SI(1,:,:))','DisplayRange',[])
colormap 
colorbar

% retrieve the coefficient of CP Region 2 AP and reshape them in a 67x67 grid
C_reg2_AP = reshape(C(:,(4489*2+1):4489*3),[],67,67);
C_p2_reg2_AP = reshape(C_p2(:,(4489*2+1):4489*3),[],67,67);
C_p3_reg2_AP = reshape(C_p3(:,(4489*2+1):4489*3),[],67,67);

figure;
imshow(squeeze(C_reg1_AP(1,:,:))','DisplayRange',[])
colormap 
colorbar

% retrieve the coefficient of CP Region 2 SI and reshape them in a 67x67 grid
C_reg2_SI = reshape(C(:,(4489*3+1):4489*4),[],67,67);
C_p2_reg2_SI = reshape(C_p2(:,(4489*3+1):4489*4),[],67,67);
C_p3_reg2_SI = reshape(C_p3(:,(4489*3+1):4489*4),[],67,67);


cpg1_nii = load_untouch_nii('0100_cpp_region1.nii.gz');
cpg2_nii = load_untouch_nii('0100_cpp_region2.nii.gz');
figure;
subplot(1,2,1)
imshow(cpg1_nii.img(:,:,1,1,2)',[])
colormap 
colorbar
title("0100 cpp region1.nii.gz")
subplot(1,2,2)
imshow(cpg2_nii.img(:,:,1,1,2)',[])
colormap 
colorbar
title("0100 cpp region2.nii.gz")

% retrieve the Y_train of CP: Region 1 AP and reshape them in a 67x67 grid
Y_reg1_AP = reshape(Y_train_mat(:,(4489*0+1):4489*1),[],67,67);


% retrieve the Y_train of CP: Region 1 SI and reshape them in a 67x67 grid
Y_reg1_SI = reshape(Y_train_mat(:,(4489*1+1):4489*2),[],67,67);

% retrieve the Y_train of CP: Region 2 AP and reshape them in a 67x67 grid
Y_reg2_AP = reshape(Y_train_mat(:,(4489*2+1):4489*3),[],67,67);

% retrieve the Y_train of CP: Region 2 SI and reshape them in a 67x67 grid
Y_reg2_SI = reshape(Y_train_mat(:,(4489*3+1):4489*4),[],67,67);

%% Print the estimated models at CP 44-38
x=linspace(min(x_20(1:100))-0.5,max(x_20(1:100))+0.5,100)';
% get the corresponding coefficients
i = 44;
j = 38;

% region 1, AP
% get the corresponding coefficients
c1_r1_AP = C_reg1_AP(:,i,j);
c2_r1_AP = C_p2_reg1_AP(:,i,j);
c3_r1_AP = C_p3_reg1_AP(:,i,j);

figure;
subplot(2,2,1)
plot(x_20(1:100),Y_reg1_AP(:,44,38),'bx');
hold on;
plot(x,[x,ones(100,1)]*c1_r1_AP,'r-','LineWidth',1.5);
hold on
plot(x,[x.^2,x,ones(100,1)]*c2_r1_AP,'g-','LineWidth',1.5);
hold on
plot(x,[x.^3,x.^2,x,ones(100,1)]*c3_r1_AP,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 1, AP")
hold off

% region 2, AP
% get the corresponding coefficients
i = 44;
j = 38;
c1_r2_AP = C_reg2_AP(:,i,j);
c2_r2_AP = C_p2_reg2_AP(:,i,j);
c3_r2_AP = C_p3_reg2_AP(:,i,j);

subplot(2,2,2);
plot(x_20(1:100),Y_reg2_AP(:,44,38),'bx');
hold on;
plot(x,[x,ones(100,1)]*c1_r2_AP,'r-','LineWidth',1.5);
hold on
plot(x,[x.^2,x,ones(100,1)]*c2_r2_AP,'g-','LineWidth',1.5);
hold on
plot(x,[x.^3,x.^2,x,ones(100,1)]*c3_r2_AP,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 2, AP")
hold off

% region 1, SI
c1_r1_SI = C_reg1_SI(:,i,j);
c2_r1_SI = C_p2_reg1_SI(:,i,j);
c3_r1_SI = C_p3_reg1_SI(:,i,j);

subplot(2,2,3)
plot(x_20(1:100),Y_reg1_SI(:,44,38),'bx');
hold on;
plot(x,[x,ones(100,1)]*c1_r1_SI,'r-','LineWidth',1.5);
hold on
plot(x,[x.^2,x,ones(100,1)]*c2_r1_SI,'g-','LineWidth',1.5);
hold on
plot(x,[x.^3,x.^2,x,ones(100,1)]*c3_r1_SI,'m-','LineWidth',1.5);
ylim([30,55])
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 1, SI")
hold off

% region 2, SI
% get the corresponding coefficients
i = 44;
j = 38;
c1_r2_SI = C_reg2_SI(:,i,j);
c2_r2_SI = C_p2_reg2_SI(:,i,j);
c3_r2_SI = C_p3_reg2_SI(:,i,j);
subplot(2,2,4);
plot(x_20(1:100),Y_reg2_SI(:,44,38),'bx');
hold on;
plot(x,[x,ones(100,1)]*c1_r2_SI,'r-','LineWidth',1.5);
hold on
plot(x,[x.^2,x,ones(100,1)]*c2_r2_SI,'g-','LineWidth',1.5);
hold on
plot(x,[x.^3,x.^2,x,ones(100,1)]*c3_r2_SI,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 2, SI")
hold off
%% model fitting on test data
for i=1:1400
    Y_test_mat(i,:) = [reshape(cpg1(100+i).img(:,:,1,1,1), [1,67*67]), ...
                        reshape(cpg1(100+i).img(:,:,1,1,2), [1,67*67]), ...
                        reshape(cpg2(100+i).img(:,:,1,1,1), [1,67*67]), ...
                        reshape(cpg2(100+i).img(:,:,1,1,2), [1,67*67])];
end
n_test = 1400;
x_test =linspace(min(x_20(101:1500))-0.5,max(x_20(101:1500))+0.5,1400)';

% retrieve the Y_train of CP: Region 1 AP and reshape them in a 67x67 grid
Y_reg1_AP_test = reshape(Y_test_mat(:,(4489*0+1):4489*1),[],67,67);

% retrieve the Y_train of CP: Region 1 SI and reshape them in a 67x67 grid
Y_reg1_SI_test = reshape(Y_test_mat(:,(4489*1+1):4489*2),[],67,67);

% retrieve the Y_train of CP: Region 2 AP and reshape them in a 67x67 grid
Y_reg2_AP_test = reshape(Y_test_mat(:,(4489*2+1):4489*3),[],67,67);

% retrieve the Y_train of CP: Region 2 SI and reshape them in a 67x67 grid
Y_reg2_SI_test = reshape(Y_test_mat(:,(4489*3+1):4489*4),[],67,67);

% get the corresponding coefficients
i = 44;
j = 38;

% region 1, AP
% get the corresponding coefficients
c1_r1_AP = C_reg1_AP(:,i,j);
c2_r1_AP = C_p2_reg1_AP(:,i,j);
c3_r1_AP = C_p3_reg1_AP(:,i,j);

figure;
subplot(2,2,1)
plot(x_20(101:1500),Y_reg1_AP_test(:,44,38),'bx');
hold on;
plot(x_test,[x_test,ones(1400,1)]*c1_r1_AP,'r-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^2,x_test,ones(1400,1)]*c2_r1_AP,'g-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^3,x_test.^2,x_test,ones(1400,1)]*c3_r1_AP,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 1, AP")
hold off

% region 2, AP
% get the corresponding coefficients
i = 44;
j = 38;
c1_r2_AP = C_reg2_AP(:,i,j);
c2_r2_AP = C_p2_reg2_AP(:,i,j);
c3_r2_AP = C_p3_reg2_AP(:,i,j);

subplot(2,2,2);
plot(x_20(101:1500),Y_reg2_AP_test(:,44,38),'bx');
hold on;
plot(x_test,[x_test,ones(1400,1)]*c1_r2_AP,'r-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^2,x_test,ones(1400,1)]*c2_r2_AP,'g-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^3,x_test.^2,x_test,ones(1400,1)]*c3_r2_AP,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 2, AP")
hold off

% region 1, SI
c1_r1_SI = C_reg1_SI(:,i,j);
c2_r1_SI = C_p2_reg1_SI(:,i,j);
c3_r1_SI = C_p3_reg1_SI(:,i,j);

subplot(2,2,3)
plot(x_20(101:1500),Y_reg1_SI_test(:,44,38),'bx');
hold on;
plot(x_test,[x_test,ones(1400,1)]*c1_r1_SI,'r-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^2,x_test,ones(1400,1)]*c2_r1_SI,'g-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^3,x_test.^2,x_test,ones(1400,1)]*c3_r1_SI,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 1, SI")
hold off

% region 2, SI
% get the corresponding coefficients
i = 44;
j = 38;
c1_r2_SI = C_reg2_SI(:,i,j);
c2_r2_SI = C_p2_reg2_SI(:,i,j);
c3_r2_SI = C_p3_reg2_SI(:,i,j);
subplot(2,2,4);
plot(x_20(101:1500),Y_reg2_SI_test(:,44,38),'bx');
hold on;
plot(x_test,[x_test,ones(1400,1)]*c1_r2_SI,'r-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^2,x_test,ones(1400,1)]*c2_r2_SI,'g-','LineWidth',1.5);
hold on
plot(x_test,[x_test.^3,x_test.^2,x_test,ones(1400,1)]*c3_r2_SI,'m-','LineWidth',1.5);
legend;
xlabel("surrogate value")
ylabel("control-point value")
legend('data','linear',"2nd poly","3rd poly")
title("Region 2, SI")
hold off

%% Calculate the residual fitting error.
n_train =100;
res_lin = (Y_train_mat - S_lin*C);
res_p2 = (Y_train_mat - S_p2*C_p2);
res_p3 = (Y_train_mat - S_p2*C_p2);

mse_lin = sum((Y_train_mat - S_lin*C).^2,1)./n_train;
mse_p2 = sum((Y_train_mat - S_p2*C_p2).^2,1)./n_train;
mse_p3 = sum((Y_train_mat - S_p3*C_p3).^2,1)./n_train;


% retrieve the info CP: Region 1 AP and reshape them in a 67x67 grid
mse_reg1_AP1 = reshape(mse_lin(1,(4489*0+1):4489*1),67,67);

% retrieve info of CP: Region 1 SI and reshape them in a 67x67 grid
mse_reg1_SI1 = reshape(mse_lin(1,(4489*1+1):4489*2),67,67);

% retrieve info of CP: Region 2 AP and reshape them in a 67x67 grid
mse_reg2_AP1 =  reshape(mse_lin(1,(4489*2+1):4489*3),67,67);

% retrieve info of CP: Region 2 SI and reshape them in a 67x67 grid
mse_reg2_SI1 =  reshape(mse_lin(1,(4489*3+1):4489*4),67,67);
dist_map = dist_nii.img;

figure;
subplot(1,2,1)
imshow(mse_reg1_AP1',[]);
colormap hot
colorbar
title("MSE : region 1, AP")
subplot(1,2,2)
imshow(mse_reg2_AP1',[]);
colormap hot
colorbar
title("MSE : region 2, AP")

figure;
subplot(1,2,1)
imshow(mse_reg1_SI1',[]);
colormap hot
colorbar
title("MSE : region 1, SI")
subplot(1,2,2)
imshow(mse_reg2_SI1',[]);
colormap hot
colorbar
title("MSE : region 2, SI")

figure;
imshow(mse_reg1_SI1'+mse_reg2_SI1'+mse_reg1_AP1'+mse_reg2_AP1',[]);
colormap turbo
colorbar
title("MSE: combination of the four parts")

% retrieve the info CP: Region 1 AP and reshape them in a 67x67 grid
mse_reg1_AP2 = reshape(mse_p2(1,(4489*0+1):4489*1),67,67);

% retrieve info of CP: Region 1 SI and reshape them in a 67x67 grid
mse_reg1_SI2 = reshape(mse_p2(1,(4489*1+1):4489*2),67,67);

% retrieve info of CP: Region 2 AP and reshape them in a 67x67 grid
mse_reg2_AP2 =  reshape(mse_p2(1,(4489*2+1):4489*3),67,67);

% retrieve info of CP: Region 2 SI and reshape them in a 67x67 grid
mse_reg2_SI2 =  reshape(mse_p2(1,(4489*3+1):4489*4),67,67);


figure;
subplot(1,2,1)
imshow(mse_reg1_AP2',[]);
colormap hot
colorbar
title("MSE : region 1, AP")
subplot(1,2,2)
imshow(mse_reg2_AP2',[]);
colormap hot
colorbar
title("MSE : region 2, AP")

figure;
subplot(1,2,1)
imshow(mse_reg1_SI2',[]);
colormap hot
colorbar
title("MSE : region 1, SI")
subplot(1,2,2)
imshow(mse_reg2_SI2',[]);
colormap hot
colorbar
title("MSE : region 2, SI")

figure;
imshow(mse_reg1_SI2'+mse_reg2_SI2'+mse_reg1_AP2'+mse_reg1_AP2',[]);
colormap turbo
colorbar
title("MSE: 2nd poly")


% retrieve the info CP: Region 1 AP and reshape them in a 67x67 grid
mse_reg1_AP3 = reshape(mse_p3(1,(4489*0+1):4489*1),67,67);

% retrieve info of CP: Region 1 SI and reshape them in a 67x67 grid
mse_reg1_SI3 = reshape(mse_p3(1,(4489*1+1):4489*2),67,67);

% retrieve info of CP: Region 2 AP and reshape them in a 67x67 grid
mse_reg2_AP3 =  reshape(mse_p3(1,(4489*2+1):4489*3),67,67);

% retrieve info of CP: Region 2 SI and reshape them in a 67x67 grid
mse_reg2_SI3 =  reshape(mse_p3(1,(4489*3+1):4489*4),67,67);


figure;
subplot(1,2,1)
imshow(mse_reg1_AP3',[]);
colormap hot
colorbar
title("MSE : region 1, AP")
subplot(1,2,2)
imshow(mse_reg2_AP3',[]);
colormap hot
colorbar
title("MSE : region 2, AP")

figure;
subplot(1,2,1)
imshow(mse_reg1_SI3',[]);
colormap hot
colorbar
title("MSE : region 1, SI")
subplot(1,2,2)
imshow(mse_reg2_SI3',[]);
colormap hot
colorbar
title("MSE : region 2, SI")

figure;
subplot(1,3,1)
imshow(mse_reg1_SI1'+mse_reg2_SI1'+mse_reg1_AP1'+mse_reg2_AP1',[]);
colormap turbo
colorbar
title("MSE: linear reg")
subplot(1,3,2)
imshow(mse_reg1_SI2'+mse_reg2_SI2'+mse_reg1_AP2'+mse_reg2_AP2',[]);
colormap turbo
colorbar
title("MSE: 2nd poly")
subplot(1,3,3)
imshow(mse_reg1_SI3'+mse_reg2_SI3'+mse_reg1_AP3'+mse_reg2_AP3',[]);
colormap turbo
colorbar
title("MSE: 3rd poly")
%% bootstrap 
% version one: parametric bootstap
% resampling Y_train

n_boot = 1000;
C1_boot_v1 = nan(n_boot,2);
C2_boot_v1 = nan(n_boot,3);
C3_boot_v1 = nan(n_boot,4);

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];
Y_linear = X_linear*inv(X_linear'*X_linear)*X_linear'*Y_train;
sigma = sqrt(sum((Y_linear-Y_train).^2)/(100-2));

coef_v1 = nan(n_boot,2);
for i = 1: n_boot
    T = 100; 
    Y_train_boot_v1 = Y_train + normrnd(0,sigma,[100,1]);
    coef_v1(i,:) = inv(X_linear'*X_linear)*X_linear'*Y_train_boot_v1;
end

figure;
subplot(1,2,1)
hist(coef_v1(:,1))
subplot(1,2,2)
hist(coef_v1(:,2))

% version two : residual bootstrap
% version three: wild bootstap

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];


%% % version one: parametric bootstap
% resampling Y_train

n_boot = 1000;
C1_boot_v1 = nan(n_boot,2);
C2_boot_v1 = nan(n_boot,3);
C3_boot_v1 = nan(n_boot,4);

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];
Y_linear = X_linear*inv(X_linear'*X_linear)*X_linear'*Y_train;
sigma = sqrt(sum((Y_linear-Y_train).^2)/(100-2));

coef_v1 = nan(n_boot,2);
for i = 1: n_boot
    T = 100; 
    Y_train_boot_v1 = Y_train + normrnd(0,sigma,[100,1]);
    coef_v1(i,:) = inv(X_linear'*X_linear)*X_linear'*Y_train_boot_v1;
end

figure;
subplot(1,2,1)
hist(coef_v1(:,1))
subplot(1,2,2)
hist(coef_v1(:,2))

% version two : residual bootstrap
% version three: wild bootstap

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];

% S0: mean 
% 2.5% - 97.5%
coef_v1_1 = coef_v1(:,1);
sorted_S0 = sort(coef_v1_1);

S0_025_boot3 = sorted_S0(0.025*n_boot); %4.1449e+03
S0_975_boot3 = sorted_S0(0.975*n_boot); %4.3646e+03

% 2 sigma 
S0_2sigma_left_boot3 = mean(sorted_S0)+2*std(sorted_S0); %4.3658e+03
S0_2sigma_right_boot3 = mean(sorted_S0)-2*std(sorted_S0); % 4.1453e+03

figure;
hist(coef_v1_1);
hold on;
plot([S0_025_boot3,S0_975_boot3],[270,270],"r","linewidth",2);
hold on;
plot([S0_2sigma_left_boot3,S0_2sigma_right_boot3],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
%title("Wild Bootstrap: S0")
hold off;


% version two: residual bootstrap
n_boot = 1000;
C1_boot_v2 = nan(n_boot,2);
C2_boot_v2 = nan(n_boot,3);
C3_boot_v2 = nan(n_boot,4);

Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];
Y_linear = X_linear*inv(X_linear'*X_linear)*X_linear'*Y_train;
sigma = sqrt(sum((Y_linear-Y_train).^2)/(100-2));

coef_v2 = nan(n_boot,2);
for i = 1: n_boot
    T = 100; 
    Y_train_boot_v2 = Y_train + normrnd(0,sigma,[100,1]);
    coef_v2(i,:) = inv(X_linear'*X_linear)*X_linear'*Y_train_boot_v2;
end

figure;
subplot(1,2,1)
hist(coef_v2(:,1))
subplot(1,2,2)
hist(coef_v2(:,2))

% version three : wild bootstrap
coef_v3 = nan(n_boot,2);
Y_train = SI_deform_CP_44_38(1:100);
X_linear = [ones(length(x_20(1:100)),1),x_20(1:100)];
Y_linear = X_linear*inv(X_linear'*X_linear)*X_linear'*Y_train;
noises = Y_train - Y_linear;

for i = 1: n_boot
    Y_train_boot_v3 = Y_linear + noises(ceil(rand(n_train,1)*n_train)).* normrnd(0,1,[100,1]);
    coef_v3(i,:) = inv(X_linear'*X_linear)*X_linear'*Y_train_boot_v3;
end

figure;
subplot(1,2,1)
hist(coef_v3(:,1))
subplot(1,2,2)
hist(coef_v3(:,2))