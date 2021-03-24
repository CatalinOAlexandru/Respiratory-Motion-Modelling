%% 6.1
% load the data
n = length(images);

landmark = load("landmark_pos_phys.mat");
landmark_pos = landmark.landmark_pos_phys;

target_nii = load_untouch_nii('0500.nii.gz');
cpg1_nii = load_untouch_nii('0500_cpp_region1.nii.gz');
cpg2_nii = load_untouch_nii('0500_cpp_region2.nii.gz');
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

% deform an image with a B-spline transformation 
[def_vol_nii, def_field_nii, dis_field_nii] = ...
    deformNiiWithCPGsSliding(cpg1_nii, cpg2_nii, dist_nii, source_nii, target_nii);
figure;
dispNiiSlice(def_vol_nii,"z",1);

%%
% Just display all images one after another as the tasks asks for: 
% "display all the MR images to get a feel..."

for k = 1:n
  disp(k);
  dispNiiSlice(images(k),"z",1);
  drawnow;
end

%%
% Side by side comparison for all deformations. Runs slow or might be my
% laptop. Set to only first 100 images.
% Plot titles dont work for some reason. Might be dispNiiSlice funtion.

% I also suggest to edit dispNiiSlice to have 'show_axes' as folse
% directly in the function.

for k = 1:100
  disp(k);
  [def_vol_nii, ~, ~] = ...
        deformNiiWithCPGsSliding(cpg1(k), cpg2(k), dist_nii, source_nii, images(k));
    
  t = tiledlayout(1,3);
  nexttile
  dispNiiSlice(source_nii,"z",1);
  title('Source Image')
  
  nexttile
  dispNiiSlice(def_vol_nii,"z",1);
  title('Deformed Image')
  
  nexttile
  dispNiiSlice(images(k),"z",1);
  title('Target Image')

  drawnow;
end

%% 6.2 Calculate the surrogate signal

x_20 = nan(n,1);
for k = 1:n
%     if k<=9
%         target_nii = load_untouch_nii(['000',int2str(k),'.nii.gz']);
%     elseif (10 <=k) && (k<=99)
%         target_nii = load_untouch_nii(['00',int2str(k),'.nii.gz']);
%     elseif (100 <=k) && (k <= 999)
%         target_nii = load_untouch_nii(['0',int2str(k),'.nii.gz']);
%     else
%         target_nii = load_untouch_nii([int2str(k),'.nii.gz']);
%     end
    
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
xlabel('frame');
ylabel('pixel');

%% 6.3 
% work already done for 6.3 is: fit CP image: (44,38,1,1,2)
% work has not been done yet for 6.3 is: fit a separate model to every CP deformation
% using the matrix format in the Appendix is recommended 

% fit CP image: (44,38,1,1,2)
SI_deform_CP_44_38 = nan(n,1);
for k = 1:n
%     if k<=9
%         cpg1_nii = load_untouch_nii(['000',int2str(k),'_cpp_region1.nii.gz']);
%     elseif (10 <=k) && (k<=99)
%         cpg1_nii = load_untouch_nii(['00',int2str(k),'_cpp_region1.nii.gz']);
%     elseif (100 <=k) && (k <= 999)
%         cpg1_nii = load_untouch_nii(['0',int2str(k),'_cpp_region1.nii.gz']);
%     else
%         cpg1_nii = load_untouch_nii([int2str(k),'_cpp_region1.nii.gz']);
%     end
    
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
xlabel("surrogate value")
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
% retrieve the coefficient of CP Region 1 SI and reshape them in a 67x67 grid
C_reg1_SI = reshape(C(:,(4489+1):8978),[],67,67);
C_p2_reg1_SI = reshape(C_p2(:,(4489+1):8978),[],67,67);
C_p3_reg1_SI = reshape(C_p3(:,(4489+1):8978),[],67,67);

% get the corresponding coefficients
C_44_38 = C_reg1_SI(:,44,38);
C_p2_44_38 = C_p2_reg1_SI(:,44,38);
C_p3_44_38 = C_p3_reg1_SI(:,44,38);


%% Print the estimated models at CP 44-38
figure;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),':')
hold on;
plot(x_20(1:100),SI_deform_CP_44_38(1:100),'bx');
hold on;
plot(x_20(1:100),S_lin*C_44_38,'r-');
hold on
[sorted_p2,ind] = sort(S_p2*C_p2_44_38,'descend');
plot(x_20(ind), sorted_p2,'g-');
hold on
[sorted_p3,ind] = sort(S_p3*C_p3_44_38,'descend');
plot(x_20(ind), sorted_p3,'m-');
ylim([30,55])
legend;
xlabel("surrogate value")
ylabel("control-point value")
hold off

%% Calculate the residual fitting error.

res_lin = (Y_train_mat - S_lin*C);
res_p2 = (Y_train_mat - S_p2*C_p2);
res_p3 = (Y_train_mat - S_p2*C_p2);

mse_lin = sum((Y_train_mat - S_lin*C).^2,1)./(n_train- 2);
mse_p2 = sum((Y_train_mat - S_p2*C_p2).^2,1)./(n_train- 3);
mse_p3 = sum((Y_train_mat - S_p3*C_p3).^2,1)./(n_train- 4);

% retrieve the info CP: Region 1 AP and reshape them in a 67x67 grid
mse_reg1_AP = reshape(mse_lin(1,(4489*0+1):4489*1),67,67);

% retrieve info of CP: Region 1 SI and reshape them in a 67x67 grid
mse_reg1_SI = reshape(mse_lin(1,(4489*1+1):4489*2),67,67);

% retrieve info of CP: Region 2 AP and reshape them in a 67x67 grid
mse_reg2_AP =  reshape(mse_lin(1,(4489*2+1):4489*3),67,67);

% retrieve info of CP: Region 2 SI and reshape them in a 67x67 grid
mse_reg2_SI =  reshape(mse_lin(1,(4489*3+1):4489*4),67,67);

figure;
imshow(mse_reg1_AP');
figure;
imshow(mse_reg2_AP');


figure;
imshow(mse_reg1_SI');
figure;
imshow(mse_reg2_SI');
