% 6.1
% load the data

landmark = load("landmark_pos_phys.mat")
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
    deformNiiWithCPGsSliding(cpg1_nii, cpg2_nii, dist_nii, source_nii, target_nii)
figure;
dispNiiSlice(def_vol_nii,"z",1)

% 6.2
x_20 = nan(1500,1);
for k = 1:1500
    if k<=9
        target_nii = load_untouch_nii(['000',int2str(k),'.nii.gz']);
    elseif (10 <=k) && (k<=99)
        target_nii = load_untouch_nii(['00',int2str(k),'.nii.gz']);
    elseif (100 <=k) && (k <= 999)
        target_nii = load_untouch_nii(['0',int2str(k),'.nii.gz']);
    else
        target_nii = load_untouch_nii([int2str(k),'.nii.gz']);
    end
    
    target_img = target_nii.img;
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
plot(x_20(1:300))

% 6.3 
% work already done for 6.3 is: fit CP image: (44,38,1,1,2)
% work has been done yet for 6.3 is: fit a separate model to every CP deformation
% using the matrix format in the Appendix is recommended 

% fit CP image: (44,38,1,1,2)
SI_deform_CP_44_38 = nan(1500,1);
for k = 1:1500
    if k<=9
        cpg1_nii = load_untouch_nii(['000',int2str(k),'_cpp_region1.nii.gz']);
    elseif (10 <=k) && (k<=99)
        cpg1_nii = load_untouch_nii(['00',int2str(k),'_cpp_region1.nii.gz']);
    elseif (100 <=k) && (k <= 999)
        cpg1_nii = load_untouch_nii(['0',int2str(k),'_cpp_region1.nii.gz']);
    else
        cpg1_nii = load_untouch_nii([int2str(k),'_cpp_region1.nii.gz']);
    end
    
    cpg1_img = cpg1_nii.img;
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

