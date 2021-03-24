%% Code for 2nd order Polynomial 


%% Model Images poly 2

% Deform the source image with the Model transformations

Z = 1;
for i = 1:1400
    
        [modImg2(Z), modDef2(Z), modDis2(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p2(:,i), cpg2_p2(:,i), SDMS, source_image, images(:,i+100));

        Z = Z+1;
end
%% Masking

% Deform the source mask image with the model transformations

% Variable name = Deformed Source Map Images (mDSMI)

Z = 1;
for i = 1:1400
    
        [mDSMI2(Z), mDSMI_def2(Z), mDSMI_dis2(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_p2(:,i), cpg2_p2(:,i), SDMS, source_mask, images(:,i+100));

        Z = Z+1;
end

%% Map the mask onto the deformed modelled images


for i = 1:1400
    for j = 1:160
        for k = 1:160
            if mDSMI2(:,i).img(j,k) == 0
                modImg2(:,i).img(j,k) = 0;
                modDef2(:,i).img(j,k) = 0;                
            end      
        end 
    end
end

%% Plot arbitrary figures to see if its working


figure(2);
subplot(1,3,1)
dispNiiSlice(modImg2(:,1200),"z",1)
title('Model Deformed Masked MRI Image')

subplot(1,3,2)
dispNiiSlice(mDSMI2(:,1200),"z",1)
title('Model Deformed Mask Image')

subplot(1,3,3)
dispNiiSlice(modDef2(:,1200),"z",1)
title('Model Deformed Mask Deformation Field Image')



%% Deformation Field Error Again

% Create a new struct with the same values as rDSMI for the for loop
DEF_def22 = rDSMI_def;

% Calculate the deformation field error for every voxel/pixel among the...
% ...1400 target imgs.

for i = 1:1400
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
         DEF_def22(i).img(j,k) = regIm_def(:,i).img(j,k) - modDef2(:,i).img(j,k);  
         
        end 
    end
end

%% Making the deformation with just the registration mask
DEF_def2_mask2 = DEF_def22;

for i = 1:25
    for j = 1:160
        for k = 1:160
            if rDSMI_def(:,i).img(j,k) == 0
                DEF_def2_mask2(i).img(j,k) = 0;
            end      
        end 
    end
end
%%
figure(3);
subplot(1,2,1)
dispNiiSlice(DEF_def2_mask(:,25),"z",1,[-2 2])
colorbar;
subplot(1,2,2)
dispNiiSlice(DEF_def2_mask2(:,25),"z",1,[-2 2])
colorbar;
%%
L2_norm2 = DEF_def2;
mean_Val = [];


for i = 1:10
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
            L2_norm2(i).img(j,k) = sqrt((regIm_def(:,i).img(j,k) - modDef2(:,i).img(j,k))^2);
            
            if rDSMI_def(:,i).img(j,k) == 0
                L2_norm2(i).img(j,k) = 0;
            end      
        end 
    end
end

%%
figure(4);
subplot(1,2,1)
dispNiiSlice(L2_norm(:,9),"z",1,[-2 2])
colorbar;
subplot(1,2,2)
dispNiiSlice(L2_norm2(:,9),"z",1,[-2 2])
colorbar;