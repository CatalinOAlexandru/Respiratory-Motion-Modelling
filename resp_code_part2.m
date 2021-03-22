%% 6.5 Create the Masks and Transformations and Deformation Field Errors

% Pre-requisites

% Signed Distance Map Segmentation = SDMS
SDMS = load_untouch_nii('0007_sdt.nii.gz');

% Source mask
source_mask = load_untouch_nii('0007_mask.nii.gz');

% Source image
source_image = load_untouch_nii('0007.nii.gz');

%% Registration Images

% Deform the source image with the registration transformations
% Registation images = regIm
Z = 1;
for i = 101:1500
    
        [regIm(Z), regIm_def(Z), regIm_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_image, images(:,i));

        Z = Z+1;
end

%% Masking The Registration Images

% Deform the source mask image with the registration transformations
% The deformNiiWithCPGsSliding will also calculate the Deformation Field
% Signed Distance Map Segmentation = SDMS
% Source_image = SI
% Target_image = TI
% [The_Image , deformation_Field, distortion_Field] = ...
% ...deformNiiWithCPGsSliding(registration_Region1, registration_Region2, SDMS, SI, TI);

% Creating the registration Deformed Source Map Images (rDSMI)

Z = 1;
for i = 101:1500
    
        [rDSMI(Z), rDSMI_def(Z), rDSMI_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1(:,i), cpg2(:,i), SDMS, source_mask, images(:,i));

        Z = Z+1;
end



% Map the mask onto the registration images

% Create New registration Masked images (NRMI)

NRMI = regIm;

for i = 1:1400
    for j = 1:160
        for k = 1:160
            if rDSMI(:,i).img(j,k) == 0
                NRMI(:,i).img(j,k) = 0;
                rDSMI_def(:,i).img(j,k) = 0;                
            end      
        end 
    end
end


%% Plot arbitrary figures to see if its working


figure(1);
subplot(1,3,1)
dispNiiSlice(NRMI(:,108),"z",1)
title('Registration Deformed Masked MRI Image')

subplot(1,3,2)
dispNiiSlice(regIm(:,108),"z",1)
title('Registration MRI Image')

subplot(1,3,3)
dispNiiSlice(rDSMI_def(:,108),"z",1)
title('Registration Deformed Mask Deformation Field Image')

%% Copy the header from all of region 1 and region 2

Z = 1;
for i = 1:1400
    hdr1(Z) = cpg1(:,i).hdr;
    hdr2(Z) = cpg2(:,i).hdr;
    Z = Z+1;
end 

%% Retrieve the signal for the 1400 images and input into linear model

% Re-label the variables and matrices to understand whats going on

% linear
linearsurrogateSignals = [x_20(101:1500),ones(length(x_20(101:1500)),1)];
linearcoefficients = C;
lineartransformations = linearsurrogateSignals*linearcoefficients;

% retrieve cpg's
[cpg1_lin, cpg2_lin] = transForm(lineartransformations,hdr1,hdr2);

%% Polynomial of 2nd order

p2surrogateSignals = [x_20(101:1500).^2,x_20(101:1500),ones(length(x_20(101:1500)),1)];
p2coefficients = C_p2;
p2transformations = p2surrogateSignals*p2coefficients;

% retrieve cpg's
[cpg1_p2, cpg2_p2] = transForm(p2transformations,hdr1,hdr2);


%% Polynomial of 3rd order

p3surrogateSignals = [x_20(101:1500).^3, x_20(101:1500).^2,x_20(101:1500),ones(length(x_20(101:1500)),1)];
p3coefficient = C_p3;
p3transformations = p3surrogateSignals*p3coefficient;

% retrieve cpg's
[cpg1_p3, cpg2_p3] = transForm(p3transformations,hdr1,hdr2);
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
                   % Linear Deformation Fielf Error %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Model Images

% Deform the source image with the Model transformations

Z = 1;
for i = 1:1400
    
        [modImg(Z), modDef(Z), modDis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_lin(:,i), cpg2_lin(:,i), SDMS, source_image, images(:,i+100));

        Z = Z+1;
end
%% Masking

% Deform the source mask image with the model transformations

% Variable name = Deformed Source Map Images (mDSMI)

Z = 1;
for i = 1:1400
    
        [mDSMI(Z), mDSMI_def(Z), mDSMI_dis(Z)] = ...
            deformNiiWithCPGsSliding(cpg1_lin(:,i), cpg2_lin(:,i), SDMS, source_mask, images(:,i+100));

        Z = Z+1;
end

%% Map the mask onto the deformed modelled images


for i = 1:1400
    for j = 1:160
        for k = 1:160
            if mDSMI(:,i).img(j,k) == 0
                modImg(:,i).img(j,k) = 0;
                modDef(:,i).img(j,k) = 0;                
            end      
        end 
    end
end

%% Plot arbitrary figures to see if its working


figure(2);
subplot(1,3,1)
dispNiiSlice(modImg(:,1200),"z",1)
title('Model Deformed Masked MRI Image')

subplot(1,3,2)
dispNiiSlice(mDSMI(:,1200),"z",1)
title('Model Deformed Mask Image')

subplot(1,3,3)
dispNiiSlice(modDef(:,1200),"z",1)
title('Model Deformed Mask Deformation Field Image')

%% Deformation Field Error

% Create a new struct with the same values as rDSMI for the for loop
DEF_def = rDSMI_def;

% Calculate the deformation field error for every voxel/pixel among the...
% ...1400 target imgs.

for i = 1:1400
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
         DEF_def(i).img(j,k) = rDSMI_def(:,i).img(j,k) - modDef(:,i).img(j,k);  
         
        end 
    end
end
    
%% Plot arbitrary figures to see if its working

figure(3);
subplot(1,3,1)
dispNiiSlice(L2_norm(:,1),"z",1)

subplot(1,3,2)
dispNiiSlice(DEF_def(:,500),"z",1)

subplot(1,3,3)
dispNiiSlice(DEF_def(:,1000),"z",1)

%% Calculating the L2 norm of each image

L2_norm = DEF_def;
mean_Val = [];

% histo=zeros(10,160,160);

for i = 1:1400
    for j = 1:160
        for k = 1:160
            % Registration Deformation field minus Model Deformation field
            L2_norm(i).img(j,k) = sqrt((rDSMI_def(:,i).img(j,k) - modDef(:,i).img(j,k))^2);
            
            % tried to remove all the NaN values and failed miserably
            if isnan(L2_norm(i).img(j,k))
                L2_norm(i).img(j,k)=0;
            end
            
         
        end 
    end
end
%% Other statistical measures

% measure the average error value per image

mean_Val = [];

for i =1:1400
    
    SimpleArray = L2_norm(:,i).img;
    NewSimpleArray = SimpleArray(SimpleArray ~= 0);
    mean_Val(i) = nanmean(NewSimpleArray);
    
end

% Measure the error in the Superior Inferior (SI) position:


% Measure the error in the Anterior Posterior (AP) position

%% Plot Average value against corresponding signal

figure(3);
subplot(1,2,1)
plot(mean_Val);

subplot(1,2,2)
plot(x_20(100:1500,1));

