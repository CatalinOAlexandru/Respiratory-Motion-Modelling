%  ========================== other bootstrap methods (non-parametric)===
% 
%  ========================== classical bootstrap (bad performance, not used) ==========================

% step 1 obtain the residuals

% the parameters S0, d and f in the ball and stick model. try other voxels
Avox = dwis(:,92,65,72); % 108*1

% step 2: generate bootstrap data using sampling with return
K = size(Avox,1);
bootstrap_iters = 1000;
bootstrap_sample = nan(K,bootstrap_iters);

n_iters = 30;
for i = 1:bootstrap_iters
    for j =1:K 
        bootstrap_sample(j,i) = Avox(ceil(rand*K));
    end
end

% step 3: obtain parameter estimations from bootstrap samples
bootstrap_results_voxel1 = nan(bootstrap_iters,6);
h=optimset('Display','off',...
           'MaxFunEvals',20000,...
           'Algorithm','quasi-newton',... 
           'TolX',1e-10,... 
           'TolFun',1e-10);
tic
for k = 1:bootstrap_iters
        bs_result = nan(n_iters,6);
        for i = 1:n_iters
            try
            Avox = bootstrap_sample(:,k);
            if (min(Avox) > 0)
                % add normally distributed random numbers to the one we used above; 
                % the scale of each random number should reflect the scale of each individual parameter to ensure large enough perturbations.
    
                %S0 = x(1)^2; diff = x(2)^2;f = exp(x(3))/(1+exp(x(3))); theta = x(4); phi = x(5);
 
                x1_start = normrnd(40,10);  % S0 = x(1)^2 in [0,4000], 
                x2_start = normrnd(0,0.5);  % diff = x(2)^2 in [0,1], so x2 in N(0,1) 
                x3_start = normrnd(0,5); % f = exp(x(3))/(1+exp(x(3))) in [0,1] so x3 in [0,6]
                x4_start = normrnd(0,5); % theta = 2*pi*exp(x(4))/(1+exp(x(4))) in [0,2*pi] so x4 in [0,6]
                x5_start = normrnd(0,5); % phi = 2*pi*exp(x(5))/(1+exp(x(5))) in [0,2*pi] so x5 in [0,6]
    
                % Inverse transformation to maintain the same starting point
                startx = [x1_start x2_start x3_start x4_start x5_start];
                % Define various options for the non-linear fitting % algorithm.

                % Now run the fitting
                [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transform',startx,h,Avox,bvals,qhat);
            
                % parameter_hat: fitted parameter values S0, d, f, theta, phi
                x = parameter_hat;
                S0 = x(1)^2;
                diff = x(2)^2;
                f = exp(x(3))/(1+exp(x(3)));
                % f + f*exp(x(3)) = exp(x(3)) => exp(x(3)) = f/(1-f) => x(3) = log(f/(1-f))
                theta = x(4);%2*pi*exp(x(4))/(1+exp(x(4)));
                phi = x(5);%2*pi*exp(x(5))/(1+exp(x(5)));
    
                bs_result(i,:) = [S0,diff,f,theta,phi,RESNORM];
            end
            end
         end
        round_min = round(min(bs_result(:,6)),3);
        resnorm_all = bs_result(:,6); 
        %percent_minimum = length(find(round(resnorm_all,3) == round_min ))/n_iters
        min_outcomes = bs_result(find(round(bs_result(:,6),3)==round_min),:);
        opt_index = find(min_outcomes(:,2)==min(min_outcomes(:,2)));
        opt_outcomes = min_outcomes(opt_index,:);
        if size(opt_outcomes,1) ~= 0
            S0 = opt_outcomes(1);diff= opt_outcomes(2);f= opt_outcomes(3);theta= opt_outcomes(4); phi= opt_outcomes(5);RESNORM = opt_outcomes(6);
            bootstrap_results_voxel1(k,:) = [S0,diff,f,theta,phi,RESNORM];
        else
            bootstrap_results_voxel1(k,:) = [nan,nan,nan,nan,nan,nan];
        end
end
toc  
% 304.490796 seconds. n_iters = 30; bootstrap_iters = 1000;
% order of parameters: [S0,diff,f,theta,phi,RESNORM]

% S0:
% 2.5% - 97.5%
sorted_S0 = sort(bootstrap_results_voxel1(:,1));
S0_025 = sorted_S0(0.025*bootstrap_iters);
S0_975 = sorted_S0(0.975*bootstrap_iters);

% 2 sigma 
S0_2sigma_left = mean(sorted_S0)+2*std(sorted_S0);
S0_2sigma_right = mean(sorted_S0)-2*std(sorted_S0);

figure;
hist(bootstrap_results_voxel1(:,1));
hold on;
plot([S0_025,S0_975],[270,270],"r","linewidth",2);
hold on;
plot([S0_2sigma_left,S0_2sigma_right],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Bootstrap: S0")
hold off;

% d 
% 2.5% - 97.5%
sorted_d = sort(bootstrap_results_voxel1(:,2));
d_025 = sorted_d(0.025*bootstrap_iters);
d_975 = sorted_d(0.975*bootstrap_iters);

% 2 sigma 
d_2sigma_left = mean(sorted_d)+2*std(sorted_d);
d_2sigma_right = mean(sorted_d)-2*std(sorted_d);

figure;
hist(bootstrap_results_voxel1(:,2));
hold on;
plot([d_025,d_975],[270,270],"r","linewidth",2);
hold on;
plot([d_2sigma_left,d_2sigma_right],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Bootstrap: d")
hold off;

% f
figure;
% 2.5% - 97.5%
sorted_f = sort(bootstrap_results_voxel1(:,3));
f_025 = sorted_f(0.025*bootstrap_iters);
f_975 = sorted_f(0.975*bootstrap_iters);

% 2 sigma 
f_2sigma_left = mean(sorted_f)+2*std(sorted_f);
f_2sigma_right = mean(sorted_f)-2*std(sorted_f);

figure;
hist(bootstrap_results_voxel1(:,3));
hold on;
plot([f_025,f_975],[280,280],"r","linewidth",2);
hold on;
plot([f_2sigma_left,f_2sigma_right],[290,290],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Bootstrap: f")
hold off;

% ========================== residual bootstrap ==========================
% step 1 obtain the residuals

% the parameters S0, d and f in the ball and stick model. try other voxels
Avox = dwis(:,92,65,72); % 108*1

% step 1: fit the model (several iterations for global minimiser) 
n_iters = 20;

bs_result = nan(n_iters,(6));
Avox = dwis(:,92,65,72);
h=optimset('Display','off',...
           'MaxFunEvals',20000,...
           'Algorithm','quasi-newton',... 
           'TolX',1e-10,... 
           'TolFun',1e-10);

for i = 1:n_iters
    try
    % add normally distributed random numbers to the one we used above; 
    % the scale of each random number should reflect the scale of each individual parameter to ensure large enough perturbations.
    
    %S0 = x(1)^2; diff = x(2)^2;f = exp(x(3))/(1+exp(x(3))); theta = x(4); phi = x(5);
 
    x1_start = normrnd(40,10); % S0 = x(1)^2 in [0,4000], 
    x2_start = normrnd(0,0.5); % diff = x(2)^2 in [0,1], so x2 in N(0,1) 
    x3_start = normrnd(0,5); % f = exp(x(3))/(1+exp(x(3))) in [0,1] so x3 in [0,6]
    x4_start = normrnd(0,5); % theta = 2*pi*exp(x(4))/(1+exp(x(4))) in [0,2*pi] so x4 in [0,6]
    x5_start = normrnd(0,5); % phi = 2*pi*exp(x(5))/(1+exp(x(5))) in [0,2*pi] so x5 in [0,6]
    
    % Inverse transformation to maintain the same starting point
    startx = [x1_start x2_start x3_start x4_start x5_start];
    % Define various options for the non-linear fitting % algorithm.

    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transform',startx,h,Avox,bvals,qhat);

    % parameter_hat: fitted parameter values S0, d, f, theta, phi
    x = parameter_hat;
    S0 = x(1)^2;
    diff = x(2)^2;
    f = exp(x(3))/(1+exp(x(3)));
    % f + f*exp(x(3)) = exp(x(3)) => exp(x(3)) = f/(1-f) => x(3) = log(f/(1-f))
    theta =x(4); %2*pi*exp(x(4))/(1+exp(x(4)));
    phi = x(5);%2*pi*exp(x(5))/(1+exp(x(5)));

    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
    fitted = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(- bvals*diff));
    
    bs_result(i,:) = [S0,diff,f,theta,phi,RESNORM];
    end
end


resnorm_all = bs_result(:,6);
min_res = min(resnorm_all);
percent_minimum = length(find(abs(resnorm_all-min_res)<0.0001))/n_iters;
%percent_voxel(:,m,n) = [m,n,percent_minimum];

opt_outcomes = bs_result(find(abs(resnorm_all-min_res)<0.0001),:);
opt_outcomes = opt_outcomes(1,:)

S0 = opt_outcomes(1);diff= opt_outcomes(2);f= opt_outcomes(3);theta= opt_outcomes(4); phi= opt_outcomes(5);

fibdir_opt = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad_opt = sum(qhat.*repmat(fibdir_opt, [length(qhat) 1])');
fitted_opt = S0*(f*exp(-bvals*diff.*(fibdotgrad_opt.^2)) + (1-f)*exp(- bvals*diff));

noises = fitted_opt - Avox'; % 1*108

% step 2: use noises obtained before to form new samples and estimate x_hat
K = size(fitted_opt,2);
bootstrap_iters = 1000;
bootstrap_sample = nan(K,bootstrap_iters);

n_iters = 30;
for i = 1:bootstrap_iters
    for j =1:K 
        bootstrap_sample(j,i) = fitted_opt(j)+ noises(ceil(rand*K));
        if isnan(bootstrap_sample(j,i))
            bootstrap_sample(j,i) = fitted_opt(j)+ noises(ceil(rand*K));
        end
    end
end

% step 3: obtain parameter estimations from bootstrap samples
bootstrap_results_voxel1 = nan(bootstrap_iters,6);
tic
for k = 1:bootstrap_iters
    %for n = 1:174
        bs_result = nan(n_iters,6);
        for i = 1:n_iters
            try
            Avox = bootstrap_sample(:,k);
            if (min(Avox) > 0)
                % add normally distributed random numbers to the one we used above; 
                % the scale of each random number should reflect the scale of each individual parameter to ensure large enough perturbations.
    
                %S0 = x(1)^2; diff = x(2)^2;f = exp(x(3))/(1+exp(x(3))); theta = x(4); phi = x(5);
 
                x1_start = normrnd(40,10);  % S0 = x(1)^2 in [0,4000], 
                x2_start = normrnd(0,0.5);  % diff = x(2)^2 in [0,1], so x2 in N(0,1) 
                x3_start = normrnd(0,5); % f = exp(x(3))/(1+exp(x(3))) in [0,1] so x3 in [0,6]
                x4_start = normrnd(0,5); % theta = 2*pi*exp(x(4))/(1+exp(x(4))) in [0,2*pi] so x4 in [0,6]
                x5_start = normrnd(0,5); % phi = 2*pi*exp(x(5))/(1+exp(x(5))) in [0,2*pi] so x5 in [0,6]
    
                % Inverse transformation to maintain the same starting point
                startx = [x1_start x2_start x3_start x4_start x5_start];
                % Define various options for the non-linear fitting % algorithm.

                % Now run the fitting
                [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transform',startx,h,Avox,bvals,qhat);
            
                % parameter_hat: fitted parameter values S0, d, f, theta, phi
                x = parameter_hat;
                S0 = x(1)^2;
                diff = x(2)^2;
                f = exp(x(3))/(1+exp(x(3)));
                % f + f*exp(x(3)) = exp(x(3)) => exp(x(3)) = f/(1-f) => x(3) = log(f/(1-f))
                theta = x(4);%2*pi*exp(x(4))/(1+exp(x(4)));
                phi = x(5);%2*pi*exp(x(5))/(1+exp(x(5)));
    
                bs_result(i,:) = [S0,diff,f,theta,phi,RESNORM];
            end
            end
            round_min = round(min(bs_result(:,6)),3);
            resnorm_all = bs_result(:,6); 
            %percent_minimum = length(find(round(resnorm_all,3) == round_min ))/n_iters
            min_outcomes = bs_result(find(round(bs_result(:,6),3)==round_min),:);
            opt_index = find(min_outcomes(:,2)==min(min_outcomes(:,2)));
            opt_outcomes = min_outcomes(opt_index,:);
            if size(opt_outcomes,1) ~= 0
                %S0 = opt_outcomes(1);diff= opt_outcomes(2);f= opt_outcomes(3);theta= opt_outcomes(4); phi= opt_outcomes(5);RESNORM = opt_outcomes(6);
                bootstrap_results_voxel1(k,:) = opt_outcomes(1,:);
            else
                bootstrap_results_voxel1(k,:) = [nan,nan,nan,nan,nan,nan];
            end
        end
    %end
end
toc  
% 273.436136 seconds. n_iters = 30; bootstrap_iters = 1000;
% order of parameters: [S0,diff,f,theta,phi,RESNORM]

% S0:
% 2.5% - 97.5%
sorted_S0 = sort(bootstrap_results_voxel1(:,1));
S0_025_boot2 = sorted_S0(round(0.025*bootstrap_iters)); %4.1535e+03
S0_975_boot2 = sorted_S0(round(0.025*bootstrap_iters)); %4.3606e+03

% 2 sigma  mean: 4.2591e+03 std 54.3671
S0_2sigma_left_boot2 = mean(sorted_S0)-2*std(sorted_S0); %4.1504e+03
S0_2sigma_right_boot2 = mean(sorted_S0)+2*std(sorted_S0); %4.3678e+03

figure;
hist(bootstrap_results_voxel1(:,1));
hold on;
plot([S0_025_boot2,S0_975_boot2],[270,270],"r","linewidth",2);
hold on;
plot([S0_2sigma_left_boot2,S0_2sigma_right_boot2],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Residual Bootstrap: S0")
hold off;

% d mean 0.001143558507082  std  3.6743e-05
% 2.5% - 97.5% 
sorted_d = sort(bootstrap_results_voxel1(:,2));
d_025_boot2 = sorted_d(round(0.025*bootstrap_iters)); %0.001071270409988
d_975_boot2 = sorted_d(round(0.975*bootstrap_iters)); %0.001219559196781

% 2 sigma % std 3.6743e-05
d_2sigma_left_boot2 = mean(sorted_d)-2*std(sorted_d); % 0.001070072685124
d_2sigma_right_boot2 = mean(sorted_d)+2*std(sorted_d); % 0.001217044329040


figure;
hist(bootstrap_results_voxel1(:,2));
hold on;
plot([d_025_boot2,d_975_boot2],[270,270],"r","linewidth",2);
hold on;
plot([d_2sigma_left_boot2,d_2sigma_right_boot2],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Residual Bootstrap: d")
hold off;

% f mean 0.3578 std  0.0217
% 2.5% - 97.5%
sorted_f = sort(bootstrap_results_voxel1(:,3));
f_025_boot2 = sorted_f(round(0.025*bootstrap_iters)); %0.3158
f_975_boot2 = sorted_f(round(0.025*bootstrap_iters)); %0.4000

% 2 sigma 
f_2sigma_left_boot2 = mean(sorted_f)-2*std(sorted_f); %0.3144
f_2sigma_right_boot2 = mean(sorted_f)+2*std(sorted_f); %0.4012

figure;
hist(bootstrap_results_voxel1(:,3));
hold on;
plot([f_025_boot2,f_975_boot2],[280,280],"r","linewidth",2);
hold on;
plot([f_2sigma_left_boot2,f_2sigma_right_boot2],[290,290],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Residual Bootstrap: f")
hold off;

% ============================ wild bootstrap ===========================
% step one: obtain noises (already done in the residual bootstrap)
Avox = dwis(:,92,65,72); % 108*1

% step 1: fit the model (several iterations for global minimiser) 
n_iters = 100;

bs_result = nan(n_iters,6);
Avox = dwis(:,92,65,72);
h=optimset('Display','off',...
           'MaxFunEvals',20000,...
           'Algorithm','quasi-newton',... 
           'TolX',1e-10,... 
           'TolFun',1e-10);

for i = 1:n_iters
    try
    % add normally distributed random numbers to the one we used above; 
    % the scale of each random number should reflect the scale of each individual parameter to ensure large enough perturbations.
    
    %S0 = x(1)^2; diff = x(2)^2;f = exp(x(3))/(1+exp(x(3))); theta = x(4); phi = x(5);
 
    x1_start = normrnd(40,10); % S0 = x(1)^2 in [0,4000], 
    x2_start = normrnd(0,0.5); % diff = x(2)^2 in [0,1], so x2 in N(0,1) 
    x3_start = normrnd(0,5); % f = exp(x(3))/(1+exp(x(3))) in [0,1] so x3 in [0,6]
    x4_start = normrnd(0,5); % theta = 2*pi*exp(x(4))/(1+exp(x(4))) in [0,2*pi] so x4 in [0,6]
    x5_start = normrnd(0,5); % phi = 2*pi*exp(x(5))/(1+exp(x(5))) in [0,2*pi] so x5 in [0,6]
    
    % Inverse transformation to maintain the same starting point
    startx = [x1_start x2_start x3_start x4_start x5_start];
    % Define various options for the non-linear fitting % algorithm.

    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transform',startx,h,Avox,bvals,qhat);

    % parameter_hat: fitted parameter values S0, d, f, theta, phi
    x = parameter_hat;
    S0 = x(1)^2;
    diff = x(2)^2;
    f = exp(x(3))/(1+exp(x(3)));
    % f + f*exp(x(3)) = exp(x(3)) => exp(x(3)) = f/(1-f) => x(3) = log(f/(1-f))
    theta =x(4); %2*pi*exp(x(4))/(1+exp(x(4)));
    phi = x(5);%2*pi*exp(x(5))/(1+exp(x(5)));

    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
    fitted = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(- bvals*diff));
    
    bs_result(i,:) = [S0,diff,f,theta,phi,RESNORM];
    end
end

round_min = round(min(bs_result(:,6)),3);
min_outcomes = bs_result(find(round(bs_result(:,6),3)==round_min),:);
opt_index = find(min_outcomes(:,2)==min(min_outcomes(:,2)));
opt_outcomes = min_outcomes(opt_index,:);


S0 = opt_outcomes(1);diff= opt_outcomes(2);f= opt_outcomes(3);theta= opt_outcomes(4); phi= opt_outcomes(5);

fibdir_opt = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad_opt = sum(qhat.*repmat(fibdir_opt, [length(qhat) 1])');
fitted_opt = S0*(f*exp(-bvals*diff.*(fibdotgrad_opt.^2)) + (1-f)*exp(- bvals*diff));

noises = fitted_opt - Avox'; % 1*108

% step 2: use noises obtained before to form new samples and estimate x_hat
K = size(fitted_opt,2);
bootstrap_iters = 1000;
bootstrap_sample = nan(K,bootstrap_iters);

n_iters = 30;
for i = 1:bootstrap_iters
    for j =1:K 
        bootstrap_sample(j,i) = fitted_opt(j)+ noises(j)*normrnd(0,1);
        if isnan(bootstrap_sample(j,i))
            bootstrap_sample(j,i) = fitted_opt(j)+ noises(ceil(rand*K))*normrnd(0,1);
        end
    end
end

% step 3: obtain parameter estimations from bootstrap samples
bootstrap_results_voxel1 = nan(bootstrap_iters,6);
tic
rng(20071865);
for k = 1:bootstrap_iters
    %for n = 1:174
        bs_result = nan(n_iters,6);
        for i = 1:n_iters
            try
            Avox = bootstrap_sample(:,k);
            if (min(Avox) > 0)
                % add normally distributed random numbers to the one we used above; 
                % the scale of each random number should reflect the scale of each individual parameter to ensure large enough perturbations.
    
                %S0 = x(1)^2; diff = x(2)^2;f = exp(x(3))/(1+exp(x(3))); theta = x(4); phi = x(5);
 
                x1_start = normrnd(60,10);  % S0 = x(1)^2 in [0,4000], 
                x2_start = normrnd(0,0.5);  % diff = x(2)^2 in [0,1], so x2 in N(0,1) 
                x3_start = normrnd(0,5); % f = exp(x(3))/(1+exp(x(3))) in [0,1] so x3 in [0,6]
                x4_start = normrnd(0,5); % theta = 2*pi*exp(x(4))/(1+exp(x(4))) in [0,2*pi] so x4 in [0,6]
                x5_start = normrnd(0,5); % phi = 2*pi*exp(x(5))/(1+exp(x(5))) in [0,2*pi] so x5 in [0,6]
    
                % Inverse transformation to maintain the same starting point
                startx = [x1_start x2_start x3_start x4_start x5_start];
                % Define various options for the non-linear fitting % algorithm.

                % Now run the fitting
                [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_transform',startx,h,Avox,bvals,qhat);
            
                % parameter_hat: fitted parameter values S0, d, f, theta, phi
                x = parameter_hat;
                S0 = x(1)^2;
                diff = x(2)^2;
                f = exp(x(3))/(1+exp(x(3)));
                % f + f*exp(x(3)) = exp(x(3)) => exp(x(3)) = f/(1-f) => x(3) = log(f/(1-f))
                theta = x(4);%2*pi*exp(x(4))/(1+exp(x(4)));
                phi = x(5);%2*pi*exp(x(5))/(1+exp(x(5)));
    
                bs_result(i,:) = [S0,diff,f,theta,phi,RESNORM];
            end
            end
        end
        resnorm_all = bs_result(:,6);
        min_res = min(resnorm_all);
        percent_minimum = length(find(abs(resnorm_all-min_res)<0.0001))/n_iters;
        %percent_voxel(:,m,n) = [m,n,percent_minimum];

        opt_outcomes = bs_result(find(abs(resnorm_all-min_res)<0.0001),:);

        if size(opt_outcomes,1) ~= 0
            %S0 = opt_outcomes(1);diff= opt_outcomes(2);f= opt_outcomes(3);theta= opt_outcomes(4); phi= opt_outcomes(5);RESNORM = opt_outcomes(6);
            bootstrap_results_voxel1(k,:) = opt_outcomes(1,:);%[S0,diff,f,theta,phi,RESNORM];
        else
            bootstrap_results_voxel1(k,:) = [nan,nan,nan,nan,nan,nan];
        end
        
end
toc  % 308.540427 seconds.

% S0: mean 4.2556e+03 std 55.6637
% 2.5% - 97.5%
sorted_S0 = sort(bootstrap_results_voxel1(:,1));
S0_025_boot3 = sorted_S0(0.025*bootstrap_iters); %4.1449e+03
S0_975_boot3 = sorted_S0(0.975*bootstrap_iters); %4.3646e+03

% 2 sigma 
S0_2sigma_left_boot3 = mean(sorted_S0)+2*std(sorted_S0); %4.3658e+03
S0_2sigma_right_boot3 = mean(sorted_S0)-2*std(sorted_S0); % 4.1453e+03

figure;
hist(bootstrap_results_voxel1(:,1));
hold on;
plot([S0_025_boot3,S0_975_boot3],[270,270],"r","linewidth",2);
hold on;
plot([S0_2sigma_left_boot3,S0_2sigma_right_boot3],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Wild Bootstrap: S0")
hold off;

% d mean 0.001143734209097 (roughtly, to little) std 2.9203e-05
% 2.5% - 97.5%
sorted_d = sort(bootstrap_results_voxel1(:,2));
d_025_boot3 = sorted_d(0.025*bootstrap_iters); % 0.001090580970035
d_975_boot3 = sorted_d(0.975*bootstrap_iters); % 0.001201666572979

% 2 sigma 
d_2sigma_left_boot3 = mean(sorted_d)-2*std(sorted_d); %0.001085327823738
d_2sigma_right_boot3 = mean(sorted_d)+2*std(sorted_d); % 0.001202140594457

figure;
hist(bootstrap_results_voxel1(:,2));
hold on;
plot([d_025_boot3,d_975_boot3],[270,270],"r","linewidth",2);
hold on;
plot([d_2sigma_left_boot3,d_2sigma_right_boot3],[280,280],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Wild Bootstrap: d")
hold off;

% f :mean 0.3580 std 0.0224
% 2.5% - 97.5%
sorted_f = sort(bootstrap_results_voxel1(:,3));
f_025_boot3 = sorted_f(0.025*bootstrap_iters); %0.3114
f_975_boot3 = sorted_f(0.975*bootstrap_iters); %0.4000

% 2 sigma % 
f_2sigma_left_boot3 = mean(sorted_f)-2*std(sorted_f); % 0.3138
f_2sigma_right_boot3 = mean(sorted_f)+2*std(sorted_f); %0.4023

figure;
hist(bootstrap_results_voxel1(:,3));
hold on;
plot([f_025_boot3,f_975_boot3],[280,280],"r","linewidth",2);
hold on;
plot([f_2sigma_left_boot3,f_2sigma_right_boot3],[290,290],"g","linewidth",2);
legend("P(x|A)","95%","2 sigma")
title("Wild Bootstrap: f")
hold off;

% Comparisons 
figure;
subplot(1,3,1);
plot([S0_2sigma_left_boot1,S0_2sigma_right_boot1],[50,50],"r","linewidth",2);
hold on;
plot([S0_2sigma_left_boot2,S0_2sigma_right_boot2],[40,40],"g","linewidth",2);
hold on;
plot([S0_2sigma_left_boot3,S0_2sigma_right_boot3],[30,30],"b","linewidth",2);
hold on;
plot([S0_2sigma_left_MCMC,S0_2sigma_right_MCMC],[20,20],"k","linewidth",2);
hold on;
plot([S0_2sigma_left_lap,S0_2sigma_right_lap],[10,10],"c","linewidth",2);
legend("Parametric Bootstrap","Residual Bootstrap","Wild Bootstrap","MCMC","Laplace")
ylim([-10,60])
title("2-sigma range of SO")

subplot(1,3,2);

plot([d_2sigma_left_boot1,d_2sigma_right_boot1],[50,50],"r","linewidth",2);
hold on;
plot([d_2sigma_left_boot2,d_2sigma_right_boot2],[40,40],"g","linewidth",2);
hold on;
plot([d_2sigma_left_boot3,d_2sigma_right_boot3],[30,30],"b","linewidth",2);
hold on;
plot([diff_2sigma_left_MCMC,diff_2sigma_right_MCMC],[20,20],"k","linewidth",2);
hold on;
plot([diff_2sigma_left_lap,diff_2sigma_right_lap],[10,10],"c","linewidth",2);
legend("Parametric Bootstrap","Residual Bootstrap","Wild Bootstrap","MCMC","Laplace")
ylim([-10,60])
title("2-sigma range of diff")

subplot(1,3,3);
plot([f_2sigma_left_boot1,f_2sigma_right_boot1],[50,50],"r","linewidth",2);
hold on;
plot([f_2sigma_left_boot2,f_2sigma_right_boot2],[40,40],"g","linewidth",2);
hold on;
plot([f_2sigma_left_boot3,f_2sigma_right_boot3],[30,30],"b","linewidth",2);
hold on;
plot([f_2sigma_left_MCMC,f_2sigma_right_MCMC],[20,20],"k","linewidth",2);
hold on;
plot([f_2sigma_left_lap,f_2sigma_right_lap],[10,10],"c","linewidth",2);
legend("Parametric Bootstrap","Residual Bootstrap","Wild Bootstrap","MCMC","Laplace")
ylim([-10,60])
title("2-sigma range of f")
hold off;
