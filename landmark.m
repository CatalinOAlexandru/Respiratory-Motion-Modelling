%{
By Qiyue Liu
%This is a whole individual code, which only depend on libraries given by
%teacher. The last several parts are for landmark estimation. But I'm not
%confident about them. Maybe there are some mistakes... I will check them.

2021.3.23

The results are:
1 is linear model. 2 is square model. 3 is pow3 model
First value for AP, second for SI

lme1=  1.679388062822844  15.278224875826789
lme2=  1.765073648229051  14.482136842195599
lme3=  1.573436758278711  10.187688073524885

aic1 =

   1.0e+03 *

   0.729801270237259   3.821000045518490


aic2 =

   1.0e+03 *

   0.801469383179333   3.748082327262624


aic3 =

   1.0e+03 *

   0.642567142948090   3.257651915442061


bic1 =

   1.0e+03 *

   0.740289725268466   3.831488500549696


bic2 =

   1.0e+03 *

   0.747533952784069   3.838732728065299


bic3 =

   1.0e+03 *

   0.754778180299672   3.845976955580903

%}
%%
 
addpath './CMBI-project/NIfTI_20140122'
addpath './CMBI-project/display'
addpath './CMBI-project/transforms'

%% surrogate

seg_nii=load_untouch_nii('./segmentation/0007_sdt.nii');
src_nii=load_untouch_nii('./images/0007.nii');

surrogate=zeros(100,1);
history2=zeros(100,1);
history3=zeros(100,1);

for n=1:100
cpg1_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region1.nii',n));
cpg2_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region2.nii',n));


tar_nii=load_untouch_nii(sprintf('./images/%04d.nii',n));

%[res_img,def_field,dis_field]=deformNiiWithCPGsSliding(cpg1_nii,cpg2_nii,seg_nii,src_nii,tar_nii);

% image=dispNiiSlice(res,'z',1);
% pause(0.1);

mark=1;
res2=tar_nii.img';
for i =1:length(res2(1,:))
    if (res2(150,i)>20)
        
        mark=(i-1)+(20-res2(150,i-1))/(res2(150,i)-res2(150,i-1));
        break
    end
end

surrogate(n)=mark;
end
plot(surrogate);

%% Q6.3 one CP only


SI=zeros(100,1);


for n=1:100
cpg1_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region1.nii',n));
%cpg2_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region2.nii',filenum));
%tar_nii=load_untouch_nii(sprintf('./images/%04d.nii',filenum));
%[res,~,~]=deformNiiWithCPGsSliding(cpg1_nii,cpg2_nii,seg_nii,src_nii,tar_nii);
SI(n)=cpg1_nii.img(44,38,1,1,2);  
end

plot(surrogate(1:100),SI,'x');


A=[ones(100,1),surrogate(1:100,:)];
x=A\SI;

y1=x(1)+50*x(2);
y2=x(1)+52*x(2);

plot(surrogate(1:100),SI,'x'); hold on
plot([50,52],[y1,y2],'-');

predict=A*x;
residual=sum((predict-SI).^2)/length(SI);

%% For whole map 
%load data

X1=zeros(100,67*67); %region1 AP
X2=zeros(100,67*67); %region1 SI
X3=zeros(100,67*67); %region2 AP
X4=zeros(100,67*67); %region2 SI

for n=1:100
    
    cpg1_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region1.nii',n));
    cpg2_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region2.nii',n));
    
    X1(n,:)=reshape(cpg1_nii.img(:,:,1,1,1),1,[]);
    X2(n,:)=reshape(cpg1_nii.img(:,:,1,1,2),1,[]);
    X3(n,:)=reshape(cpg2_nii.img(:,:,1,1,1),1,[]);
    X4(n,:)=reshape(cpg2_nii.img(:,:,1,1,2),1,[]);
 
    
end

%% linear
S=[surrogate(1:100),ones(100,1)];
C1=S\X1;
C2=S\X2;
C3=S\X3;
C4=S\X4;

C_1=[C1;C2;C3;C4];
%% for pow2
S=[surrogate(1:100).^2,surrogate(1:100),ones(100,1)];

C1=S\X1;
C2=S\X2;
C3=S\X3;
C4=S\X4;

C_2=[C1;C2;C3;C4];
%% for pow 3

S=[surrogate(1:100).^3,surrogate(1:100).^2,surrogate(1:100),ones(100,1)];

C1=S\X1;
C2=S\X2;
C3=S\X3;
C4=S\X4;

C_3=[C1;C2;C3;C4];
%% test for residual in 44,38,2
m=44;
n=38;
y=reshape(X2,100,67,67);
y=y(:,m,n);

c1=reshape(C_1(3:4,:),[],67,67);
c1=c1(:,m,n);
c2=reshape(C_2(4:6,:),[],67,67);
c2=c2(:,m,n);
c3=reshape(C_3(5:8,:),[],67,67);
c3=c3(:,m,n);

x=linspace(min(surrogate(1:100))-0.5,max(surrogate(1:100))+0.5,100)';
plot(surrogate(1:100),y,'x'); hold on
plot(x,[x,ones(100,1)]*c1,'r-'); hold on
plot(x,[x.^2,x,ones(100,1)]*c2,'g-'); hold on
plot(x,[x.^3,x.^2,x,ones(100,1)]*c3,'b-');

residual1=mean((y-[surrogate(1:100),ones(100,1)]*c1).^2)
residual2=mean((y-[surrogate(1:100).^2,surrogate(1:100),ones(100,1)]*c2).^2)
residual3=mean((y-[surrogate(1:100).^3,surrogate(1:100).^2,surrogate(1:100),ones(100,1)]*c3).^2)

%% Evaluation

%find surrogate signals

surrogate_test=zeros(1400,1);

for n=101:1500

tar_nii=load_untouch_nii(sprintf('./images/%04d.nii',n));

mark=1;
res2=tar_nii.img';
for i =1:length(res2(1,:))
    if (res2(150,i)>20)
        
        mark=(i-1)+(20-res2(150,i-1))/(res2(150,i)-res2(150,i-1));
        break
    end

end

surrogate_test(n-100)=mark;

end
plot(surrogate_test);
%% generate deformation maps with pow3 model and save them in ./eve3
tmp=load_untouch_nii('./registrations/0001_cpp_region1.nii');
tmp2=load_untouch_nii('./registrations/0001_cpp_region2.nii');


S=[surrogate_test.^3,surrogate_test.^2,surrogate_test,ones(1400,1)];
reg1_AP=S*C_3(1:4,:);
reg1_SI=S*C_3(5:8,:);
reg2_AP=S*C_3(9:12,:);
reg2_SI=S*C_3(13:16,:);

for n=101:1500
    tmp.img(:,:,1,1,1)=reshape(reg1_AP(n-100,:),[67,67]);
    tmp.img(:,:,1,1,2)=reshape(reg1_SI(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,1)=reshape(reg2_AP(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,2)=reshape(reg2_SI(n-100,:),[67,67]);
    
    save_untouch_nii(tmp, sprintf('./eve3/%04d_cpp_region1.nii',n));
    save_untouch_nii(tmp2, sprintf('./eve3/%04d_cpp_region2.nii',n));
    
end

%% generate deformation maps with pow2 model and save them in ./eve2

tmp=load_untouch_nii('./registrations/0001_cpp_region1.nii');
tmp2=load_untouch_nii('./registrations/0001_cpp_region2.nii');


S=[surrogate_test.^2,surrogate_test,ones(1400,1)];
reg1_AP=S*C_2(1:3,:);
reg1_SI=S*C_2(4:6,:);
reg2_AP=S*C_2(7:9,:);
reg2_SI=S*C_2(10:12,:);

for n=101:1500
    tmp.img(:,:,1,1,1)=reshape(reg1_AP(n-100,:),[67,67]);
    tmp.img(:,:,1,1,2)=reshape(reg1_SI(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,1)=reshape(reg2_AP(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,2)=reshape(reg2_SI(n-100,:),[67,67]);
    
    save_untouch_nii(tmp, sprintf('./eve2/%04d_cpp_region1.nii',n));
    save_untouch_nii(tmp2, sprintf('./eve2/%04d_cpp_region2.nii',n));
    
end


%% generate deformation maps with pow1 model and save them in ./eve1

tmp=load_untouch_nii('./registrations/0001_cpp_region1.nii');
tmp2=load_untouch_nii('./registrations/0001_cpp_region2.nii');


S=[surrogate_test,ones(1400,1)];
reg1_AP=S*C_1(1:2,:);
reg1_SI=S*C_1(3:4,:);
reg2_AP=S*C_1(5:6,:);
reg2_SI=S*C_1(7:8,:);

for n=101:1500
    tmp.img(:,:,1,1,1)=reshape(reg1_AP(n-100,:),[67,67]);
    tmp.img(:,:,1,1,2)=reshape(reg1_SI(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,1)=reshape(reg2_AP(n-100,:),[67,67]);
    tmp2.img(:,:,1,1,2)=reshape(reg2_SI(n-100,:),[67,67]);
    
    save_untouch_nii(tmp, sprintf('./eve1/%04d_cpp_region1.nii',n));
    save_untouch_nii(tmp2, sprintf('./eve1/%04d_cpp_region2.nii',n));
    
end

%% Evaluate with landmarks
load('./images/landmark_pos_phys.mat');

dist_nii=load_untouch_nii('./segmentation/0007_sdt.nii');
pointsGT=landmark_pos_phys(:,:,7);
target_nii=load_untouch_nii('./images/0007.nii');

for n=101:1500
cpg1_nii3=load_untouch_nii(sprintf('./eve3/%04d_cpp_region1.nii',n));
cpg2_nii3=load_untouch_nii(sprintf('./eve3/%04d_cpp_region2.nii',n));

cpg1_nii2=load_untouch_nii(sprintf('./eve2/%04d_cpp_region1.nii',n));
cpg2_nii2=load_untouch_nii(sprintf('./eve2/%04d_cpp_region2.nii',n));

cpg1_nii1=load_untouch_nii(sprintf('./eve1/%04d_cpp_region1.nii',n));
cpg2_nii1=load_untouch_nii(sprintf('./eve1/%04d_cpp_region2.nii',n));

% cpg1_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region1.nii',n));
% cpg2_nii=load_untouch_nii(sprintf('./registrations/%04d_cpp_region2.nii',n));

points_orig=landmark_pos_phys(:,:,n);

points3=transPointsWithCPGsSliding(cpg1_nii3, cpg2_nii3, dist_nii, points_orig, target_nii);
lme3(n-100,:)=mean((points3-pointsGT).^2);

points2=transPointsWithCPGsSliding(cpg1_nii2, cpg2_nii2, dist_nii, points_orig, target_nii);
lme2(n-100,:)=mean((points2-pointsGT).^2);

points1=transPointsWithCPGsSliding(cpg1_nii1, cpg2_nii1, dist_nii, points_orig, target_nii);
lme1(n-100,:)=mean((points1-pointsGT).^2);

end
%%
figure;
plot(lme1(:,1)); hold on
plot(lme1(:,2));
disp(mean(lme1));
figure;
plot(lme2(:,1)); hold on
plot(lme2(:,2));
disp(mean(lme2));
figure;
plot(lme3(:,1)); hold on
plot(lme3(:,2));
disp(mean(lme3));
%%
k=1400;

aic1=2*2+k*log(mean(lme1))
aic2=2*3+k*log(mean(lme2))
aic3=2*4+k*log(mean(lme3))

bic1=log(k)*2+k*log(mean(lme1))
bic2=log(k)*3+k*log(mean(lme1))
bic3=log(k)*4+k*log(mean(lme1))
