function [cpg1_new, cpg2_new] = transForm(transformations, hdr1, hdr2)

% Function converts the transformation values into a 5D vector

%Split into four columns, transpose and reshape
cpg1_1 = reshape(transformations(1:1400,1:4489)', [67,67,1400,1,1]);
cpg1_2 = reshape(transformations(1:1400,4490:8978)', [67,67,1400,1,1]);
cpg2_1 = reshape(transformations(1:1400,8979:13467)', [67,67,1400,1,1]);
cpg2_2 = reshape(transformations(1:1400,13468:17956)', [67,67,1400,1,1]);

% Create empty 5D arrays for the for-loop below
cpg1reshapefinal = zeros(67,67,1400,1,2);
cpg2reshapefinal = zeros(67,67,1400,1,2);

% Place the four columns into the empty 5D arrays
for i = 1:1400
    cpg1reshapefinal(:,:,i,1,1)=cpg1_1(:,:,i);
    cpg1reshapefinal(:,:,i,1,2)=cpg1_2(:,:,i);
    cpg2reshapefinal(:,:,i,1,1)=cpg2_1(:,:,i);
    cpg2reshapefinal(:,:,i,1,2)=cpg2_2(:,:,i);
end

% Create the transformations structs

% Create structs like the old cpg1 and cpg2 for the deformation function
Z = 1;
for i = 1:1400
    cpg1_new(Z).img = cpg1reshapefinal(:,:,i,:,:);
    cpg1_new(Z).hdr = hdr1(:,i);
    cpg2_new(Z).img = cpg2reshapefinal(:,:,i,:,:);
    cpg2_new(Z).hdr = hdr2(:,i);
    Z=Z+1;
end



end 