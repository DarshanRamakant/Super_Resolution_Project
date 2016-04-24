function [y W]=new_get_img_params(imSeq, num_imges,n1,n2,m1,m2,ref_img)

%ref_img = imread('test1.jpg');
%ref_img=imresize(ref_img,[n1 n2]);
disp('getting_image_parameters....')


N1=n1;
N2=n2;
M1=m1;
M2=m2;
M=M1*M2;
N=N1*N2;
A=0.04;
r=1.0;
beta = 0.05;

W=zeros(num_imges*M,N);
y=zeros(num_imges*M,1);

ref_img1 = imread('test1.jpg');
%ref_img=imresize(ref_img1,[N1 N2]);

f=dir('test*.jpg');
files={f.name};
for k=1:numel(files)
  Im{k}=imread(files{k});
end

for i =1:num_imges
     temp = Im{i};
     tmp=temp';
     tmp = tmp(:);
     %whos tmp;
     %M
     %y((i-1)*M+1:i*M)=tmp;
        
        startRow = M *(i - 1) + 1;
        endRow = startRow + M - 1;
        
        % New system matrix entry.
       % W(startRow:endRow, :) = composeSystemMatrix(size(imSeq(:,:,k)), model.magFactor, model.psfWidth, model.motionParams(k));
        
        % Stack input LR images together.
        startRow
        endRow
        whos imSeq(:,:,i)
        LRImages(startRow:endRow) = imageToVector(imSeq(:,:,i));
        

end
delta_est = marcel_shift(Im,2);

for i =1:num_imges
     W((i-1)*M+1:i*M,:) = new_get_w_new(delta_est(i),N1,N2,M1,M2);   
end
    % Normalize the row sums to one.
    W = spdiags( sum(abs(W),2) + eps, 0, size(W,1), size(W,1) ) \ W;

%y=y';



end