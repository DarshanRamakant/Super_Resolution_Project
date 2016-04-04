function [y W]=get_img_params(num_imges,n1,n2,m1,m2,ref_img)

%ref_img = imread('test1.jpg');
%ref_img=imresize(ref_img,[n1 n2]);
disp('getting_image_parameters....')
M=m1*m2;
N=n1*n2;
W=zeros(num_imges*M,N);
y=zeros(num_imges*M,1);
for i =1:num_imges
   
     str = strcat('test',int2str(i));
     name = strcat(str,'.jpg');
     img = imread(name);
     img1=imresize(img,[n1 n2]);
     sk = estimate_shift(ref_img,img1);
     tmp = img';
     tmp = tmp(:);
     W((i-1)*M+1:i*M,:) = get_w_original(sk,n1,n2,m1,m2);
     y((i-1)*M+1:i*M)=tmp;
     
end


end