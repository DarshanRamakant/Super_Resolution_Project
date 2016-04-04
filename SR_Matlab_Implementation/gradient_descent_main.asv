clc;
clear all;
alpha =0.8;
N1=128;N2=128;
M1=32;M2=32;
num_img =16;
ref_img = imread('test1.jpg');
ref_img=imresize(ref_img,[N1 N2]);
[y W]=get_img_params(num_img,N1,N2,M1,M2,ref_img);
%d = compute_d(128,128);
ref_temp = ref_img';
z = ref_temp(:);
disp('starting gradiant descent...')
for i=1:10
    grad=compute_gradiant(z,W,y);
    z = z-alpha*grad;
    ref_img = reshape(z,N2,N1)';
    [y W]=get_img_params(num_img,N1,N2,M1,M2,ref_img);
    i
end