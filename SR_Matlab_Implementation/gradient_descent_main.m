clc;
clear all;
alpha =1.4;
N1=64;N2=64;
M1=16;M2=16;
num_img =16;
ref_img = imread('test1.jpg');
ref_img=imresize(ref_img,[N1 N2]);
[y W]=get_img_params(num_img,N1,N2,M1,M2,ref_img);
%d = compute_d(128,128);
ref_temp = ref_img';
z = ref_temp(:);
z=double(z);
disp('starting gradiant descent...')
for i=1:20
    grad=compute_gradiant(z,W,y,N1,N2);
    z = z-alpha*grad;
    ref_img = reshape(z,N2,N1)';
    [y W]=get_img_params(num_img,N1,N2,M1,M2,ref_img);
    i
end

z_int = uint8(z);
img_out = reshape(z_int,64,64)';
imshow(img_out)