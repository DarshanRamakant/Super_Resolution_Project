clc;
clear all;
img = imread('lena30.jpg');
img_original = imresize(img,[128 128]);
img_gray = rgb2gray(img_original);
imwrite(img_gray,'original.jpg');
num_images =16;
generate_test_imgaes(img_gray,num_images);
