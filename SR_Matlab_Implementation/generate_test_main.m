clc;
clear all;
img = imread('lena30.jpg');
img_gray = rgb2gray(img);
num_images =16;
generate_test_imgaes(img_gray,num_images);