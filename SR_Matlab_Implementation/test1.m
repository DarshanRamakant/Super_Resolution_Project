clc;
clear all;
tic

sigma= 10;
lambda= 150;
d=[1 2 ;3 4];
z=[0.2 ; 0.3];
W=[1 3 ; 5 7 ];
y=[4 ; 6];
temp=((W'*(W*z-y))/(sigma^2))+((d'*(d*z))/(lambda));
t =toc