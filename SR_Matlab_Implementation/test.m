clc;
clear all;
tic
for i =1:16
    sk=randi(15,2,1);
w_new = get_w(sk,128,128,32,32);
end

t =toc