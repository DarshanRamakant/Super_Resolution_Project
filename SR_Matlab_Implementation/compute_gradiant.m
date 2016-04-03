function grad = compute_gradiant(z,W,y)
p =16;
[M N] = size(W);

sigma = 10;
lambda = 150;
grad = zeros(N,1);
d = compute_d(128,128);
for k =1:N
sum_1=0;
for m =1:p*M
    sum_1 = sum_1+W(m,k)*(W(m,:)*z - y(m));
end
sum_1 = (1/sigma^2)*sum_1;
sum_2 = 0;
for i =1:N
    sum_2 = sum_2+d(i,k)*(d(i,:)*z);
end
sum_2 = (1/lambda)*sum_2;
grad(k) = sum_1+sum_2;
end
end