function w = get_w_original(sk,N1,N2,M1,M2)

gamma = 2;
M = M1*M2;
N = N1*N2;
w = zeros(M,N);
for i =1:M
    for j =1:N
        row_i = 1+fix((i-1)/M2);
        col_i = 1+rem((i-1),M2);
        v_i = [row_i;col_i];
        
        row_j = 1+fix((j-1)/N2);
        col_j = 1+rem((j-1),N2);
        v_j = [row_j;col_j];
        
        sk=sk(:);
        u_j = v_j + sk;
        
        w(i,j) = exp(-(norm(v_i-u_j)^2)/gamma^2);
    end
end
end