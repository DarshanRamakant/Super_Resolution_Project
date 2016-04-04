function w = get_w_original(sk,N1,N2,M1,M2)

gamma = 2;
M = M1*M2;
N = N1*N2;
w = zeros(M,N);
for j =1:M
    for i =1:N
        row_j = 1+fix((j-1)/N2);
        col_j = 1+rem((j-1),N2);
        sk=sk(:);
        v_j = [row_j;col_j];
        u_j = v_j +sk;
        row_i = 1+fix((i-1)/N2);
        col_i = 1+rem((i-1),N2);
        v_i = [row_i;col_i];
        w(j,i) = exp(-(norm(v_i-u_j)^2)/gamma^2);
    end
    %w(j,:) = w(j,:)/sum(w(j,:));
end
end