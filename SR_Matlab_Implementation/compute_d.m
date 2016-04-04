function d = compute_d(N1,N2)
disp('computing d...')
N = N1*N2;
d = eye(N);
for i = 1:N
        row_i = 1+fix((i-1)/N2);
        col_i = 1+rem((i-1),N2);
        
        if((row_i)>1)
            row_j = row_i -1;
            col_j = col_i;
            j = (row_j-1)*N2+col_j;
            d(i,j) = -0.25;
        end
        if(row_i<N1)
            row_j = row_i +1;
            col_j = col_i;
            j = (row_j-1)*N2+col_j;
            d(i,j) = -0.25;
        end  
        if(col_i > 1)
            row_j = row_i;
            col_j = col_i-1;
            j = (row_j-1)*N2+col_j;
            d(i,j) = -0.25;
        end  
        if(col_i < N2)
            row_j = row_i;
            col_j = col_i+1;
            j = (row_j-1)*N2+col_j;
            d(i,j) = -0.25;
        end  
end
end