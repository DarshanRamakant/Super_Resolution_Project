function w = get_w_new(sk,N1,N2,M1,M2)
disp('getting w....')
gamma = 2.5;
M = M1*M2;
N = N1*N2;
w = zeros(M,N);
sk=sk(:);
v_j = [M1/2;M2/2];
for j =1:M
        row_j = 1+fix((j-1)/N2);
        col_j = 1+rem((j-1),N2);
        %v_j = [row_j;col_j];
        u_j = v_j +sk;
        
        offset = ceil(8+norm(sk));
        
        if(col_j > offset)
            col_i_start = col_j -offset;
        else
            col_i_start =1;
        end
         if(col_j < N2-offset)
             col_i_end = col_j+offset;
         else
             col_i_end = N2;
         end
         
        if(row_j > offset)
            row_i_start = row_j -offset;
        else
            row_i_start =1;
        end
         if(row_j < N1-offset)
             row_i_end = row_j+offset;
         else
             row_i_end = N1;
         end
            
         for row_i = row_i_start:row_i_end
             for col_i = col_i_start:col_i_end
        v_i = [row_i;col_i];
        if(row_i<1 | col_i <0)
            row_i
            col_i
        end
        i = (row_i-1)*N2+col_i;
        w(j,i) = exp(-(norm(v_i-u_j)^2)/gamma^2);
             end
         end
  %w(j,:)=w(j,:)/sum(w(j,:));
end
end