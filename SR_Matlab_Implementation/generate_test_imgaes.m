function success = generate_test_imgaes(img,num_imges)
success =0;
rand_shift = randi(20,num_imges,2);
ds = 4;
for i = 1:num_imges
   img_out = downsample_shift(img,ds,rand_shift(i,:));
   str = strcat('test',int2str(i));
   name = strcat(str,'.jpg');
   imwrite(img_out,name);
end
success =1;
end