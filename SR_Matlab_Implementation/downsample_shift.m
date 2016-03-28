function img_out = downsample_shift(img,ds,shift)
[m n l] = size(img);
if(l==3)
shift_img1 = shift_image(img(:,:,1),shift);
shift_img2 = shift_image(img(:,:,2),shift);
shift_img3 = shift_image(img(:,:,3),shift);
shift_img(:,:,1) = shift_img1;
shift_img(:,:,2) = shift_img2;
shift_img(:,:,3) = shift_img3;
else
  shift_img = shift_image(img,shift);  
end 
m_r = m/ds; n_r = n/ds;
img_out=imresize(shift_img,[m_r n_r]);
end
