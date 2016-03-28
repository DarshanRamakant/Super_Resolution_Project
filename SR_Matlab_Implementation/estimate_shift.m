function sh=estimate_shift(cim1,cim2)
[m n l]= size(cim1);
f=zeros([size(cim1,1) size(cim1,2)]);
if (l==3)
for c=1:3
    im1=double(cim1(:,:,c));
    im2=double(cim2(:,:,c));
    tmp=fft2(im1).*fft2(flipud(fliplr(im2)));
    f=f+fftshift(abs(ifft2(tmp)));
end
else
    im1=double(cim1);
    im2=double(cim2);
    tmp=fft2(im1).*fft2(flipud(fliplr(im2)));
    f=f+fftshift(abs(ifft2(tmp)));
end
[dummy,shid]=max(f(:));
sh(1)=(mod(shid-1,size(im1,1))+1)-size(im1,1)/2;
sh(2)=ceil((shid-1)/size(im1,1))-size(im1,2)/2;

sh(1) = -sh(1);
sh(2) = -sh(2);
end