%function ReadData

load CTReconPhantom
%%%%%  A is the geometric matrix, 

figure; %%show the measurment, the corresponding vector version is p
imshow(p_img,[]);  
%title('measurement','FontSize',20)
figure; %%%% Show the true image, resolution 256 by 256
imshow(TrueImage,[]);
%title('True image','FontSize',20);


%%%%%%   Wavelet transformation  W*u,, W^T*W = I
coef = swt2(TrueImage,1,1);
%%% Soft Thresholding needs to apply on coef(:,:,2:end)

%%%%%  Inverse wavelet transformation W^T*u
newI = iswt2(coef,1,1);

figure;
imshow(newI,[]);