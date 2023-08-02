clc;
clear all;
close all;

load CTReconPhantom

figure;
imshow(p_img,[]);  
title('p')
figure;
imshow(TrueImage,[]);
title('True image');

max_iteration=10000;
max_iteration_a=100;
tol_1=1e-3;
tol_2=tol_1;
tol_a_1=1e-3;
a=zeros(256,256,9);
v=zeros(30720,1);

tau=0.001;
lam=0.009;
rho=0.0005;
a_pre=a;
v_pre=v;
mse_plot=[];
cons_plot=[];
Lag_plot=[];
Recon=iswt2(a,1,1);
constr=p-(A*Recon(:));
Lag_plot(1,1)=norm(a(:),1)+(v'*constr)+((lam/2)*(norm(constr))^2);
mse_plot(1,1)=immse(Recon,TrueImage);
cons_plot(1,1)=norm((A*Recon(:))-p)/norm(p);
for j=1:max_iteration
    a_pre_one=a_pre;
    a_pre=a;
    v_pre_one=v_pre;
    v_pre=v;
    [a,i]=inn_loop(a,tau,A,p,v,lam,max_iteration_a,tol_a_1);
    w_a=iswt2(a,1,1);
    v=v+(rho*(p-(A*w_a(:))));
    Recon=iswt2(a,1,1);
    mse_plot(1,j+1)=immse(Recon,TrueImage);
    cons_plot(1,j+1)=norm((A*Recon(:))-p)/norm(p);
    constr=p-(A*Recon(:));
    Lag_plot(1,j+1)=norm(a(:),1)+(v'*constr)+((lam/2)*(norm(constr))^2);
    if (norm(a(:)-a_pre(:))/max(1,norm(a_pre(:)-a_pre_one(:)))<tol_1)&&(norm(v-v_pre)/max(1,norm(v_pre-v_pre_one))<tol_2)
        break
    end
    fprintf('Iteration:%i\n',j);
end

fprintf('Iteration required for convergence:%i\n',j);

figure;
imshow(Recon,[]);
title("Reconstructed CT Image")
x=1:j+1;
figure;
plot(x,mse_plot,'LineWidth',2)
xlabel('iteration')
ylabel('mse')
figure;
plot(x,cons_plot,'LineWidth',2)
xlabel('iteration')
ylabel(' ||Au-p||/||p||')
figure;
plot(x,Lag_plot,'LineWidth',2)
xlabel('iteration')
ylabel('Augmented Lagrangian value')

Recon=iswt2(a,1,1);

Immed=medfilt2(Recon,[3 3]);
figure;
imshow(Immed,[])
title("Median Filtered image")

cv_p=norm((A*Recon(:))-p)/norm(p);
d_raw=norm(Recon-TrueImage,'fro')/norm(TrueImage,'fro');
mse=immse(Recon,TrueImage);
fprintf('recon image ||Au-p||/||p||:%f\n',cv_p);
fprintf('Normalized Frobenius distance between true image and reconstructed image:%f\n',d_raw);
fprintf('Mean squared error of the reconstructed image:%f\n',mse);

mse2=immse(Immed,TrueImage);
fprintf('Mean squared error of the medfilt image:%f\n',mse2);
d_med_normalized=norm(Immed-TrueImage,'fro')/norm(TrueImage,'fro');
fprintf('Normalized Frobenius distance between true image and Median filteres image:%f\n',d_med_normalized);
cv_p_med=norm((A*Immed(:))-p)/norm(p);
fprintf('median filtered image ||Au-p||/||p||:%f\n',cv_p_med);



function [a,i]=inn_loop(a_s, tau, A, p, v, mew,max_iter_a,tol_a_loop)
a_t_prev=a_s;
for i=1:max_iter_a
    a_s_prev_prev=a_t_prev;
    a_t_prev=a_s;
    W_a=iswt2(a_s,1,1);
    c=A'*(mew*((A*W_a(:))-p)-v);
    a_s=a_s-(tau*(swt2(reshape(c,[256,256]),1,1)));
    a_s=cat(3,a_s(:,:,1),wthresh(a_s(:,:,2:end),'s',tau));
    
    if norm(a_s(:)-a_t_prev(:))/max(1,norm(a_t_prev(:)-a_s_prev_prev(:)))<tol_a_loop
        break;
    end
end
a=a_s;
end






 


