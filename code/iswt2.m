%inverse stationary wavelet transform for 2D image
%s---signal
%order---spline order: 1 or 3
%level---decomposition level
function rs=iswt2(ds,order,level)

%get the dimension
[height, width, arrayno] = size(ds);

%get the mask
allmask = splmask(order,1,'U');
[maskno, masklen] = size(allmask);

%space allocation
%|  |    |   |...| ......|   |    |...|
%app(le)  det(le)             det(1)
%ie, first array is approx part, followed by detail part on level le,
%then detail part on level (le-1) all the way up to detail on level 1
for le=level:-1:1
    %approx part on current level
    if le==level
        att=reshape(ds(:,:,1),height,width);
    else
        att=rs;
    end
    %mask with padding zeros
    allmask = splmask(order,le,'U');
    allmask = allmask(:,end:-1:1); %flip the masks
    [maskno, masklen] = size(allmask);
    %reconstruction one level up
    %storage of wavelets on current level starts from windex+2
    windex = (level-le)*(maskno^2-1);
    rs = zeros(size(att));
    for fi=1:maskno
        %first reconstruct on Y
        tmpy = zeros(size(rs'));
        for fj=1:maskno
            if(fi==1)&&(fj==1)
                tmpim = att;
            else
                tmpim = reshape(ds(:,:,windex+maskno*(fi-1)+fj),height,width);
            end
            if mod(fj,2)==0
                es = MirrorExtY(tmpim',(masklen-1)/2,2);
            else
                es = MirrorExtY(tmpim',(masklen-1)/2,1);
            end
            tmpy = tmpy+conv2(es,allmask(fj,:)','valid');
        end
        %then reconstruct on X
        if mod(fi,2)==0
            es = MirrorExtY(tmpy',(masklen-1)/2,2);
        else
            es = MirrorExtY(tmpy',(masklen-1)/2,1);
        end
        rs = rs+conv2(es,allmask(fi,:)','valid');
    end
    %rs=rs/4;
end
%================================================%


%%=========     SCALING FACTOR CONSIDERATION    ============%%
%its purpose is to energy preserving by L2-norm coefficients
%that is: normalization should be the one satisfying Bessel inequality

%====   the downsampled version ============================%%
%the normalization factor is 1/sqrt(2)
%with respect to the mask in the refinement or wavelet equation

%====   the undecimated version ============================%%
%the normalization factor is 1/2
%with respect to the mask in the refinement or wavelet equation

%%%=========================================================%%
%please refer to p149 <<A wavelet tour of signal processing>>

%the mask for linear or cubic spline wavelets, and Harr wavelet
    function m = splmask(order, level, version)
        
        pad=zeros(1,2^(level-1)-1);
        if order == 0
            %linear spline frame masks
            if version == 'U' %undecimated
                phi=1/2*[1 pad 0 pad 1];
                psi=1/2*[1 pad 0 pad -1];
            else  %version == 'D'
                phi=1/2*[1 pad 0 pad 1];
                psi=1/2*[1 pad 0 pad -1];
            end
            m=[phi;psi];
        elseif order==1
            %linear spline frame masks
            if version == 'U' %undecimated
                phi=1/4*[1 pad 2 pad 1];       %refinement mask
                psi1=sqrt(2)/4*[1 pad 0 pad -1];    %first wavelet mask
                psi2=1/4*[-1 pad 2 pad -1];         %second wavelet mask
            else  %version == 'D'
                phi=sqrt(2)/4*[1 pad 2 pad 1];
                psi1=1/2*[1 pad 0 pad -1];
                psi2=sqrt(2)/4*[-1 pad 2 pad -1];
            end
            m=[phi;psi1;psi2];
        elseif order==3
            %cubic spline frame masks
            if version == 'U' %undecimated
                phi=[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];            %refinement mask
                psi1=[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];             %first wavelet mask
                psi2=sqrt(3)/sqrt(2)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];%second wavelet mask
                psi3=[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];             %third wavelet mask
                psi4=[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];         %fourth wavelet mask
            else
                phi=sqrt(2)*[1/16 pad 1/4 pad 3/8 pad 1/4 pad 1/16];
                psi1=sqrt(2)*[-1/8 pad -1/4 pad 0 pad 1/4 pad 1/8];
                psi2=sqrt(3)*[1/8 pad 0 pad -1/4 pad 0 pad 1/8];
                psi3=sqrt(2)*[-1/8 pad 1/4 pad 0 pad -1/4 pad 1/8];
                psi4=sqrt(2)*[1/16 pad -1/4 pad 3/8 pad -1/4 pad 1/16];
            end
            m=[phi;psi1;psi2;psi3;psi4];
        else
            error('wrong spline order');
        end
        
    end



%mirror extension for a signal or image on X axis
    function MX=MirrorExtX(M,sn,lr)
        
        [height, width] = size(M);
        
        if height == 1 %1D signal
            if(lr==1)%symmetric extension
                MX=[M(sn:-1:1) M M(end:-1:end+1-sn)];
            end
            if(lr==2)%asymmetric extension
                MX=[-M(sn:-1:1) M -M(end:-1:end+1-sn)];
            end
        else %2D image
            if(lr==1)%symmetric extension
                MX=[M(:,sn:-1:1) M M(:,end:-1:end+1-sn)];
            end
            if(lr==2)%asymmetric extension
                MX=[-M(:,sn:-1:1) M -M(:,end:-1:end+1-sn)];
            end
        end
    end
        
        
        function MY=MirrorExtY(M,sn,lr)
            %mirror extension for a signal or image on Y axis
            
            [height, width] = size(M);
            
            if height == 1 %1D signal
                if(lr==1)%symmetric extension
                    MY=[M(sn:-1:1)'; M; M(end:-1:end+1-sn)'];
                end
                if(lr==2)%asymmetric extension
                    MY=[-M(sn:-1:1)'; M; -M(end:-1:end+1-sn)'];
                end
            else %2D image
                if(lr==1)%symmetric extension
                    MY=[M(sn:-1:1,:); M; M(end:-1:end+1-sn,:)];
                end
                if(lr==2)%asymmetric extension
                    MY=[-M(sn:-1:1,:); M; -M(end:-1:end+1-sn,:)];
                end
            end
            
        end
    



end