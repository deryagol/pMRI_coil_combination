% Demo for Subspace based coil combination approach for parallel MRI

% Derya Gol Gungor
% golgungor.1@osu.edu
% The Ohio State University
% Updated: May 2014


clear all; close all;
disp('~~~ Demo for subspace based coil combination approach ~~~')

%% Load data

%Spine dataset from PULSAR
disp('Loading data ...')
load('data_spine')
DATA_full=full_kspace_data; clear full_kspace_data
%% Parameters
samp=logical(abs(DATA_full)); %sampling pattern
[M,N,K,T]=size(DATA_full);
kSize=[5 5]; % 2D Kernel size
apply_to_ACS=1;
nCalib=60; % ACS region size
alpha=[1 0.5 0]; % Exponent value of singular values

%% SoS coil combine
sos_coil=sos_kspace(DATA_full,3);
sos_coil=sos_coil/norm(sos_coil);

if apply_to_ACS
    %Mask
    mask=zeros(size(sos_coil)); mask(sos_coil>0.015*max(sos_coil(:)))=1;
    mask=imfill(mask);
    figure; imagesc(rot90(mask,3)); axis image off; title('Mask')
end

%% Subspace based coil combine (using fully sampled data)

if T>1 % For dynamic imaging, apply to average
    yc=sum(DATA_full,4)./sum(samp,4); isnan(yc)=0;
else
    yc=DATA_full;
end
[MM,NN,~]=size(yc);
if apply_to_ACS
    yc=crop(yc,[nCalib nCalib K]);
    [MM,NN,~]=size(yc);
    hamm1=hamming(MM); hamm2=hamming(NN);
    hamm=hamm1*hamm2.';
    yc=yc.*repmat(hamm,[1 1 K]);
end
tic;
%Calibration matrix
disp('Construting calibration matrix...')
Y=ptl_convmtx2(yc(:,:,1),[1 MM],[1 NN],kSize(1),kSize(2));
for i=2:K
    Y=[Y ptl_convmtx2(yc(:,:,i),[1 MM],[1 NN],kSize(1),kSize(2))]; %2D partial convolution matrix
end

%Singular value decomposition
disp('Calculating SVD ... ')
[U,S,~]=svd(full(Y),'econ');
s=diag(S);
diff=s(1)-s(end);
threshold=0.05*diff;
r=length(find(s>=threshold)); %estimated rank of Y

if apply_to_ACS
    % SoS low-resolution image
    x=ifftshift(ifftshift(ifft(ifft(yc,M,1),N,2),1),2);
    sos_coil_lr = (sum(abs(x.^2),3)).^(1/2);
    sos_coil_lr=sos_coil_lr/norm(sos_coil_lr);
end

%Coil combination
h=figure;   set(h,'Position',[134  618  1347 375]); title_tmp='SoS, '; title_tmp2='';
img_disp0=rot90(sos_coil,3);
if apply_to_ACS
    img_disp0=[];
    img_disp00=[]; img_disp2=rot90(sos_coil,3); mmap=[];
end
combined_img=zeros(M,N,length(alpha));
combined_img_acs=zeros(M,N,length(alpha));
for j=1:length(alpha); %exponent of the singular values
    fun=s.^alpha(j);
    clear aa
    for i=1:r
        aa(:,:,i)=ifftshift(ifftshift(ifft(ifft((fun(i))*reshape(U(:,i),MM,NN),M,1),N,2),1),2);
        combined_img(:,:,j)=sos(aa,3);
        combined_img(:,:,j)=combined_img(:,:,j)/norm(combined_img(:,:,j));
        
        img_disp=[img_disp0 rot90(combined_img(:,:,j),3)];
        imagesc(abs(img_disp),[0 0.4*max(abs(img_disp(:)))]); colormap gray; axis image off; title(num2str(i));  drawnow;
        %         pause(0.01)
    end
    img_disp0=[img_disp0 rot90(combined_img(:,:,j),3)];
    if apply_to_ACS
        tmp=combined_img(:,:,j);
        modulation_map=zeros(size(sos_coil));
        modulation_map(logical(mask))=sos_coil_lr(logical(mask))./tmp(logical(mask));
        modulation_map=modulation_map/norm(modulation_map);
        Iunmodulated=sos_coil;
        Iunmodulated(logical(mask))=sos_coil(logical(mask))./modulation_map(logical(mask));
        Iunmodulated=Iunmodulated/norm(Iunmodulated);
        combined_img_acs(:,:,j)=Iunmodulated;
        mmap=[mmap rot90(modulation_map,3)];
        if j==1
            title_tmp='';
        end
        title_tmp=strcat(title_tmp,'  S\alphaS_C (\alpha=',num2str(alpha(j)),') ,');
        title_tmp2=strcat(title_tmp2,' \alpha=',num2str(alpha(j)),', ');
        img_disp2=[img_disp2 rot90(combined_img_acs(:,:,j),3)];
    else
        title_tmp=strcat(title_tmp,'  S\alphaS (\alpha=',num2str(alpha(j)),') ,');
    end
end

disp('Done!')
toc;

%%  Display results
if apply_to_ACS
    title(strcat('Low resolution S\alphaS : ',title_tmp2));
    h=figure;   set(h,'Position',[ 129  136 1347  375]);
    imagesc(abs(mmap),[0 0.4*max(abs(mmap(:)))]); colormap gray; axis image off;title(strcat('Modulation map : ',title_tmp2));
    
    title_tmp=strcat('SoS,  ',title_tmp);
    h=figure;   set(h,'Position',[134  618  1347 375]);
    imagesc(abs(img_disp2),[0 0.4*max(abs(img_disp2(:)))]); colormap gray; axis image off;title(title_tmp);
else
    title(title_tmp)
end
