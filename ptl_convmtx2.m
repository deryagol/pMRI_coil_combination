function Cchop=ptl_convmtx2(H,M,N,mx,nx)
%Finds partial convolution matrix as described in Eq. 5 
% G.Harikumar, Y.Bresler, "Perfect Blind Restoration of Images Blurred by Multiple Filters: Theory
% and Efficient Algorithms", IEEE Transacions on Image Processing, Vol.8,
% No:2, 1999

%It is also called 'valid convolution matrix' in Matlab.
%It finds a convolution matrix for rows mp:ml  and columns for np:nl 
% Notice: There exist bugs in the code.. 

% Derya Gol Gungor, 
% May, 2014
% The Ohio State University

% sh=size(H); 
% if strcmp(M(2),'end')
%     ml=sh(1);
% else
%     ml=M(2);
% end
% if strcmp(N(2),'end')
%     nl=sh(2);
% else
%     nl=N(2);
% end
ml=M(2);  nl=N(2);  
mp=M(1);  np=N(1);


[mh,nh]=size(H); 

for i=1:nh;  
    tmp=H(:,i);
    C1D=convmtx(tmp,mx);
    Cc(:,:,i)=C1D(mp:ml,:);
end


Cchop=[]; 
ind=0; 

for i=1:nx
    Cc_column=[];
    for j=np:nl
        if (j-ind>0)&&(j-ind<=nh)
            tmp= Cc(:,:,j-ind);
        else
            tmp=zeros(size(Cc(:,:,1)));  
        end
        Cc_column=[Cc_column; tmp];
    end
    Cchop=[Cchop Cc_column];
    ind=ind+1; 
end