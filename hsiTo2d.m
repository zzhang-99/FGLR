function [image2d] = hsiTo2d (HSIimage,M,MC,MB,stepszie,stepszieb)
[m,n,p] = size(HSIimage);  
R         =   m-M+1;       
C         =   n-MC+1; 
B         =   p- MB+1;
rr        =   [1:stepszie:R];   
rr        =   [rr rr(end)+1:R]; 
cc        =   [1:stepszie:C];   
cc        =   [cc cc(end)+1:C]; 
bb        =   [1:stepszieb:B];   
bb        =   [bb bb(end)+1:B]; 
row       =   length(rr);   
column    =   length(cc);
band        =   length(bb);

image2d=[];
iter = 0;
for bandnumber = 1:band
    for columnnumber = 1:column
        for   rownumber =1:row
            iter = iter + 1;
            i = rr(rownumber);
            j = cc(columnnumber);
            k = bb(bandnumber);
            patch_reference = HSIimage(i:i+M-1,j:j+MC-1,k:k+MB-1);
            image2d(:,iter) = patch_reference(:);
        end
    end
end