function [ HSIimage] = R2dTohsi(image2d,orihsi,M,MC,MB,stepszie,stepszieb)
[m,n,p] = size(orihsi);  
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
HSIimage=[];
iter = 0;
for bandnumber = 1:band
    for columnnumber = 1:column
        for   rownumber =1:row
            iter = iter + 1;
            i = rr(rownumber);
            j = cc(columnnumber);
            k = bb(bandnumber);
           HSIimage(i:1:i+M-1,j:1:j+MC-1,k:1:k+MB-1) = reshape(image2d(:,iter),M,MC,MB);
        end 
    end
end