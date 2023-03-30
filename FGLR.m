function [output_image]= FGLR(oriData3_noise,mu,lambda,beta)
%% %%%%%%%%%%%%%%%%%%%%%%%% Construction Graph Structure %%%%%%%%%%%%%%%%%%%%%%%%%%
[M, N, p] = size(oriData3_noise);
MS = M;MC = N;MB = 1;stepsize  = 1;stepsizeb  = 1;
image2d = hsiTo2d (oriData3_noise,MS,MC,MB,stepsize,stepsizeb);
[q,o]=size(image2d);
L  = FGLR_graph(image2d); 
%% %%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%
Y = oriData3_noise;
tmpY = image2d;
YS = zeros(M, N, p); B = zeros(M, N, p);
tmpYX = zeros(q,o);tmpB = zeros(q,o);tmpS = zeros(q,o);
rho = 0.25;
sv = 10; r2=2; tol=1e-6; max_rho1 = 1e6; rho1 = 1.5;
iter = 0; 
maxit = 4; 
%% %%%%%%%%%%%%%%%%%%%%%%%% Main Function%%%%%%%%%%%%%%%%%%%%%%%%%%
while iter<maxit
    iter = iter + 1;
    %% X subproblem
    eye_lenX = eye(o);
    temp = rho * (tmpY-tmpB-tmpS)+tmpYX;
    tmpX=((rho*eye_lenX+2*mu*L)\ temp')';
    X = R2dTohsi(tmpX,oriData3_noise,MS,MC,MB,stepsize,stepsizeb);
    L = FGLR_graph(tmpX);
    %% B subproblem
    temp_B = tmpY - tmpX -tmpS + tmpYX/rho;
    B_hat = max(temp_B - lambda/rho, 0);
    tmpB = B_hat+min(temp_B + lambda/rho, 0);
    B = R2dTohsi(tmpB,oriData3_noise,MS,MC,MB,stepsize,stepsizeb);
    %% S subproblem
    temp_S=zeros(M,N,p);
    for ii =1:p
        temp_S(:,:,ii) = (Y(:,:,ii) - X(:,:,ii)-B(:,:,ii))+ YS(:,:,ii)/(rho);
        temp2 = temp_S(:,:,ii);
        if  choosvd(p,sv) ==1
            [U1, sigma, V1] = lansvd(temp2, sv, 'L');
        else
            [U1,sigma,V1] = svd(temp2,'econ');
        end
        sigma = diag(sigma);
        svp = min(length(find(sigma>beta/(rho))),r2);
        if svp<sv
            sv = min(svp + 1, p);
        else
            sv = min(svp + round(0.05*p), p);
        end
        S(:,:,ii) = U1(:, 1:svp) * diag(sigma(1:svp) - beta/(rho)) * V1(:, 1:svp)';
    end
    tmpS= hsiTo2d(S,MS,MC,MB,stepsize,stepsizeb);
    %%
    leq1 = image2d - tmpX - tmpB -tmpS;
    leq2 = Y - X - S-B;
    stopC = max(abs(leq2(:)));
%     disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ',stopALM=' num2str(stopC,'%2.3e')]); 
    C_value(iter) = stopC;
    if stopC < tol
        break;
    else
    tmpYX = tmpYX + rho * leq1;
    YS = YS + rho * leq2;
    rho  = min(max_rho1,rho*rho1);
    end
end
output_image = X;
for i=1:p
    output_image(:,:,i)=output_image(:,:,i)+ mean(mean(S(:,:,i)));
end
end