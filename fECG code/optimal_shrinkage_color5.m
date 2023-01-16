function [Y_os,eta,r_p,k] = optimal_shrinkage_color5(Y,loss,method)
%Optimal singular value shrinkage over color noise
%=========input====================
% Y : Noisy data matrix;
% loss: = 'fro', 'op', ,'op2', 'nuc', 'rank'
% k_l : lower bound for search of k
% k_h : upper bound for search of k
%=========output===================
% Y_os_out: denoised matrix
% eta_out: shrinked singular values
% r_p_out: estimated rank
% Pei-Chun Su, 09/2022

[p,n] = size(Y);
transpose = 0;
if p>n
    Y = Y';
    transpose = 1;
end
[p,n] = size(Y);
[U,s,V] = svd(Y);
s = diag(s);



u = eig(Y'*Y); u = sort(u,'descend');
lab = eig(Y*Y'); lab = sort(lab,'descend');
k = floor(p^(1/2.01));
lambda_p = lab(k+1)+ 1.1/(2^(2/3)-1) * (lab(k+1)-lab(2*(k)+1));

r_p = sum(lab>lambda_p+n^(-1/3));
fZ = createPseudoNoisemod(s,2*r_p,r_p,'i');
%fZ = createPseudoNoise(s,2*r_p,'i');
ov = lab(1:r_p);
eta = zeros(1,length(lab));
eta_a = zeros(1,length(lab));

for j = 1:r_p

    if method == "cut"
        lab(1:p) = s.^2; u(1:p) = s.^2;
        m1 = (1/(p-r_p) *sum(1./(lab((r_p+1):end)-ov(j))));
        dm1 = (1/(p-r_p) *sum(1./(lab((r_p+1):end)-ov(j)).^2));
        %m2 = (1/(n-r_p) *sum(1./(u((r_p+1):end)-ov(j))));
        %dm2 = (1/(n-r_p) *sum(1./(u((r_p+1):end)-ov(j)).^2));
        m2 = -(1-p/n)/ov(j) + m1*p/n;
        dm2 = (1-p/n)/(ov(j)^2) + dm1*p/n;
        
    elseif method == "imp"
        lab(1:p) = fZ.^2; u(1:p) = fZ.^2;
        m1 = (1/p *sum(1./(lab(1:end)-ov(j))));
        dm1 = (1/p *sum(1./(lab(1:end)-ov(j)).^2));
        %m2 = (1/n *sum(1./(u(1:end)-ov(j))));
        %dm2 = (1/n *sum(1./(u(1:end)-ov(j)).^2));
        m2 = -(1-p/n)/ov(j) + m1*p/n;
        dm2 = (1-p/n)/(ov(j)^2) + dm1*p/n;
    end
    Tau = ov(j)*m1*m2; dTau = m1*m2 + ov(j)*dm1*m2 + ov(j)*m1*dm2;
    d = 1/sqrt(ov(j)*m1*m2);
    a1 = abs(m1/(d^2*dTau)); a2 = (m2/(d^2*dTau));

    if loss == "fro"
        eta(j) = d*sqrt(a1*a2);
    elseif loss == "op"
        eta(j) = d;%*sqrt(min(a1,a2)/max(a1,a2));
    elseif loss == "op2"
        eta(j) = d*sqrt(min(a1,a2)/max(a1,a2));
    elseif loss == "nuc"
        eta(j) = abs(d*(sqrt(a1*a2)- sqrt((1-a1)*(1-a2))));
    elseif loss == "rank"
        eta(j) = s(j);
    end
    %eta_a(j) = sqrt(a1*a2);
    %if sqrt(a1*a2)<w
    %    eta(j) = 0;
    %    eta_a(j) =0;
    %end

end

Y_os = U*diag(eta)*V(:,1:p)';
if transpose
    Y_os = Y_os';
end
end

