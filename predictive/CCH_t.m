function [t_stat] = CCH_t(X,dY,part)

%Name:      CCH_t: Cauchy Estimation Robust t-test for Predictability Regression
%Purpose:   Get the estimates and its t-statstistics
%Format:    [b_hat, t_stat] = CCH_t(X,dY,part)
%Input:     const: regressors corresponding to constant
%           X: predictor (X_0,X_1, ..., X_N)
%           dY: return (dY_1, ..., dY_N)
%Output:    b_hat: estimates(2 by 1 vector)
%           t_stat: t-statistics(2 by 1 vector)
%
%2017-08-21             By Anton Skrobotov and Yongok Choi
%Version 1.1 

 

N           = length(X);        

XX       = X;% - cumsum(X)./(1:N)'; % de-recursive_mean for X       
YY       = dY;% - cumsum(dY)./(1:N)'; % de-recursive_mean for X

t_cauchi = zeros(part,1);
for j=1:part;
    %floor(N/part)*(j-1)+1
    %floor(N/part)*j
    if j == part
     	dYY = YY(floor(N/part)*(j-1)+1:floor(N/part)*j-1);
    	XXl = XX(floor(N/part)*(j-1)+1:floor(N/part)*j-1,:);
    else
    	dYY = YY(floor(N/part)*(j-1)+1:floor(N/part)*j);
    	XXl = XX(floor(N/part)*(j-1)+1:floor(N/part)*j,:);
    end; 
   % disp(N/part);
    
  %  XXl       = XXl - cumsum(XXl)./(1:(N/part))'; % de-recursive_mean for X       
%disp(XXl);
    bb = XXl\dYY;
    t_cauchi(j)   = bb(2);%sum(sign(XXl).*dYY);%Z'*dYY;
   % disp(t_cauchi(j));
end;

%disp(mean(t_cauchi));


dif = t_cauchi-mean(t_cauchi);
s = sum(dif.^2)/(part-1);	 
t_beta =  sqrt(part)*(mean(t_cauchi))/sqrt(s);	   

% Final Results

t_stat      = [t_beta];

end
