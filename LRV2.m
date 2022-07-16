function [LRV] = LRV2(mU, cLag)

ct = size(mU,1); 
cp = size(mU,2); 


    
    
	% Andrews (1991), eqn (6.4) alpha2 with AR(1) for each component
	rho_upper = 1 - 1 / sqrt(ct);
    anumer = 0;
    adenom = 0;
	for i=1:cp
		y = mU(2:ct,i);  
    %    y(1,1) = 0;  
        x = mU(1:ct-1,i);
		
        rho1 = x\y;
        rho = min(rho1, rho_upper);
        res = y - rho * x;
		% ct-2: one lost for lag, one for AR parameter
		sigma2 =res'*res / (ct - 2);
		anumer = anumer + 4 * (rho * sigma2) ^ 2 / (1 - rho) ^ 8;
		adenom = adenom + sigma2 ^ 2 / (1 - rho) ^ 4;
    end
% 	if (cp != cq || mU != mV)
% 	{
% 		for (i = 0; i < cq; ++i)
% 		{
% 			y = mV[][i];  y[0][0] = 0;  x = lag0(mV[][i], 1);
% 			rho = min(double(y'x ./ (x'x)), rho_upper);
% 			% ct-2: one lost for lag, one for AR parameter
% 			sigma2 = double(sumsqrc(y - rho * x) ) / (ct - 2);
% 			anumer += 4 * (rho * sigma2) ^ 2 / (1 - rho) ^ 8;
% 			adenom += sigma2 ^ 2 / (1 - rho) ^ 4;
% 		}
% 	}
	% (6.2) S_T for quadratic spectral kernel
	st = 1.3221 * (ct * anumer / adenom) ^ 0.2;
	
	if cLag < 0
		% determine cutoff; is at kernel argument 50, k(50)~1e-5
		cLag = floor(50 * st);
    end
%	if (cLag >= ct) cLag = ct - 1;

	% (2.5) quadratic spectral kernel loop to compute HAC matrix
	hac = mU'*mU;
	pi65 = 6 * pi / 5; 
    pisqr = 25 / (12 * pi^2);

	for i = 1:cLag
		x = i / st;
		if x <= 0.0001
			kern = 1;
		else
			pi65x = pi65 * x;
			kern = (sin(pi65x) / pi65x - cos(pi65x)) * pisqr / (x^2);
        end
		hac = hac + kern * (mU(1:ct-i,:)'*mU(1+i:ct,:) + mU(1+i:ct,:)'*mU(1:ct-i,:));
    end
	LRV = hac / ct;
end