clc;clear;
tic;
warning off all;
seed        = 47408;
stream      = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(stream);

%% Data for BQ by Campbell and Yogo(2006)
global DFGLSbelt DFbelt tconf tquantL tquantR Qconf Qconf1;
DFGLSbelt   = importdata('DFGLSbelt.txt');
DFbelt      = importdata('DFbelt.txt');
Qconf       = importdata('Qconf.txt');
Qconf1      = importdata('Qconf1.txt');
tconf       = importdata('tconf.txt');
tquantL     = importdata('tquantL.txt');
tquantR     = importdata('tquantR.txt');


data_covid = readtable('covid_rate_confirmed_cases_columns_reordered_March_2021.csv');

sig = 0.05
cv = 1.96

deaths = []
for i=1:26
y = data_covid(:,i*4-1:i*4+2);
    y = table2array(y) ;
    r=find(y(:,1)==0);
    y(r,1)=NaN;
    %y = rmmissing(y);
    deaths = [deaths y(2:end,1)]
end

%%
ress = [];
for i=1:26
i
    y = data_covid(:,i*4-1:i*4+2);
    y = table2array(y) ;
    r=find(y(:,1)==0);
    y(r,1)=NaN;
    y = rmmissing(y);

Y = log(y(2:end,3))-log(y(1:end-1,3))-y(2:end,4);
X = log(y(2:end,1))-log(y(1:end-1,1)); % deaths
%X = X(2:end,1)-X(1:end-1,1); % d(d(d))
X = [ones(size(X,1)-1,1) X(1:end-1)];
Y = Y(2:end);%3:end

%% Usual HAC-based tests

T = size(Y,1);
beta = X\Y;
u = Y - X*beta;
g = u.*X;

S = LRV2(g,-2);
V = inv(X'*X/T)*S*inv(X'*X/T);

Cauchy_hat = sqrt(T)*beta/sqrt(V(2,2));

        %% Cauchy Estimation Robust t-test
        Cauchy_t = zeros(4,1);
        for ii = 1:4;
            part = ii*4;
            Cauchy_t(ii)   = CCH_t(X,Y,part);
        end;
        I_Cauchy_t = zeros(4,1);
        for ii=1:4;
            I_Cauchy_t(ii)     = Cauchy_t(ii)%abs(Cauchy_t(ii));% > tinv(1-sig/2, ii*4-1);
        %disp(tinv(0.95, ii*4-1));
        end;
        
        I_Cauchy_hat   = Cauchy_hat(2);%abs(Cauchy_hat(2));% > cv;
        res = [size(Y,1);I_Cauchy_t;I_Cauchy_hat]
        ress = [ress res];
end

toc;