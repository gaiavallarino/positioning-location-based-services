%%% EX.05 -- CYCLE SLIPS

clear
close all
format longg

%% STEP 1 - IMPORT DATA
[newdata] = importdata('CycleSlipsData.txt', '	', 1)
idata = newdata.data;

threshold = 3.8*1e-2; % [m]
lam = 19*1e-2;

%% STEP 2 - DIFFERENCES BETWEEN CALCULATED AND APPROX

epochs = idata(:,1);
dd_obs = idata(:,2);
dd_app = idata(:,3);

plot(epochs, dd_obs)

delta_l = dd_obs - dd_app;

plot(epochs, delta_l)

%% STEP 3 - DIFFERENCES BETWEEN EPOCHS

delta_d_L = zeros(size(epochs));
i=1;

while (i < length(epochs))
    delta_d_L(i) = delta_l(i+1) - delta_l(i);
    i = i+1;
end

plot(epochs, delta_d_L)

%% STEP 4 - ITERATIVE APPROACH

cslips = zeros(100, 1);
indices = zeros(100, 1);
k = 1;

while (k < length(delta_d_L))
    if (delta_d_L(k) > threshold)
        cslips(k) = delta_d_L(k);
    end
    k=k+1;
end

%% STEP 5 - REPAIRING CYCLE SLIP

x = delta_d_L(49)/lam;
n_s=round(x);

abs_val = lam*abs(n_s-x);

p_m = dd_obs(1:48, 1);
s_m = zeros(52, 1);
p=49;
while (p < 100)
        b = p-48;
        s_m(b) = dd_obs(p)-lam*n_s;
        p=p+1;
end

vettore_finale=[p_m; s_m(1:51, 1)];

plot(1:99, vettore_finale)

app_corretto = [dd_app(1:48, 1); dd_app(50:100, 1)];
residuals = vettore_finale - dd_app(1:99);
residuals(49)=(residuals(48)+residuals(50))/2;
plot(1:99, residuals)

