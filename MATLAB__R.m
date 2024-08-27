% Initialize parameters
n = 250;
k = 2;
m = 200;
p = 0.3;
mu = 1.5;
sigma = 1;
sigma2 = 1;

% Initialize result storage variables
mhde1prop = zeros(m, 1);
mhde1sig = zeros(m, 1);
mhde1mu = zeros(m, 1);
mhde2prop = zeros(m, 1);
mhde2sig = zeros(m, 1);
mhde2mu = zeros(m, 1);
symmprop = zeros(m, 1);
symmsig = zeros(m, 1);
symmmu = zeros(m, 1);
semiemprop = zeros(m, 1);
semiemsig = zeros(m, 1);
semiemmu = zeros(m, 1);
semipiprop = zeros(m, 1);
semipisig = zeros(m, 1);
semipimu = zeros(m, 1);
semiemtrueprop = zeros(m, 1);
semiemtruesig = zeros(m, 1);
semiemtruemu = zeros(m, 1);
mhde = [];
sym = [];
numitersemi = [];

% Simulation loop
for i = 1:m
    % Generate data
    n1 = binornd(n, p);
    x1 = normrnd(0, sigma, 1, n1);
    x2 = normrnd(mu, sigma2, 1, n - n1);
    x = [x1, x2]';

    % Calculate MHDE estimations
    temp = mixonekn(x);
    mhdeest = mhdem1(x, temp, p, sigma, mu);
    mhde1prop(i) = mhdeest.pi;
    mhde1sig(i) = mhdeest.sigma;
    mhde1mu(i) = mhdeest.mu;
    mhde(i, :) = [mhdeest.initialtrue, mhdeest.numiter];

    % Calculate symmetrization estimations
    symmest = symm2(x, temp, p, sigma, mu);
    symmprop(i) = symmest.pi;
    symmsig(i) = symmest.sigma;
    symmmu(i) = symmest.mu;
    sym(i, :) = [symmest.initialtrue, symmest.numiter];

    % Calculate semi-parametric estimations
    semiest = semisong(x, temp, p, sigma, mu);
    semiemprop(i) = semiest.emprop;
    semiemsig(i) = semiest.emsig;
    semiemmu(i) = semiest.emmu;
    semiemtrueprop(i) = semiest.emtrueprop;
    semiemtruesig(i) = semiest.emtruesig;
    semiemtruemu(i) = semiest.emtruemu;
    semipiprop(i) = semiest.piprop;
    semipisig(i) = semiest.pisig;
    semipimu(i) = semiest.pimu;
    numitersemi(i, :) = [semiest.emnumiter, semiest.emtruenumiter];
end

% Calculate results for MHDE1
resmhde1.prop = [mean(mhde1prop), sqrt(var(mhde1prop)), mean((mhde1prop - p).^2)];
resmhde1.sig = [mean(mhde1sig), sqrt(var(mhde1sig)), mean((mhde1sig - sigma).^2)];
resmhde1.mu = [mean(mhde1mu), sqrt(var(mhde1mu)), mean((mhde1mu - mu).^2)];

% Calculate results for semi-parametric estimators
ressemipi.prop = [mean(semipiprop), sqrt(var(semipiprop)), mean((semipiprop - p).^2)];
ressemipi.sig = [mean(semipisig), sqrt(var(semipisig)), mean((semipisig - sigma).^2)];
ressemipi.mu = [mean(semipimu), sqrt(var(semipimu)), mean((semipimu - mu).^2)];

% Calculate results for semi-em estimators
ressemiem.prop = [mean(semiemprop), sqrt(var(semiemprop)), mean((semiemprop - p).^2)];
ressemiem.sig = [mean(semiemsig), sqrt(var(semiemsig)), mean((semiemsig - sigma).^2)];
ressemiem.mu = [mean(semiemmu), sqrt(var(semiemmu)), mean((semiemmu - mu).^2)];

% Calculate results for semi-em true estimators
ressemiemtrue.prop = [mean(semiemtrueprop), sqrt(var(semiemtrueprop)), mean((semiemtrueprop - p).^2)];
ressemiemtrue.sig = [mean(semiemtruesig), sqrt(var(semiemtruesig)), mean((semiemtruesig - sigma).^2)];
ressemiemtrue.mu = [mean(semiemtruemu), sqrt(var(semiemtruemu)), mean((semiemtruemu - mu).^2)];

% Calculate results for symmetrization
ressymm.prop = [mean(symmprop), sqrt(var(symmprop)), mean((symmprop - p).^2)];
ressymm.sig = [mean(symmsig), sqrt(var(symmsig)), mean((symmsig - sigma).^2)];
ressymm.mu = [mean(symmmu), sqrt(var(symmmu)), mean((symmmu - mu).^2)];


