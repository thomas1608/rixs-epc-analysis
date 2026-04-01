%fine alignment of scans within a few meV, and subtraction of elastic peak
clear
close all


%% load the spectra
folder = '.\';
filename = 'KTO110_Qdep_529p4eV_aligned';
aux = load([filename '.txt']);
qvec = [150, 127.5, 105, 82, 60];
normalize = 0;
mu = 0.475; w = 0.023; %mu=0.25 best 0.50    w=0.025 best 0.023


%% for model with one phonon, DHO for magnon and a linear function
variableVec = qvec;


par0 = [20, ...                                %cupling constant of phonon 1 - k(1)
        10, ...                                %cupling constant of phonon 2 - k(2)
        5, ...                                 %cupling constant of phonon 3 - k(3)
        6, ...                                 %elastic peak intensity - k(4)
        0, ...                                 %elastic peak energy - k(5)
        0.05];                                 %wheight parameter for multiplying in with intensitys in fitting function - k(6)


lb = [0, ...
      0, ...
      0, ...
      0, ...
      -0.02,...
      0];

ub = [50, ...
      40, ...
      20, ...
      inf, ...
      0.02, ...
      inf];


norm_residuals = zeros(size(variableVec));

indEnergyInitial = find(aux(:,end-1)>aux(1,1),1,'first');
energy0 = aux(indEnergyInitial:end,end-1);
matrix(:,1) = energy0;
noElastic = zeros(length(energy0), length(variableVec)+1);

ParMatrix = zeros(length(variableVec), length(par0));
ErrorParMatrix = zeros(length(variableVec), length(par0));


for ii= 1:length(variableVec)-2

    varValue = variableVec(ii);
    energy = aux(:,2*ii-1);
    rixs = aux(:,2*ii)/1000;

    %normalize to dd?
    indd1 = find(energy>1,1,'first');
    indd2 = find(energy>4,1,'first');
    if normalize
        integraldd(ii) = trapz(energy(indd1:indd2), rixs(indd1:indd2));
        rixs = rixs/trapz(energy(indd1:indd2), rixs(indd1:indd2));
    else
        integraldd(ii) = trapz(energy(indd1:indd2), rixs(indd1:indd2));
    end

    %fit all spectrum with 3 peaks (+overtones)
    quasiElMin = find(energy > -0.1,1,'first');
    quasiElMax = find(energy > 0.35,1,'first');

    figure(1)
    hold on

    plot(energy, rixs,'linewidth',2)
    [fitting,cleanScan,par,ci,normRes] = completeFit_threePhonons(energy(quasiElMin:quasiElMax),rixs(quasiElMin:quasiElMax),par0,variableVec(ii), lb, ub, w, mu);

    % ---- Print EPC (g1..g3) with 95% CI, and M=omega*sqrt(g) ----
    g = par(1:3).';           % column [g1 g2 g3]
    ci_g = ci(1:3,:);         % their CIs [lower, upper]
    err_minus = g - ci_g(:,1);
    err_plus  = ci_g(:,2) - g;

    fprintf('[%6.1f°] e-ph couplings (95%% CI):\n', varValue);
    for j = 1:3
        fprintf('   g%d = %.4g  (+%.3g / -%.3g)\n', j, g(j), err_plus(j), err_minus(j));
    end

    % Same phonon energies used in completeFit_threePhonons (in eV):
    omegas = [0.018; 0.051; 0.107];   % w1,w2,w3
    M_eV   = omegas .* sqrt(max(g,0));
    M_meV  = 1000 * M_eV;

    % Propagate asymmetric errors: dM ≈ (ω/(2√g)) * dg
    eps = 1e-300;               % avoid division by zero if g ~ 0
    coef = omegas ./ (2*sqrt(max(g,eps)));
    dM_minus_meV = 1000 * coef .* err_minus;
    dM_plus_meV  = 1000 * coef .* err_plus;

    fprintf('[%6.1f°] M = ω√g (meV):\n', varValue);
    for j = 1:3
        fprintf('   M%d = %.2f  (+%.2f / -%.2f)\n', j, M_meV(j), dM_plus_meV(j), dM_minus_meV(j));
    end

    norm_residuals(ii) = normRes;

    aligned(:,2*ii-1:2*ii) = [energy-par(2),rixs];
    matrix(:,ii+1) = interp1(energy-par(2),rixs,energy0);
    noElastic(:,ii+1) = interp1(energy-par(2),[rixs(1:quasiElMin-1); cleanScan; rixs(quasiElMax+1:end)],energy0, 'linear', 0);

    ParMatrix(ii,:) = par;
    ErrorParMatrix(ii,:) = ci(:,2)-par(:);

end
noElastic(:,1) = energy0;

MeanNormRes = mean(norm_residuals);


for ii=1:6
    titles = {'cupling constant of phonon 1 vs. angle in 110 direction', ...
        'cupling constant of phonon 2 vs. angle in 110 direction','cupling constant of phonon 3 vs. angle in 110 direction', ...
        'Intensity vs. Angle for elastic peak in 110 direction','Energy vs. Angle for elastic peak in 110 direction', ...
        'wheigt parameter for multiplying in with intensitys in fitting function vs. angle in 110 direction'};
    yaxistitles = {'cupling constant g','cupling constant g','cupling constant g','Intensity (a.u.)','Energy (eV)','constant that multiplies fitting'};
    figure(ii+6)
    errorbar([150, 127.5, 105],ParMatrix(1:3,ii),ErrorParMatrix(1:3,ii),'-o','linewidth',1.5)
    title(titles{ii})
    xlabel('angle (degree)')
    ylabel(yaxistitles{ii})
    xlim([60 150])
    savefig(['fitting_KTO_110_' num2str(ii) '.fig'])

end

ParMatrix110KTO_newfitting = ParMatrix;
ErrorParMatrix110KTO_newfitting = ErrorParMatrix;
cleanScan110KTO = noElastic(:,2);
energy110KTO = noElastic(:,1);

writematrix(ParMatrix110KTO_newfitting)
writematrix(ErrorParMatrix110KTO_newfitting)
writematrix(cleanScan110KTO)
writematrix(energy110KTO)



function [fitting,cleanScan,par,ci,resnorm] = completeFit_threePhonons(energy,scan,par0,variable,lb, ub, w, mu)

w1 = 0.018;    %took this values from Powerpoint slide 16 upper three images for Sto 110 direction and angle 150 deg
w2 = 0.051;
w3 = 0.107;
ms1 = 35; % ms = number of middlestates
ms2 = 35;
ms3 = 35;

indd = find(energy<0.05, 1, 'last');
mm = max(scan(1:indd));
sigma = w/(2*sqrt(2*log(2)));
norm_factor = I(1,0,0,ms1,ms2,ms3,w1,w2,w3,15,15,15);

baseFunction = @(k,x) pseudoVoigtAsymmetric([k(1),k(2),w+0.001,w-0.001,mu],x);% +-0.005

%phonons
g = @(k,x) pseudoVoigtAsymmetric([k(1),k(2),w+0.001,w-0.001,mu],x);

%overtones
v = @(k,x)  k(1)/2.*exp(-0.5*((x-k(2))./(sigma*2)).^2);

modelFun = @(k,x) k(6)/norm_factor*(I(1,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*g([1;w1+k(5)],x) + ...% 1 ph
                        I(0,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*g([1;w2+k(5)],x) + ...
                        I(0,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*g([1;w3+k(5)],x) + ...
                        I(2,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w1+k(5)],x) + ...% 2 ph
                        I(0,2,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w2+k(5)],x) + ...
                        I(0,0,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w3+k(5)],x) + ...
                        I(1,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w1+w2+k(5)],x) + ...
                        I(1,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w1+w3+k(5)],x) + ...
                        I(0,1,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w2+w3+k(5)],x) + ...
                        I(3,0,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,3*w1+k(5)],x) + ... % 3 ph
                        I(0,3,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,3*w2+k(5)],x) + ...
                        I(0,0,3,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,3*w3+k(5)],x) + ...
                        I(1,1,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w1+w2+w3+k(5)],x) + ...
                        I(2,1,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w1+w2+k(5)],x) + ...
                        I(2,0,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w1+w3+k(5)],x) + ...
                        I(0,2,1,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,2*w2+w3+k(5)],x) + ...
                        I(1,2,0,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w1+2*w2+k(5)],x) + ...
                        I(1,0,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w1+2*w3+k(5)],x) + ...
                        I(0,1,2,ms1,ms2,ms3,w1,w2,w3,k(1),k(2),k(3)).*v([1,w2+2*w3+k(5)],x)) + ...
                        g([k(4),k(5)],x);


exitflag=0;
while exitflag==0
    [par,resnorm,residual,exitflag,outputStats,~,jacobian] = lsqcurvefit(modelFun, par0, energy, scan, lb, ub);
    par0 = par;
end


fitting = modelFun(par,energy);
ci = nlparci(par,residual,'jacobian',jacobian);

k1 = [par(1)];
k2 = [par(2)];
k3 = [par(3)];
k4 = [par(4)];
k5 = [par(5)];
k6 = [par(6)]/norm_factor;

I100 = k6*I(1,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I010 = k6*I(0,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I001 = k6*I(0,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I200 = k6*I(2,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I020 = k6*I(0,2,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I002 = k6*I(0,0,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I110 = k6*I(1,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I101 = k6*I(1,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I011 = k6*I(0,1,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I300 = k6*I(3,0,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I030 = k6*I(0,3,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I003 = k6*I(0,0,3,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I111 = k6*I(1,1,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I210 = k6*I(2,1,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I201 = k6*I(2,0,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I021 = k6*I(0,2,1,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I120 = k6*I(1,2,0,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I102 = k6*I(1,0,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);
I012 = k6*I(0,1,2,ms1,ms2,ms3,w1,w2,w3,k1,k2,k3);


%% %%%%%%%%%%%%%%% elastic subtraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleanScan = scan - g([1,k5],energy);

%% %%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
energy2 = linspace(energy(1),energy(end),5*length(energy))-par(5);
energy2 = energy2';
plot(energy-par(5),scan ,'ko-','markersize',4,'linewidth',0.7,'DisplayName',[num2str(variable), ' r.l.u.'])
hold on
plot(energy-par(5), fitting  , 'r-','linewidth',2,'displayname','Fitting')

plot(energy2-par(5),g([k4,k5],energy2),'-','color',[0.9290, 0.6940, 0.1250],'linewidth',2,'DisplayName','Elastic')
plot(energy2-par(5),I100*g([1;w1+k5],energy2),'-','color','#E45618','linewidth',2,'DisplayName','Phonon 1')
plot(energy2-par(5),I010*g([1;w2+k5],energy2),'-','color',[0.6,0.1,1],'linewidth',2,'DisplayName','Phonon 2')
plot(energy2-par(5),I001*g([1;w3+k5],energy2),'-','color',[0.3010, 0.7450, 0.9330],'linewidth',2,'DisplayName','Phonon 3')
plot(energy2-par(5),I200*v([1,2*w1+k5],energy2)+I300*v([1,3*w1+k5],energy2),'--','color','#E45618','linewidth',2,'DisplayName','2*Ph1+ 3*Ph1 overtone')
plot(energy2-par(5),I020*v([1,2*w2+k5],energy2)+I030*v([1,3*w2+k5],energy2),'--','color',[0.6,0.1,1],'linewidth',2,'DisplayName','2*Ph2+3*Ph2 overtone')
plot(energy2-par(5),I002*v([1,2*w3+k5],energy2)+I003*v([1,3*w3+k5],energy2),'--','color',[0.3010, 0.7450, 0.9330],'linewidth',2,'DisplayName','2*Ph3+3*Ph3 overtone')
plot(energy2-par(5),I110*v([1,w1+w2+k5],energy2)+I101*v([1,w1+w3+k5],energy2)+I011*v([1,w2+w3+k5],energy2),'--','color',[0.3010, 0.7450, 0.9330],'linewidth',2,'DisplayName','all other degree 2 overtones') %1 Ph1 + 1 Ph2 overtone

plot(energy2-par(5),I111*v([1,w1+w1+w3+k5],energy2)+I210*v([1,2*w1+w2+k5],energy2)+I201*v([1,2*w1+w3+k5],energy2)+I021*v([1,2*w2+w3+k5],energy2)+I120*v([1,w1+2*w2+k5],energy2)+I102*v([1,w1+2*w3+k5],energy2)+I012*v([1,w2+2*w3+k5],energy2),'--','color','#A11937','linewidth',2,'DisplayName','all other degree 3 overtones')  %1 Ph1 +1 Ph2 +1 Ph3 overtone

plot(energy-par(5), scan-fitting, 'color',[0.7,0.7,0.7],'linewidth',1.5,'displayname','residual')

xlabel ('Energy Loss (eV)')
xlim([-0.1 0.35])
a = ylim;
plot([0,0],[a(1), a(2)],'k','linewidth',1,'handlevisibility','off');
ylabel('Intensity (a.u.)')
legend
savefig(['fitting_KTO110_' num2str(variable) '.fig'])
end



function y = pseudoVoigtAsymmetric(par,x)
%Voigt function (with convolution)
amp = par(1);                            %intensity of peak (NOT integral)
xc = par(2);                             %position of zero
sigma1 = par(3)/(2*sqrt(2*log(2)));      %sigma of Gauss 1
Gamma1 = par(3);                         %FWHM of Lorentzian 1
sigma2 = par(4)/(2*sqrt(2*log(2)));      %sigma of Gauss 2
Gamma2 = par(4);                         %FWHM of Lorentzian 2
W = par(5);                              %Relative weight of Lorentzian

%I want to work with column vector
if isrow(x)
    x = x';
end

gaussian1 =exp(-0.5*((x-xc)./sigma1).^2);
lorentzian1 = (Gamma1/2).^2./((x-xc).^2+(Gamma1/2).^2);

gaussian2 =exp(-0.5*((x-xc)./sigma2).^2);
lorentzian2 = (Gamma2/2).^2./((x-xc).^2+(Gamma2/2).^2);

y1 = amp.*((1-W) .* gaussian1 + W.*lorentzian1);
y2 = amp.*((1-W) .* gaussian2 + W.*lorentzian2);

y = [y1(x<xc) ; y2(x>=xc)];

end


function intensity = I(n1,n2,n3,m1_values,m2_values,m3_values,w1,w2,w3,g1,g2,g3)

Intensity_number = 0;
delta = 0;     % detuning
gamma = 0.150; % inverse lifetime (keep exactly as you had it)
z = delta + 1i*gamma;

% ---------- PRECOMPUTE all B factors used in the loops ----------
% m>=0 always here, so B(max(m,0),min(m,0),g) == B(m,0,g)
B1_n_m = zeros(m1_values+1,1);
B1_m_0 = zeros(m1_values+1,1);
for m1 = 0:m1_values
    B1_n_m(m1+1) = B_nm_g(max(n1,m1), min(n1,m1), g1);
    B1_m_0(m1+1) = B_nm_g(m1, 0, g1);
end

B2_n_m = zeros(m2_values+1,1);
B2_m_0 = zeros(m2_values+1,1);
for m2 = 0:m2_values
    B2_n_m(m2+1) = B_nm_g(max(n2,m2), min(n2,m2), g2);
    B2_m_0(m2+1) = B_nm_g(m2, 0, g2);
end

B3_n_m = zeros(m3_values+1,1);
B3_m_0 = zeros(m3_values+1,1);
for m3 = 0:m3_values
    B3_n_m(m3+1) = B_nm_g(max(n3,m3), min(n3,m3), g3);
    B3_m_0(m3+1) = B_nm_g(m3, 0, g3);
end
% ---------------------------------------------------------------

for m1 = 0:m1_values
    B_nm_g_first1 = B1_n_m(m1+1);
    B_nm_g_sec1   = B1_m_0(m1+1);

    for m2 = 0:m2_values
        B_nm_g_first2 = B2_n_m(m2+1);
        B_nm_g_sec2   = B2_m_0(m2+1);

        for m3 = 0:m3_values
            % exactly your original summe:
            summe = w1.*(m1-g1) + w2.*(m2-g2) + w3.*(m3-g3);

            B_nm_g_first3 = B3_n_m(m3+1);
            B_nm_g_sec3   = B3_m_0(m3+1);

            Intensity_number = Intensity_number + ...
                (B_nm_g_first1.*B_nm_g_first2.*B_nm_g_first3.* ...
                 B_nm_g_sec1 .*B_nm_g_sec2 .*B_nm_g_sec3) ./ (z - summe);
        end
    end
end

intensity = abs(Intensity_number).^2;

end



function [B_n_m_of_g] = B_nm_g(n,m,g)

sum=0;
for l = 0:m
    sum = sum + ((-g)^l*(sqrt(g).^(n-m)))/(factorial(m-l)*factorial(l)*factorial(n-m+l));
end
B_n_m_of_g = (-1)^n*sqrt(exp(-g).*factorial(n)*factorial(m))*sum;

end
