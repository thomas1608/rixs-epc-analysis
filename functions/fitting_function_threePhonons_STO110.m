

%fine alignment of scans within a few meV, and subtraction of elastic peak
clear
close all


%% load the spectra
folder = '.\';
filename = 'STO_110_Qdep_aligned';
aux = load([filename '.txt']);
qvec = [150, 127.5, 105, 82, 60];
normalize = 0;
mu = 0.30; w = 0.027; %mu=0.25 best 0.50    w=0.025 best 0.023

%% load the spectra
% folder = '.\';
% filename = 'KTO_110_Qdep_aligned';
% aux = load([filename '.txt']);
% qvec = [150, 127.5, 105, 82, 60];
% normalize = 0;
% mu = 0.475; w = 0.027; %mu=0.25 best 0.50    w=0.025 best 0.023

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
      14, ...
      -0.003,...
      0];
  
ub = [50, ...
      40, ...
      20, ...
      inf, ...
      0.003, ...
      inf];
  

totalpar = zeros(length(variableVec), length(par0));
totalerror = totalpar;
norm_residuals = zeros(size(variableVec));

indEnergyInitial = find(aux(:,end-1)>aux(1,1),1,'first');
energy0 = aux(indEnergyInitial:end,end-1);
matrix(:,1) = energy0;
noElastic = zeros(length(energy0), length(variableVec)+1);

ParMatrix = zeros(length(variableVec), length(par0));
ErrorParMatrix = zeros(length(variableVec), length(par0));

%CleanScanMatrix = zeros(length(variableVec),180);


%for ii= 1:length(variableVec)-2  
for ii=3
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


    norm_residuals(ii) = normRes;
    

    

    %use determined parameters as input for next scan
%     par0 = par;
%     par0(2) = 0;
    
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
    savefig(['fitting_STO_110_' num2str(ii) '.fig'])

end

ParMatrix110STO_newfitting = ParMatrix;
ErrorParMatrix110STO_newfitting = ErrorParMatrix;
cleanScan110STO = noElastic(:,2);
energy110STO = noElastic(:,1);

writematrix(ParMatrix110STO_newfitting)
writematrix(ErrorParMatrix110STO_newfitting)
writematrix(cleanScan110STO)
writematrix(energy110STO)



function [fitting,cleanScan,par,ci,resnorm] = completeFit_threePhonons(energy,scan,par0,variable,lb, ub, w, mu)

w1 = 0.022;    %took this values from Powerpoint slide 16 upper three images for Sto 110 direction and angle 150 deg
w2 = 0.062;
w3 = 0.100;
ms1 = 35; % ms = number of middlestates
ms2 = 35;
ms3 = 35;

%fits:

%lineshape for phonon 1
%lineshape for phonon 2
%lineshape for phonon 3

%gaussians for second overtones
%no mixed overtones

%requires the energyLoss vector and the scan vector as an input

indd = find(energy<0.05, 1, 'last');
mm = max(scan(1:indd)); 
sigma = w/(2*sqrt(2*log(2)));
norm_factor = I(1,0,0,ms1,ms2,ms3,w1,w2,w3,15,15,15);

gaussian = @(k,x) k(1)./(sqrt(2*pi).*k(3)).*exp(-0.5*((x-k(2))./k(3)).^2);
    % baseFunction = @(k,x) k(1).*exp(-0.5*((x-k(2))/(sigma)).^2);
    % baseFunction = @(k,x) pseudoVoigt([k(1),k(2),w,mu],x);
baseFunction = @(k,x) pseudoVoigtAsymmetric([k(1),k(2),w+0.001,w-0.001,mu],x);% +-0.005

% elastic
%g_el = @(k,x) baseFunction([k(1),k(2)], x);

%go = @(k,x)  gaussian([k(1),k(2),sigma*2], x);
% v = @(k,x)  k(1)./(sqrt(2*pi).*sigma*2).*exp(-0.5*((x-k(2))./sigma*2).^2);

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
%plot(energy2-par(5),,'--','color','#69CE3F','linewidth',2,'DisplayName','1 Ph1 + 1 Ph3 overtone')
%plot(energy2-par(5),,'--','color','#399F0F','linewidth',2,'DisplayName','1 Ph2 + 1 Ph3 overtone')

%plot(energy2-par(5),,'--','color','#D15124','linewidth',2,'DisplayName','3 Ph1 overtone')
%plot(energy2-par(5),,'--','color','#D19F24','linewidth',2,'DisplayName','3 Ph2 overtone')
%plot(energy2-par(5),,'--','color','#AD1689','linewidth',2,'DisplayName','3 Ph3 overtone')

plot(energy2-par(5),I111*v([1,w1+w1+w3+k5],energy2)+I210*v([1,2*w1+w2+k5],energy2)+I201*v([1,2*w1+w3+k5],energy2)+I021*v([1,2*w2+w3+k5],energy2)+I120*v([1,w1+2*w2+k5],energy2)+I102*v([1,w1+2*w3+k5],energy2)+I012*v([1,w2+2*w3+k5],energy2),'--','color','#A11937','linewidth',2,'DisplayName','all other degree 3 overtones')  %1 Ph1 +1 Ph2 +1 Ph3 overtone
%plot(energy2-par(5),,'--','color','#169DAD','linewidth',2,'DisplayName','2 Ph1 + 1 Ph2 overtone')
%plot(energy2-par(5),,'--','color','#167BAD','linewidth',2,'DisplayName','2 Ph1 + 1 Ph3 overtone')
%plot(energy2-par(5),,'--','color','#AD16A9','linewidth',2,'DisplayName','2 Ph2 + 1 Ph3 overtone')
%plot(energy2-par(5),,'--','color','#8DAD16','linewidth',2,'DisplayName','1 Ph1 + 2 Ph2 overtone')
%plot(energy2-par(5),,'--','color','#33AA7C','linewidth',2,'DisplayName','1 Ph1 + 2 Ph3 overtone')
%plot(energy2-par(5),,'--','color','#D65370','linewidth',2,'DisplayName','1 Ph2 + 2 Ph3 overtone')

plot(energy-par(5), scan-fitting, 'color',[0.7,0.7,0.7],'linewidth',1.5,'displayname','residual')

xlabel ('Energy Loss (eV)')
xlim([-0.1 0.35])
%ylim([0 10])
a = ylim;
% ylim([-10 a(2)])
plot([0,0],[a(1), a(2)],'k','linewidth',1,'handlevisibility','off');
ylabel('Intensity (a.u.)')
legend
savefig(['fitting_STO110_' num2str(variable) '.fig'])
end






%% technical stuff: complex fitting functions
function [semiVoigt] = semiVoigt(k,x)
f1 = voigtLeo([k(1);k(2);k(4);k(5)],x);
f2 = k(1).*exp(-0.5*((x-k(2))/(k(3)/(2*sqrt(2*log(2))))).^2);

ind1 = find(x<k(2),1,'last');
semiVoigt = [f1(1:ind1); f2(ind1+1:end)];

end

function y = pseudoVoigt(par,x)
%Voigt function (with convolution)
amp = par(1);
xc = par(2);
sigma = par(3)/(2*sqrt(2*log(2)));      %sigma of Gauss
Gamma = par(3); 
W = par(4);

%I want to work with column vector
if isrow(x)
    x = x';
end

gaussian =1/(sqrt(2*pi)*sigma)* exp(-0.5*((x-xc)./sigma).^2);
lorentzian = (Gamma/(2*pi))./((x-xc).^2+(Gamma/2).^2);

y = amp.*((1-W) .* gaussian + W.*lorentzian);

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

function y = convAsymmetricLorentzian(par,x)

%I want to work with column vector
if isrow(x)
    x = x';
end

%I enlarge the intevral to make sure convolution is right a extremes
amp = par(1);

sigma = par(5)/(2*sqrt(2*log(2)));                             %FWHM lorentzian
xRight = x(end)+(x(2)-x(1)):(x(2)-x(1)):max(par(2),x(end))+par(3)+10*(par(4))^2;
xRight = xRight';
% xLeft = x(1)-(x(2)-x(1)):-(x(2)-x(1)):min(par(2),x(1))-3*(par(4));
% xLeft = xLeft';
% xLeft = flip(xLeft);

% xTot = [xLeft;x;xRight];
xTot = [x;xRight];

% par = 1;
%I add more pointsto minimize numerical error: I will take back the
%original intevral later;
xTot = linspace(xTot(1),xTot(end), length(xTot)*3-2);
xTot = xTot';

xm = 0.5*(xTot(1)+xTot(end));

gaussian =1/(sqrt(2*pi)*sigma)* exp(-0.5*((xTot-xm)./sigma).^2);
asymmetricLorentz1 =  par(4).*(xTot-par(3))./((xTot.^2-par(2).^2).^2 + 4*(par(4).*xTot).^2);
% asymmLorentz = @(k,x)  asymmetricLorentz1(k3,x) + abs(asymmetricLorentz1(k,x));
g6 =  asymmetricLorentz1 + abs(asymmetricLorentz1);

% lorentzian = (w/(2*pi))./((xTot-xc).^2+(w/2).^2);

y = amp.*(conv(gaussian, g6, 'same')./max(conv(gaussian, g6, 'same')));

y = y(1:3:end);
y = y(1:end-length(xRight));

% y = y(1+length(xLeft):end-length(xRight));

end

function [intensity] = amentFunction(nf,g,omega0,gamma,detuningMatrix)
%detuning can be a vector or a mtrix
 
global Bmatrix
global flag

nmax = 150;
Bmatrix = zeros(nmax+1,1);
int = zeros(size(detuningMatrix));
intensity = int;



for n=0:1:nmax
    B1 = B(max(nf,n),min(nf,n),g);
    B1_2 = B(n,0,g);
    
    Bmatrix(n+1) = B1.*B1_2;
    
    int = int + Bmatrix(n+1)./(detuningMatrix+1j.*gamma+omega0.*(g-n));
    
end

intensity = abs(int).^2;


if max(Bmatrix(end-5:end)/max(Bmatrix))>10^(-5)
    flag = 1;
end

end

function [B]=B(a,b,g)
%Calculates FC factor B
%a,b must be integers
factor = sqrt(exp(-g)*factorial(a)*factorial(b));
B=0;
for l=0:1:b
    B = B + (-1)^a * (-g)^l * g^(a/2-b/2) / (factorial(b-l)*factorial(l)*factorial(a-b+l));
end
B = factor*B;
end









function [intensity] = I(n1,n2,n3,m1_values,m2_values,m3_values,w1,w2,w3,g1,g2,g3)
Intensity_number = 0;
delta = 0; %detuning between the energy of the electronic excitation and the incident photon energy
gamma = 0.150; %the inverse lifetime of the core hole in meV
z = delta + 1i*gamma;

for m1 = 0:m1_values

    B_nm_g_first1 = B_nm_g(max(n1,m1),min(n1,m1),g1);
    B_nm_g_sec1 = B_nm_g(max(m1,0),min(m1,0),g1);

    for m2 = 0:m2_values

        B_nm_g_first2 = B_nm_g(max(n2,m2),min(n2,m2),g2);
        B_nm_g_sec2 = B_nm_g(max(m2,0),min(m2,0),g2);

        for m3 = 0:m3_values

            summe = w1.*(m1-g1) + w2.*(m2-g2) + w3.*(m3-g3);
            B_nm_g_first3 = B_nm_g(max(n3,m3),min(n3,m3),g3);
            B_nm_g_sec3 = B_nm_g(max(m3,0),min(m3,0),g3);
           
            %Intensity_number = Intenstiy_number + (D_n1n2n3m1m2m3_g1g2g3(n1,n2,n3,m1,m2,m3,g1,g2,g3).*D_n1n2n3m1m2m3_g1g2g3(m1,m2,m3,0,0,0,g1,g2,g3))./(z-summe);
            Intensity_number = Intensity_number + (B_nm_g_first1.*B_nm_g_first2.*B_nm_g_first3.*B_nm_g_sec1.*B_nm_g_sec2.*B_nm_g_sec3)./(z-summe);
        end
    end
end

intensity = abs(Intensity_number).^2;

end



function [D_g1_g2_g3] = D_n1n2n3m1m2m3_g1g2g3(n1,n2,n3,m1,m2,m3,g1,g2,g3)

D_g1_g2_g3 = B_nm_g(max(n1,m1),min(n1,m1),g1).*B_nm_g(max(n2,m2),min(n2,m2),g2).*B_nm_g(max(n3,m3),min(n3,m3),g3);

end



function [B_n_m_of_g] = B_nm_g(n,m,g)

sum=0;
for l = 0:m
    sum = sum + ((-g)^l*(sqrt(g).^(n-m)))/(factorial(m-l)*factorial(l)*factorial(n-m+l));
end
B_n_m_of_g = (-1)^n*sqrt(exp(-g).*factorial(n)*factorial(m))*sum;

end