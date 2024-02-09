%% Dispersion calculation in a stretched soft strip
% Strip with rectangular cross-section b x h.
%
% Acoustoelasticity describes the effect of pre-stress on the wave propagation,
% which we account for using a compressible Mooney-Rivlin hyperelastic material
% model. Viscoelastic losses are included with a fractional Kelvin-Voigt model.
% The difficulty relies in the interdependence of both effects.
% 
% Depends on the DMSUITE toolbox by Weideman and Reddy:
% https://mathworks.com/matlabcentral/fileexchange/29-dmsuite
% 
% see also: 
% A. Delory, D. A. Kiefer, M. Lanoy, A. Eddi, C. Prada, and F. Lemoult,
% “Viscoelastic dynamics of a soft strip subject to a large deformation,” 
% Soft Matter, Jan. 2024, doi: 10.1039/D3SM01485A.
%
% 2023 - Alexandre Delory and Daniel A. Kiefer
% Institut Langevin, ESPCI Paris | PSL, France

% Select stretch ratios and frequencies
lambdas1 = 1:.5:2;                   % list of axial stretch ratios
lambdas3 = lambdas1.^(-.5);          % uniaxial tension
colors   = parula(length(lambdas1)); % for plotting
FREQ     = 5:5:200;                  % list of frequencies

% Global parameters
h    = 3e-3;     % thickness in m (ey-dimension)
b    = 40e-3;    % width in m (ez-dimension)
N    = 6;        % discretization in y: polynomial order of interpolants
P    = 10;       % discretization in z: polynomial order of interpolants
rho  = 1070;     % density in kg/m^3
vT   = sqrt(26); % sqrt(mu0/rho) in m/s
tau  = 260e-6;   % rheological parameter  in s (relaxation time)
n    = 0.27;     % rheological parameter (fractional derivative order)

% Additional parameters
alpha     = .29; % for MooneyRivlin: C01/(C10+C01)
betaprime = .29; % for viscoelastic model: beta=betaprime.mu0.tau^n, and nu = (1-betaprime).mu0.tau^n

% allocate memory to save results:
Nmode = 45;      % number of modes to keep
KS    = nan(length(lambdas1),Nmode*length(FREQ));
KA    = nan(length(lambdas1),Nmode*length(FREQ));
TS    = nan(length(lambdas1),Nmode*length(FREQ));
TA    = nan(length(lambdas1),Nmode*length(FREQ));

% Start fig
FIG = figure('Position',[100 100 1000 400],'Units','pixels');
subplot(121), hold on, grid on, box on, axis([-50 200 0 200]), title('Free Strip - Antisymmetric')
subplot(122), hold on, grid on, box on, axis([-50 200 0 200]), title('Free Strip - Symmetric')

% For each value of lambda, derive the equivalent elasticity tensor and compute guided waves using a Spectral Collocation Method
for i_lambda=1:length(lambdas1)
   
   lambda1=lambdas1(i_lambda);
   lambda3=lambdas3(i_lambda);
   
   ks=nan(length(FREQ),Nmode);
   ka=nan(length(FREQ),Nmode);
   ts=nan(length(FREQ),Nmode);
   ta=nan(length(FREQ),Nmode);
   for freq=FREQ
      % derive the equivalent elasticity tensor
      c = customTensor(lambda1,lambda3,vT,freq,tau,n,betaprime,alpha);
      
      % renormalization of the coefficients
      c0   = c(1,1,1,1);
      h0   = h/(lambda1*lambda3); % ie h*lambda2
      b0   = b*lambda3;
      rho0 = rho;
      fh0  = sqrt(c0/rho0);
      rhon = rho/rho0;
      cn   = c/c0;
      hn   = h0/h0;
      bn   = b0/h0;
      % relevant material matrices:
      udof = 1:3; % use a coupled Lamb+SH polarization
      cxx  = squeeze(cn(1,udof,udof,1));
      cxy  = squeeze(cn(1,udof,udof,2));
      cxz  = squeeze(cn(1,udof,udof,3));
      cyx  = squeeze(cn(2,udof,udof,1));
      cyy  = squeeze(cn(2,udof,udof,2));
      cyz  = squeeze(cn(2,udof,udof,3));
      czx  = squeeze(cn(3,udof,udof,1));
      czy  = squeeze(cn(3,udof,udof,2));
      czz  = squeeze(cn(3,udof,udof,3));
      I    = eye(size(cxx));
      % discretization
      [~, Dy_dash] = chebdif(N, 2);
      D10 = -2/hn*Dy_dash(:,:,1);  % differentiation on normalized domain
      D20 = 4/hn^2*Dy_dash(:,:,2); % second order derivative
      [~, Dz_dash] = chebdif(P, 2);
      D01 = -2/bn*Dz_dash(:,:,1);  % differentiation on normalized domain
      D02 = 4/bn^2*Dz_dash(:,:,2); % second order derivative
      % diff-matrices on 2d-grid:
      Dyz = kron(D01,D10); % = Dzy
      Dy  = kron(eye(P),D10);
      Dyy = kron(eye(P),D20);
      Dz  = kron(D01,eye(N));
      Dzz = kron(D02,eye(N));
      Id  = eye(N*P);
      % define wave operators:
      L2  = kron(cxx, Id);
      L1  = kron(cyx + cxy, Dy) + kron(czx + cxz, Dz);
      L0  = kron(cyy, Dyy) + kron(czy, Dyz) + kron(cyz, Dyz) + kron(czz, Dzz);
      M   = kron(rhon*I, Id);
      % define boundary operators:
      % dofBC = kron([1 P], [1 N]);
      Ind   = reshape(1:length(Dy), N, P); % node enumeration
      dofty = [Ind(1,:), Ind(N,:)]; % where to impose traction ty = ey.T = 0
      doftz = [Ind(:,1); Ind(:,P)]; % where to impose traction tz = ez.T = 0
      % boundary operators for traction ty:
      By1   = kron(cyx, Id(dofty, :));
      By0   = kron(cyy, Dy(dofty, :)) + kron(cyz, Dz(dofty, :));
      % boundary operators for traction tz:
      Bz1   = kron(czx, Id(doftz, :));
      Bz0   = kron(czy, Dy(doftz, :)) + kron(czz, Dz(doftz, :));
      % incorporate BCs:
      doftyg = (N*P)*(0:length(udof)-1) + dofty.'; % boundary-dofs for all three displ. ux,uy,uz
      doftzg = (N*P)*(0:length(udof)-1) + doftz;
      L2(doftyg, :) = 0; L1(doftyg, :) = By1; L0(doftyg, :) = By0; M(doftyg, :) = 0;
      L2(doftzg, :) = 0; L1(doftzg, :) = Bz1; L0(doftzg, :) = Bz0; M(doftzg, :) = 0;
      
      % solve for frequency, keep only the real part of wavenumbers
      fh       = freq*h0;
      whn      = 2*pi*fh/fh0;
      [un,khi] = polyeig(L0+whn^2*M,L1,L2);
      
      % First sort : get rid of high-wavenumber solutions and of purely evanescent solutions
      sel = logical( (abs(angle(-1i*khi))<.49*pi).*(abs(real(-1i*khi/h0))<250) );
      un    = un(:,sel);
      khi   = khi(sel);
      
      % Second sort : get rid of out-of-plane strip modes by comparing total energy in each displacement components
      [~,I] = max([sum(abs(un(1:N*P,:)).^2) ; sum(abs(un((1+N*P):2*N*P,:)).^2) ; sum(abs(un((1+2*N*P):3*N*P,:)).^2)]);
      % I==2 correspond to out-of-plane displacement (direction of thickness)
      khi   = [khi(I==1); khi(I==3)];
      un    = [un(:,I==1) un(:,I==3)];
      symm  = zeros(length(khi),1);
      
      % Third sort : upon symmetry by extracting (y,z) maps of each component of the displacement
      for imode=1:length(khi)
         ux          = reshape(squeeze(un(1:N*P,imode)),[N P]);
         uz          = reshape(squeeze(un((1+2*N*P):3*N*P,imode)),[N P]);
         symm(imode) = (sum(abs(sum(ux)+fliplr(sum(ux))).^2) > sum(abs(sum(ux)-fliplr(sum(ux))).^2));
      end
      
      khis = khi(symm==1); khis_bis = khis(~isnan(khis)); Ns = min(length(khis_bis),Nmode);
      khia = khi(symm==0); khia_bis = khia(~isnan(khia)); Na = min(length(khia_bis),Nmode);
      ks(FREQ==freq,1:Ns) =       real(-1i*khis_bis(1:Ns))/h0;
      ka(FREQ==freq,1:Na) =       real(-1i*khia_bis(1:Na))/h0;
      ts(FREQ==freq,1:Ns) = 1-abs(imag(-1i*khis_bis(1:Ns)))/h0/(100);
      ta(FREQ==freq,1:Na) = 1-abs(imag(-1i*khia_bis(1:Na)))/h0/(100);
   end
   ts(ts<0)=0;
   ta(ta<0)=0;
   
   % plot
   subplot(122)
   scatter( ks(:),repmat(FREQ,1,Nmode),15,colors(i_lambda,:),'filled','AlphaData',ts(:).^2,'MarkerFaceAlpha','flat')
   scatter(-ks(:),repmat(FREQ,1,Nmode),15,colors(i_lambda,:),         'AlphaData',ts(:).^2,'MarkerEdgeAlpha','flat')
   subplot(121)
   scatter( ka(:),repmat(FREQ,1,Nmode),15,colors(i_lambda,:),'filled','AlphaData',ta(:).^2,'MarkerFaceAlpha','flat')
   scatter(-ka(:),repmat(FREQ,1,Nmode),15,colors(i_lambda,:),         'AlphaData',ta(:).^2,'MarkerEdgeAlpha','flat')
   drawnow
   
   KS(i_lambda,:)=ks(:);
   KA(i_lambda,:)=ka(:);
   TS(i_lambda,:)=ts(:);
   TA(i_lambda,:)=ta(:);
end
