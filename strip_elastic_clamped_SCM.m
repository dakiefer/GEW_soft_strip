%% Dispersion calculation in a generally anisotropic elastic strip
% Strip with rectangular cross-section b x h.
% 
% Implements Dirichlet BCs at the lateral edges by eliminating the corresponding
% dofs (doftzg) of the matrices.
%
% ATTENTION:
% - The BCs prescribed at the CORNERs of the strip affect the solution.
% - When N is even, this affects mostly the out-of-plane waves.
% - When N is odd, this affects mostly the in-plane waves.
% - For high values of N, the BCs set at the corners become unimportant but this
% "convergence" is very slow.
%
% Use a parallel pool to do fast computations (3s).
% 
% Depends on the DMSUITE toolbox by Weideman and Reddy:
% https://mathworks.com/matlabcentral/fileexchange/29-dmsuite
% 
% see also: 
% A. Delory, D. A. Kiefer, M. Lanoy, A. Eddi, C. Prada, and F. Lemoult,
% “Viscoelastic dynamics of a soft strip subject to a large deformation,” 
% Soft Matter, Jan. 2024, doi: 10.1039/D3SM01485A.
%
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France

warning('Use with care! Depending on whether the number of collocation points N and P are even or odd, only the in-plane or the out-of-plane modes are well represented. The following example is for the in-plane waves.');

h = 3e-3;   % thickness in m (ey-dimension)
b = 40e-3;  % width in m (ez-dimension)
N = 8;      % discretization in y: polynomial order of interpolants
P = 14;     % discretization in z: polynomial order of interpolants
lbd = 999928000; % first Lamé parameter 
mu = 36000;      % second Lamé parameter
rho = 1000;      % density

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = 1e-2; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;
hn = h/h0; bn = b/h0;

% relevant material matrices: 
udof = 1:3; % coupled Lamb and SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cxz = squeeze(cn(1,udof,udof,3));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
cyz = squeeze(cn(2,udof,udof,3));
czx = squeeze(cn(3,udof,udof,1));
czy = squeeze(cn(3,udof,udof,2));
czz = squeeze(cn(3,udof,udof,3));
I = eye(size(cxx)); 

%% discretization 
[~, Dy_dash] = chebdif(N, 2); 
D10 = -2/hn*Dy_dash(:,:,1); % differentiation on normalized domain
D20 = 4/hn^2*Dy_dash(:,:,2);  % second order derivative
[~, Dz_dash] = chebdif(P, 2); 
D01 = -2/bn*Dz_dash(:,:,1); % differentiation on normalized domain
D02 = 4/bn^2*Dz_dash(:,:,2);  % second order derivative

% diff matrices on 2d-grid:
Dyz = kron(D01,D10); % = Dzy
Dy  = kron(eye(P),D10);
Dyy = kron(eye(P),D20);
Dz  = kron(D01,eye(N));
Dzz = kron(D02,eye(N));
Id =  eye(N*P);

% define wave operators:
L2 = kron(cxx, Id); 
L1 = kron(cyx + cxy, Dy) + kron(czx + cxz, Dz);
L0 = kron(cyy, Dyy) + kron(czy, Dyz) + kron(cyz, Dyz) + kron(czz, Dzz);
M  = kron(rhon*I, Id);
% define boundary operators:
Ind = reshape(1:length(Dy), N, P); % node enumeration
dofty = [Ind(1,:), Ind(N,:)]; % horizontal boundaries with normal ey
doftz = [Ind(:,1); Ind(:,P)]; % vertical boundaries with normal ez (with corners)
% boundary operators for traction ty:
By1 = kron(cyx, Id(dofty, :)); 
By0 = kron(cyy, Dy(dofty, :)) + kron(cyz, Dz(dofty, :));
% boundary operators for traction tz:
Bz1 = kron(czx, Id(doftz, :)); 
Bz0 = kron(czy, Dy(doftz, :)) + kron(czz, Dz(doftz, :));

% incorporate BCs:
doftyg = (N*P)*(0:length(udof)-1) + dofty.'; % dofs for ux,uy,uz on horizontal boundaries
doftzg = (N*P)*(0:length(udof)-1) + doftz;   % dofs for ux,uy,uz on vertical boundaries
% free boundaries on horizontal boundaries:
L2(doftyg, :) = 0; L1(doftyg, :) = By1; L0(doftyg, :) = By0; M(doftyg, :) = 0;
% fix boundaries on vertical boundaries:
L2(doftzg, :) = []; L1(doftzg, :) = []; L0(doftzg, :) = []; M(doftzg, :) = [];
L2(:, doftzg) = []; L1(:, doftzg) = []; L0(:, doftzg) = []; M(:, doftzg) = [];
gdofFree = setdiff(1:3*N*P, doftzg); % used in iterative solver

%% solve for frequency and plot:
kh = linspace(1e-2, 300, 120)*h0; % wavenumber*thickness 
whn = nan(size(M, 2), length(kh));
uvec = zeros([size(whn) 3*N*P]); tic 
parfor ii = 1:length(kh)
    kh0 = kh(ii);
    [ui, wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); 
    whn(:,ii) = sqrt(wh2);
    uvec(:,ii,gdofFree) = ui.';
end
fh = real(whn/2/pi*fh0); fh(abs(fh) <= 1e-3) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);
u = reshape(uvec, [size(fh), N*P, length(udof)]);

%% plot 
% extract information on the wave's polarization:
umag = sum(conj(u).*u, 4); % square L2-norm of displacements
umag = sum(umag, 3)/(N*P); % pseudo-integral on guide cross-section 
uxz = u; uxz(:,:,:,2)=0;   % [ux, 0, uz]-displacement field
uxzmag = sum(conj(uxz).*uxz, 4); % square L2-norm of displacements
uxzmag = sum(uxzmag, 3)/(N*P);   % pseudo-integral on guide cross-section 
udev = (umag-uxzmag)./umag;      % rel. diff. between ux-uz-displ. and uy-displ.

load("comsol_dispersion_fixed.mat"); % load reference solution
figure; hold on; xlim([0, max(kh)/h0]); ylim([0, 320]); caxis([0, 1])
plot(comsol.k, comsol.f, 'o', 'Color', [.7 .7 .7])
kkh = kh.*ones(size(fh)); % expand vector into matrix for all modes
scatter(kkh(:)/h0, fh(:)/h0, 14, udev(:), 'filled')
xlabel('wavenumber $k$ in rad/m','Interpreter','latex'), 
ylabel('frequency $\omega/2\pi$ in Hz','Interpreter','latex'),
title(sprintf('soft strip $h$ = %d mm $\\times$ $b$ = %d mm', h/1e-3, b/1e-3),'Interpreter','latex')
legend({'Comsol', 'SCM'}, 'Location','southeast')
cb = colorbar;
cb.Label.Interpreter = "latex";
cb.Label.String = "in-plane $\leftarrow$ polarization $\rightarrow$ out-of-plane";
