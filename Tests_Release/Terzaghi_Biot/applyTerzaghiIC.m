function stateIn = applyTerzaghiIC(solv,mat,mesh,pL)
% Apply initial conditions for Terzaghi problem

% Get Material parameters from materials class
porosity = mat.getMaterial(1).PorousRock.getPorosity();
cf = mat.getFluid().getFluidCompressibility; % [kPa^-1] Fluid compressibility
E = mat.getMaterial(1).ConstLaw.E;
nu = mat.getMaterial(1).ConstLaw.nu;
biot = mat.getMaterial(1).PorousRock.getBiotCoefficient();

% compute material parameters for analytical formulas
lambda =(E*nu)/((1+nu)*(1-2*nu)); %[kPa] first lamè constant
G = E/(2*(1+nu)); %[kPa] second lamè constant
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(G/3) + biot^2*M;

solv.state.data.pressure = solv.state.data.pressure+(biot*M*abs(pL))/(Ku+4*G/3);
zu = mesh.coordinates(:,3);
solv.state.data.dispConv(3:3:end) = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
solv.state.data.dispCurr = solv.state.data.dispConv;
end

