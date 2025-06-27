% clear
% close all
% clc


% comparing the cost of EB algorithm and RBF

% consider mesh made of structured patches with one grid contained into the
% other

% assume hexa27 with curved interfaces

% assume 3 avg iteration per element based

ng = 2:25;

el_type = 'hexa';

rh = 1/2;         % h_slave/h_master

% costQuad8 = 0;
% costQuad9 = 75;         % minimum flops required to assemble each system

% estimate relative number of operation according to mesh size ratio
if rh >= 1
  R_slave = 0.5*(rh^2+1);
else
  R_slave = (1/rh)^2;
end


c_hexa27 = 160;   % 27+2nN (dN) + 16+nN (N) + 3^3 + 9nN
c_hexa8 = 73;     % 12 (dN) + 16 (N) + 3^3 + 9nN
switch el_type
  case 'hexa'
    c_eb = c_hexa8;
    it_avg = 3;
  case 'hexa27'
    c_eb = c_hexa27;
    it_avg = 3;
end

cost_per_gp = R_slave*(it_avg*c_eb);
eb_cost = @(x) (x.^2)*cost_per_gp;


if rh >= 1
  R_master = rh^2;
  R_slave2 = 1/rh^2; % reduced cost of evaluation after support detection
else
  R_master = 1;
  R_slave2 = 1; % reduced cost of evaluation after support detection
end

nInt = 16;
nIntSupp = 9;

c_hexa8 = R_slave*4*(2*nInt);

c_hexa27 = R_slave*3*(2*nIntSupp) + R_slave2*9*(2*nInt);

switch el_type
  case 'hexa'
    N = 4;
    c_rbf = c_hexa8;
  case 'hexa27'
    c_rbf = c_hexa27;
    N = 9;
end


% assume initial support detection

rbf_cost = @(x) R_master*((1/6)*nInt^3 + N*nInt^2) + x.^2*c_rbf;



figure(1)
eb = eb_cost(ng);
plot(ng,eb,'r-')

hold on

rbf = rbf_cost(ng);
plot(ng,rbf,'b-')


xlabel('nG')
ylabel('numb. operations')

legend('Element Based','RBF')



