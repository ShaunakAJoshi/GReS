clear
close all
clc


% comparing the cost of EB algorithm and RBF

% consider mesh made of structured patches with one grid contained into the
% other

% assume hexa27 with curved interfaces

% assume 3 avg iteration per element based

ng = 2:25;

it_avg = 3;

rh = 1;         % h_slave/h_master

costQuad8 = 0;
costQuad9 = 57;

% estimate relative number of operation according to mesh size ratio
if rh >= 1
  R_eb = 0.5*(rh^2+1);
else
  R_eb = (1/rh)^2;
end


eb_cost = @(x) (x.^2)*(it_avg*R_eb*(costQuad9+3^3));


if rh >= 1
  R_slave = 1;
  R_master = rh^2;
else
  R_slave = (1/rh)^2;
  R_master = 1;
end
nInt = 36;

costQuad9 = 13;
costQuad8 = 5;
rbf_cost = @(x) R_master*((1/3)*nInt^3 + 13*nInt^2) + R_slave*(nInt*x.^2);


figure(1)
eb = eb_cost(ng);
plot(ng,eb,'r-')

hold on

rbf = rbf_cost(ng);
plot(ng,rbf,'b-')

legend('Element Based','RBF')



