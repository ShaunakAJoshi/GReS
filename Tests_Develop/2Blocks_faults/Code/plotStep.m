function plotStep(outStruct,tStep)
% Plot solution for a desired time step
% Figure 1: converged stresses and gap along the vertical axis
% Figure 2: convergence profiele of active set loop

figure(1)
title(strcat('Converge profile - Load step_',num2str(tStep)))
str = outStruct(tStep);
numbASIter = max(str.itAS);
k = 0;
for i = 1:numbASIter
   id = str.itAS == i;
   semilogy(k+1:k+sum(id),str.rhsNorm(id),'k-s','LineWidth',1)
   hold on
   k = k+sum(id);
end
xlabel('Non-Linear Iteration')
ylabel('Residual')


coordZ = linspace(0,15,numel(str.gap));
figure(2)
subplot(1,3,1)
plot(str.s_n,coordZ,'k-s','LineWidth',1)
xlabel('\sigma_n (kPa)')
ylabel('z (m)')
xlim([-7 0])
ylim([1 15])
%
subplot(1,3,2)
plot(str.tauNorm,coordZ,'k-s','LineWidth',1)
xlabel('\tau_norm (kPa)')
ylabel('z (m)')
xlim([0 5])
ylim([1 15])
%
subplot(1,3,3)
plot(str.gap,coordZ,'k-s','LineWidth',1)
xlabel('u_z (m)')
ylabel('z (m)')




end

