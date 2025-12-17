function plot_ER_Ca_timecourses(sim)
% plot_ER_Ca_timecourses
%
% Plot Ca2+ concentration vs time for cell body and each satellite.
%
% Input:
%   sim : struct from simulate_ER_Ca_multiball

t      = sim.t_sol;
c_body = sim.C_sol(:,1);
c_sat  = sim.C_sol(:,2:end);
Nsat   = sim.geom.Nsat;

figure; hold on;
plot(t, c_body, 'k-', 'LineWidth', 2);

colors = lines(Nsat);
for i = 1:Nsat
    plot(t, c_sat(:,i), 'LineWidth', 1.5, 'Color', colors(i,:));
end

xlabel('Time (s)');
ylabel('Ca^{2+} concentration (mM)');

labels = [{'Cell body'}, arrayfun(@(i) sprintf('Satellite %d', i), 1:Nsat, 'UniformOutput', false)];
legend(labels, 'Location','best');

grid on;
title('Ca^{2+} evolution in astrocyte ER network');

end
