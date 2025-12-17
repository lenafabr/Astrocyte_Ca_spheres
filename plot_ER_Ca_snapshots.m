function plot_ER_Ca_snapshots(sim, snapshotTimes, U_rel)
% plot_ER_Ca_snapshots
%
% Make snapshot plots of ER structure (cell body + satellites)
% colored by free Ca2+ relative depletion at several time points.
%
% Inputs:
%   sim           : struct from simulate_ER_Ca_multiball
%   snapshotTimes : vector of times (in same units as sim.t)
%                   If empty or not provided, choose 4 default times.
% C_rel: total calcium concentration such that the fractional depletion
% relative to this will be shown as the color. 

if nargin < 2
    snapshotTimes = [];
end

params = sim.params;
t      = sim.t_sol;
% free calcium concentrations, relative to initial
U_body = (U_rel - sim.U_sol(:,1))./U_rel;
U_sat  = (U_rel - sim.U_sol(:,2:end))./U_rel;


Nsat   = sim.geom.Nsat;

geom2D = compute_ER_radial_geometry(params);
x_body = geom2D.x_body;
y_body = geom2D.y_body;
satGeom = geom2D.satGeom;
Rmax   = geom2D.Rmax;

% choose default snapshots: 4 times from start to end
if isempty(snapshotTimes)
    snapshotTimes = [t(1), ...
                     t(1) + 0.25*(t(end)-t(1)), ...
                     t(1) + 0.50*(t(end)-t(1)), ...
                     t(end)];
end

% indices of nearest time points
idxSnap = zeros(size(snapshotTimes));
for j = 1:numel(snapshotTimes)
    [~, idxSnap(j)] = min(abs(t - snapshotTimes(j)));
end

% global color limits across all compartments and times
cAll = [U_body, U_sat];
cmin = min(cAll(:));
cmax = max(cAll(:));

figure;
tiledlayout(1,numel(idxSnap), 'TileSpacing','compact','Padding','compact');

for j = 1:numel(idxSnap)
    k  = idxSnap(j);
    tk = t(k);

    cbk = U_body(k);
    csk = U_sat(k,:);

    nexttile(j); hold on;

    % cell body
    patch(x_body, y_body, cbk*ones(size(x_body)), 'EdgeColor','none');

    % satellites + gradient tubes
    for i = 1:Nsat
        % satellite circle
        patch(satGeom(i).x_sat, satGeom(i).y_sat, ...
              csk(i)*ones(size(satGeom(i).x_sat)), ...
              'EdgeColor','none');

        % gradient tube: body -> satellite
        c_tube = [cbk, csk(i), csk(i), cbk];
        patch(satGeom(i).x_tube, satGeom(i).y_tube, c_tube, ...
              'EdgeColor','none', 'FaceColor','interp');
    end

    axis equal;
    xlim([-Rmax*1.3, Rmax*1.3]);
    ylim([-Rmax*1.3, Rmax*1.3]);
    colormap(parula);
    caxis([cmin cmax]);

    title(sprintf('t = %.2f s', tk));
    xlabel('x (\mum)');
    ylabel('y (\mum)');
end

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb, 'fraction free Ca depleted');

end
