function geom2D = compute_ER_radial_geometry(params)
% compute_ER_radial_geometry
%
% Compute 2D radial layout for plotting:
%   - cell body centered at (0,0)
%   - satellites arranged radially with equal angle gaps
%   - tubes as rectangles with some visual thickness

R_body  = params.R_body;
R_sat   = params.R_sat(:);
L_tube  = params.L_tube(:);
N_tube  = params.N_tube(:);
r_er    = params.r_er;

Nsat = numel(R_sat);

% 2D circle for cell body
ntheta = 200;
thetaC = linspace(0, 2*pi, ntheta);
x_body = R_body * cos(thetaC);
y_body = R_body * sin(thetaC);

% radial angles for satellites
phi = linspace(0, 2*pi, Nsat+1);
phi = phi(1:Nsat);     % drop duplicate

% visual tube thickness
r_tube = N_tube .* r_er;     % physical
tubeScale = 10;              % purely visual
r_vis = tubeScale * r_tube;  % y-thickness

satGeom = struct;

for i = 1:Nsat
    dir  = [cos(phi(i)); sin(phi(i))];
    perp = [-sin(phi(i)); cos(phi(i))];
    
    % positions along radial direction
    P_body     = R_body * dir;                       % body surface
    P_tubeEnd  = (R_body + L_tube(i)) * dir;         % end of tube
    P_satCtr   = (R_body + L_tube(i) + R_sat(i)) * dir;  % satellite center
    
    % satellite circle
    x_sat_i = P_satCtr(1) + R_sat(i)*cos(thetaC);
    y_sat_i = P_satCtr(2) + R_sat(i)*sin(thetaC);
    
    % tube rectangle
    P_start = P_body;
    P_end   = P_tubeEnd;
    
    P1 = P_start + r_vis(i)*perp;
    P2 = P_end   + r_vis(i)*perp;
    P3 = P_end   - r_vis(i)*perp;
    P4 = P_start - r_vis(i)*perp;
    
    x_rect = [P1(1), P2(1), P3(1), P4(1)];
    y_rect = [P1(2), P2(2), P3(2), P4(2)];
    
    satGeom(i).x_sat  = x_sat_i;
    satGeom(i).y_sat  = y_sat_i;
    satGeom(i).x_tube = x_rect;
    satGeom(i).y_tube = y_rect;
end

Rmax = R_body + max(L_tube + R_sat);

geom2D = struct;
geom2D.x_body = x_body;
geom2D.y_body = y_body;
geom2D.satGeom = satGeom;
geom2D.Rmax = Rmax;

end
