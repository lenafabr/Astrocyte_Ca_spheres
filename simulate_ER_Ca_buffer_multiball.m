function sim = simulate_ER_Ca_buffer_multiball(params)
% simulate_ER_Ca_multiball
%
% Simulate Ca2+ concentration evolving in an ER network made of:
%   - one central ER region ("cell body" ER), composed of a 3D tubular network
%   - several small ER "balls" connected to the cell body by tubes. These
%   "balls" are assumed to be made of solid ER lumen.
%   - The satellite balls are connected to the cell body via multiple
%   parallel tubules of ER.

% below parameters should be passed in a params structure
%
% Geometry:
%
% R_body = radius of cell body
% Vfrac_body = volume fraction of ER within the cell body
% r_er = radius of ER tubules
% R_sat = column vector of satellite ball radii
% L_tube = column vector, lengths of tubes connecting the satellites to the cell body
% N_tube = column vector, number of parallel tubes in each neck connecting a particular satellite to cell body
% The cross-sectional ER area linking a satellite to the cell body is A_i = N_tube(i)*pi*r_er^2
%
% Dynamic and binding parameters
%
% Du = diffusivity of free Ca
% Db = diffusivity of bound Ca
% S = total buffer site concentration in entire system (buffers are assumed
% to be uniformly distributed throughout)
% K_D = dissociation constant for binding to each buffer site
% 
% tvals = row vector of time values at which concentrations will be
% evaluated
%
% U_init = column vector; initial free Ca concentrations in cell body followed by all
% satellites
% -----------------
%
% Output:
%   sim : struct with simulation results:
%       .tvals, .U_vals, .Nsat       
%       .params, .geom

if nargin < 1
    params = struct;
end

%% 1. Fill in defaults
Du          = getf(params, 'D',          30);          % um^2/s, diffusivity of free Ca
Db          = getf(params, 'D',          3);          % um^2/s, diffusivity of bound Ca
r_er       = getf(params, 'r_er',       0.050);       % um, radius of ER tube (can have multiple tubes in parallel in each projection)
Vfrac_body = getf(params, 'Vfrac_body', 0.0015);     % volume fraction of ER in cell body
R_body     = getf(params, 'R_body',     7);           % um, radius of cell body sphere

% satellites. List parameters for each satellite ball
% default is 2 satellites
R_sat  = getf(params, 'R_sat',  [4; 3]);    % radius of each satellite ball (um)
L_tube = getf(params, 'L_tube', [20; 40]); % length of neck connecting each satellite to cell body (um)
N_tube = getf(params, 'N_tube', [4; 4]); % within the long neck, how many ER tubules are in parallel
% make sure they are column vectors
if (size(R_sat,2)>size(R_sat,1))
    R_sat = R_sat';
end
if (size(L_tube,2)>size(L_tube,1))
    L_tube = L_tube';
end
if (size(N_tube,2)>size(N_tube,1))
    N_tube = N_tube';
end

% check inputs make sense
Nsat = numel(R_sat); % number of satellites
if numel(L_tube) ~= Nsat || numel(N_tube) ~= Nsat
    error('R_sat, L_tube, and N_tube must have the same length.');
end

% Bufer protein parameters
S = getf(params, 'S', 2.7); % concentration of buffer sites, mM
Kd = getf(params, 'Kd', 0.2); % buffer binding dissociation constant, mM

% time values at which to evaluate concentrations
tvals = getf(params, 'tvals', [0,2000]);
% let's make this a row vector
if (size(tvals,1)>size(tvals,2))
    tvals = tvals';
end
Nt = numel(tvals);

% initial concentrations of free Ca (in mM)
% default to everything full of Ca
U_init = getf(params, 'U_init', 0.5*ones(Nsat+1,1));
if (size(U_init,1)<size(U_init,2))
    U_init = U_init';
end
if numel(U_init) ~= Nsat+1
    error('U_init must have length Nsat+1.');
end

% volumes
V_ball = zeros(Nsat+1,1);
V_ball(1) = (4/3)*pi*R_body^3 * Vfrac_body;  % ER volume in cell body
V_ball(2:Nsat+1)  = (4/3)*pi * R_sat.^3;             % satellite ER volumes

A_tube = N_tube .* pi * r_er.^2;          % cross-sectional area in the neck connected to each satellite


% Effective diffusive resistance such that at steady state, the current between two balls
% is given by 1/Reff* the difference in concentration between the balls
Reff   = 1./(N_tube .* r_er) + L_tube./(A_tube);



%% Run the simulation
%     options = odeset('MaxStep', 1e-3);  % e.g., max time step = 0.001
% U_sol(i,j) is the free Ca conc at time t_sol(i), in ball j (1st ball is
% cell body)
[t_sol, U_sol] = ode45(@(t, U) compute_dUdt(U,Reff,V_ball,Du,Db,S,Kd), tvals,U_init);    


%% 5. Pack output

params_out            = params;
params_out.Du          = Du;
params_out.Db          = Db;
params_out.r_er       = r_er;
params_out.Vfrac_body = Vfrac_body;
params_out.R_body     = R_body;
params_out.R_sat      = R_sat;
params_out.L_tube     = L_tube;
params_out.N_tube     = N_tube;
params_out.tvals          = tvals;
params_out.U_init = U_init;
params_out.S = S;
params_out.Kd = Kd;

geom = struct;
geom.V_ball = V_ball;
geom.A_tube = A_tube;
geom.Reff   = Reff;
geom.N_tube = N_tube;
geom.Nsat = Nsat;

sim = struct;
sim.t_sol = t_sol;
sim.U_sol = U_sol;
% total Ca in each ball, over time
sim.C_sol = U_sol.*(1 + S./(U_sol+Kd));
sim.params  = params_out;
sim.geom    = geom;

end

% helper: get field with default
function val = getf(s, fname, default)
% Return s.(fname) if it exists and is non-empty; otherwise return default.
    if isfield(s, fname) && ~isempty(s.(fname))
        val = s.(fname);
    else
        val = default;
    end
end
