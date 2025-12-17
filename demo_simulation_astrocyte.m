%
% Demo script for Ca2+ diffusing in ER in astrocyte cell. The astrocyte has
% one central cell body ER region, and several ER satellite regions
% connected to the cell body by tubes. 
% In the central body, the ER is a (unresolved) network that occupies some
% fraction of the cell body volume. In the satellites, each ER blob is taken to
% be a sphere that represents the full volume.
%
% We assume the existence of spatially
% uniform Ca2+ buffer proteins in ER lumen, and the binding and unbinding
% of free Ca2+ and the buffer protein is very fast comparing to diffusion
% so the binding reaction is always at equilibrium. The diffusion of Ca2+
% is a combination of diffusion of free Ca2+ and bound Ca2+.
%
%
% This script:
%   1) Simulates Ca2+ diffusion between the body and satellites over time
%      using a linear ODE model with effective diffusive resistances.
%      (simulate_ER_Ca_multiball.m) 
%   2) Plots Ca2+ concentration vs time in each compartment.
%      (plot_ER_Ca_timecourses.m) 
%   3) Plots 2D snapshots of the ER structure colored by free Ca2+
%       depletion fraction compared to a reference Ca2+ concentration. (plot_ER_Ca_snapshots.m) 
%
% -------------------------
% Physical / geometric parameters:
%   Du         : diffusivity of free Ca in ER lumen (um^2/s).
%
%   Db         : diffusivity of bound Ca in ER lumen (um^2/s).
% 
%   S          : concentration of buffer sites for Ca (mM), assumed to be uniform in space
%
%   Kd         : Ca buffer binding dissociation constant (mM), Kd = U*(S-B)/B
%
%   r_er       : radius of a single ER tube (um).
%                Each connection has N_tube(i) such tubes in parallel.
%
%   Vfrac_body : fraction of the cell-body volume occupied by ER (0 to 1).
%                We assume the cell body is a sphere of radius R_body,
%                and only Vfrac_body of that is the ER network where Ca2+ diffuses.
%
%   R_body     : radius of the spherical cell body (um).
%
%   R_sat(i)   : radius of satellite spherical ER region i (um).
%
%   L_tube(i)  : length of the ER tubes connecting the body to satellite i (um).
%
%   N_tube(i)  : number of parallel tubes connecting body to satellite i.               
%
% Time and initial conditions:
%   tvals      : start and stop time (time values to be decided automatically) OR a list of desired time values (sec)
%   U_init(1)  : initial free Ca2+ concentration in the cell-body ER (mM).
%   U_init(i+1): initial free Ca2+ concentration in satellite i (mM).
%
% OUTPUTS (from simulate_ER_Ca_multiball):
%   sim.t_sol    : time vector (s).
%   sim.U_sol    : free Ca2+ in cell body and satellites over time (mM).
%   sim.C_sol    : total Ca2+ in cell body and satellites over time (mM).
%   sim.params   : parameter struct (with defaults filled in if not supplied by user).
%   sim.geom     : geometry-related volumes and resistances.

%% Define parameters
% leave empty to use defaults (defined within simulate_ER_Ca_buffer_multiball)
params = struct();

% Physical parameters
params.Du         = 30;          % um^2/s, diffusivity of free Ca in ER lumen
params.Db         = 3;           % um^2/s, diffusivity of bound Ca in ER lumen
params.r_er       = 0.050;       % um, radius of a single ER tube (~50 nm), (can have multiple tubes in parallel in each projection)
params.Vfrac_body = 0.0015;      % ER volume fraction in cell body, estimated using the 3D hexagonal ER network
params.R_body     = 7;           % um, cell body radius 
params.S          = 2.7;         % mM, concentration of buffer sites for Ca, assumed to be uniform in space
params.Kd         = 0.2;         % mM, Ca buffer binding dissociation constant, Kd = U*(S-B)/B

% Geometry: three satellite ER "balls" attached to the body
params.R_sat  = [4; 3];          % um, radii of satellites
params.L_tube = [20; 40];        % um, tube lengths body -> satellite
params.N_tube = [4; 4];          % number of parallel ER tubules in each body -> sattellite connection 

% Time and initial conditions
params.tvals     = [0,2000];      % s, start and stop time. The time steps will be decided by ode45
Nsat = numel(params.R_sat); % number of sattelites 
params.U_init    = 0.5*ones(1+Nsat,1); % initial free Ca concentration, 0.5 mM everywhere. [cell body; sat_1; sat_2;...]
params.U_init(2) = 0.5*0.2; % deplete one of the satelites by 80%

%% ----------------------- 2. Run simulation ---------------------------
sim = simulate_ER_Ca_buffer_multiball(params);

%% ----------------------- 3. Plot time courses -----------------------
plot_ER_Ca_timecourses(sim)

%% ----------------------- 4. Snapshot plots ---------------------------

snapshotTimes = [0, 10, 50, 500];   % times at which to plot structure (s)
% plot color as free Ca depletion relative to a particular starting concentration
U_rel = sim.U_sol(1,1); % use the initial cell body concentration as the reference concentration 
plot_ER_Ca_snapshots(sim, snapshotTimes,U_rel);