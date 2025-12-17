function dUdt = compute_dUdt(Uvals,Reff,V_ball,Du,Db,S,Kd)
% compute the rate of change of free calcium concentration in the cell
% body and all satellites
% inputs:
% U = column vector of free calcium (cell body then each satellite)
% Reff = column vector for the diffusive resistivity of the neck connecting
% each satellite to the cell body. Reff should not include the D
% coefficient
% Du, Db = diffusivity of free and bound Ca
% S = total buffer binding site concentration (uniform throughout)
% Kd = buffer dissociation constant

% bound Ca concentration
Bvals = Uvals*S./(Uvals+Kd);

% concentration differences
Udiff = Uvals(2:end)-Uvals(1);
Bdiff = Bvals(2:end)-Bvals(1);

% current through each neck
currents = (Du.*Udiff + Db.*Bdiff)./Reff;

% rate of change in total Ca concentration
dCdt = zeros(size(Uvals));
dCdt(1) = sum(currents)./V_ball(1);
dCdt(2:end) = -currents./V_ball(2:end);

% rate of change in free Ca concentration
dUdt = dCdt./(1 + S*Kd./(Uvals + Kd).^2);

end