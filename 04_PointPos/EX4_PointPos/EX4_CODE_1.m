%% EX.4 -- POINT POSITIONING

%% STEP 1 - Parameters & Variables
clc
clear

%Variables
run('Ex04_variables.m')

%Parameters
c=299792458; %speed of light [m/s]
itermax = 10; %max # of iterations of the loop
ths = 0.1; %threshold condition for loop
x_in = [0; 0; 0; 0]; %initial approximated values of rec. & clock offset

%% STEP 2 - Iteration & Threshold

%initialization of empty vector
x_rec = zeros(4,11);
disp('Quella giusta');
%Iteration on LS
for i=1:itermax
    
    x_approx = x_in(1:3);
    
    % the function topocent computes azimuth, elevation and distance, given
    % the sat/rec coordinates
    [azimuth, elevation, ro] = topocent(x_approx, xyz_sat);
    
    % convert into geodetic coordinates
    [phi, lambda, h] = cart2geod(x_approx(1), x_approx(2), x_approx(3));
    
    % convert from radians to degrees
    phi_rec = phi;%*180/pi; 
    lam_rec = lambda;%*180/pi;
    
    % the functions tropo_correction and iono_corrections compute the iono/
    % tropospheric effect according to the Klobuchar model and the
    % Saastamoinen model respectively
    tropo_eff = tropo_correction(h, elevation);
    iono_eff = iono_correction(phi_rec, lam_rec, azimuth, elevation, time_rx, ionoparams);
    
    % Matrix CSI of known terms and reduced term (for Least Square calc.)
    csi = ro + tropo_eff + iono_eff - c*dtS;
    deltaP = pr_C1 - csi;

    % Matrix A for Least Square calc. 
    A_1 = [(x_approx(1)-xyz_sat(:,1))./ro (x_approx(2)-xyz_sat(:,2))./ro (x_approx(3)-xyz_sat(:,3))./ro];
    % this is a matrix composed of the position of the rec. - the position
    % of the satellite, divided by the distance between the two.
    % In the first loop we use the approx value.
    A = [A_1 ones(11,1)];
    % In this second line we add the column of known constant term 1 
    
    % First estimate of x after LS
    N = A' * A;
    xcap = inv(N) * A' * deltaP;

    % Newly estimated coordinates of the receiver (that used to be 0), approx + est.: 
    x_rec(:,i) = x_in + xcap;
    
    % Delta for threshold comparison
    diff = x_rec(:,i) - x_in; %technically, this means new values - old ones
    
    % Subtitution of new estimates in initial estimate to redo the loop
    x_in = x_rec(:,i);
    
    % Saving the coordinates in case they're good
    coords_final = x_rec(1:3,i);
    
    % Convergence check. If one of the values is lower than threshold, stop the loopcheck convergence of the result and, in case exit
        
    if max(abs(diff(1:3))) < ths
        break
    end
    
    % Need to also check that we didn't overdo the # of iterations
    if i == itermax
       disp("Went over max number of iterations");
    end

    % End of this amazing LS loop
    
end

%% STEP 3 - PDOP

% Now we need the cofactor matrix, which is = N, thus:
Qxx = inv(A'*A);
% Then we get the covariance matrix of the estimates:
Qx = Qxx(1:3,1:3);

% Rotation Matrix
[R] = RotMatrix(phi, lambda);

% Rotation by covariance propagation
Q = R*Qx*R';
% This matrix contains on the diagonal the elements qee, qnn, quu, which
% are needed in the following step

% PDOP calculation
PDOP = sqrt(Q(1,1)+ Q(2,2) + Q(3,3)); % qee + qnn + quu

dtr = (x_rec(4,i)-x_rec(4,i-1))/c;

% MAIN:
disp('Output of non-COANGLE');
disp('-------------------');
iterp = sprintf('Total # of Iter: %d', i);
coordp = sprintf('Receiver Coords: %d; %d; %d;', coords_final);
offsetp = sprintf('Receiver Clock Offset: %d', dtr);
PDOPp = sprintf('PDOP: %f', PDOP);
disp(iterp);
disp(coordp);
disp(offsetp);
disp(PDOPp);

%% STEP 4 - SAME BUT WITH A CUT-OFF ANGLE - PRE-PROCESSING OF DATA

% We now need to calculate how many sats are above the cut-off angle
% The cut-off angle here is reported as the elevation at which we start
% seeing sats, therefore if a sat has elevation < than cut-off angle, we
% don't keep it in consideration
coangle = 5; % cut-off angle
nsat = length(dtS); % # of sats
empty = zeros(11,1); % vector used for storage purposes.

%Loop to check elevation and store sats
for i=1:nsat
    if elevation(i) > coangle
        empty(i)=elevation(i);
    end
end

%We use function find in order to find the non-zero elements of vector
%empty
[index] = find(empty);

% The new number of sats that pass the cut-off angle test is:
newsat = length(index);

%% STEP 5 - LS & ITERATION WITH CUT-OFF ANGLE

%New initial vectors (we use different ones because otherwise it messes up
%the code and I am not good enough to find a different way to go around
%this problem. Sorry about this. 

n_x_in = [coords_final; dtr];
n_x_rec = zeros(4,10);

for n_i=1:itermax
    new_x_approx = n_x_in(1:3);
    [n_azimuth, n_elevation, n_distance] = topocent(new_x_approx, xyz_sat(index,:));
    [n_phi, n_lambda, n_h] = cart2geod(new_x_approx(1), new_x_approx(2), new_x_approx(3));
    n_phi = n_phi*180/pi;
    n_lambda = n_lambda*180/pi;
    n_err_tropo = tropo_correction(n_h, n_elevation);
    n_err_iono = iono_correction(n_phi, n_lambda, n_azimuth, n_elevation, time_rx, ionoparams);
    n_csi = n_distance +n_err_tropo+n_err_iono -c*dtS(index);
    n_dP = pr_C1(index) - n_csi;
    n_A_1 = [(new_x_approx(1)-xyz_sat(index,1))./n_distance (new_x_approx(2)-xyz_sat(index,2))./n_distance (new_x_approx(3)-xyz_sat(index,3))./n_distance];
    n_A = [n_A_1 ones(newsat,1)];
    n_xcap=n_A\n_dP;
    n_x_rec(:, n_i) = n_x_in + n_xcap;
    n_delta = n_x_rec(:,n_i) - n_x_in;
    n_x_in = n_x_rec(:,n_i);
    n_coords_final = n_x_rec(1:3,n_i);
    
    % Threshold check, if lower, break
    if  max(abs(n_delta(1:3))) < ths
        break
    end

      % Need to also check that we didn't overdo the # of iterations
    if n_i == itermax
        print("Fool of a Took");
    end
end

%% STEP 6 - PDOP WITH COANGLE

% Now we need the cofactor matrix, which is = N, thus:
n_Qxx = inv(n_A'*n_A);
% Then we get the covariance matrix of the estimates:
n_Qx = n_Qxx(1:3,1:3);

% Rotation Matrix
[n_R] = RotMatrix(n_phi, n_lambda);

% Rotation by covariance propagation
n_Q = n_R*n_Qx*n_R';
% This matrix contains on the diagonal the elements qee, qnn, quu, which
% are needed in the following step

% PDOP calculation
n_PDOP = sqrt(n_Q(1,1)+ n_Q(2,2) + n_Q(3,3)); % qee + qnn + quu


%% STEP 7 - OUTPUTS (COANGLE AND NOT COANGLE)

%Expected output:
%1. Coordinates of the receiver and its clock offset.
%2. Variance covariance matrix of the coordinates.
%3. PDOP

%COANGLE:
disp ('-------------');
disp('Output of COANGLE');

coanglep = sprintf('CUT-OFF ANGLE: %d',coangle);
n_ip = sprintf('Total # of Iter: %d', n_i);
n_coordp = sprintf('Receiver Coords: [%d %d %d]', n_coords_final);
n_offsetp = sprintf('Receiver Clock Offset: %d', (n_x_rec(4,n_i)-n_x_rec(4,n_i-1))/c);
n_PDOPp = sprintf('PDOP: %f', n_PDOP);

disp(coanglep);
disp(n_ip);
disp(n_coordp);
disp(n_offsetp);
disp(n_PDOPp);



