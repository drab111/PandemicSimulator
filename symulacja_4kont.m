function [t_out, S, I, R, D] = symulacja_4kont(userParam)
% SYMULACJA_4KONT - SIRD model with migration for 4 continents
%
% [t_out, S, I, R, D] = symulacja_4kont(userParam)
%
% PARAMETERS (in the structure userParam):
%   .startCont  = 1,2,3,4 (continent where the pandemic starts)
%   .population = [pop1, pop2, pop3, pop4] (population numbers)
%   .beta       = [b1, b2, b3, b4]
%   .gamma      = [g1, g2, g3, g4]
%   .mu         = [m1, m2, m3, m4]
%   .alpha      = 4x4 matrix, alpha(i,j) = migration from i to j
%   .tEnd       = simulation end time
%   .birth      = [br1, br2, br3, br4] (birth rates)
%   .otherDeath = [d1, d2, d3, d4] (mortality from other causes)
%
% RETURNS:
%   t_out   - time vector
%   S, I, R, D - matrices [length(t_out) x 4], 
%               columns are continents: 
%               1 - Eurasia, 2 - Africa, 3 - America, 4 - Oceania

    % Set default parameters
    % Approximate populations:
    % 1 - Eurasia, 2 - Africa, 3 - America, 4 - Oceania
    defaultParam.startCont  = 1;  % Default pandemic starts in Eurasia
    defaultParam.population = [5.3e9, 1.4e9, 1.0e9, 0.05e9]; 
    defaultParam.beta  = [0.25, 0.30, 0.20, 0.15];  
    defaultParam.gamma = [0.05, 0.06, 0.04, 0.04];  
    defaultParam.mu    = [0.01, 0.015, 0.008, 0.005]; 

    % Migration matrix alpha (4x4). alpha(i,j) = flow from i to j.
    defaultParam.alpha = [ ...
        0,     0.5e-5,  3e-5, 0.3e-5;
        3.5e-5,  0,     1e-5, 0.1e-5;
        2.5e-5,  0.7e-5,  0,    1e-5;
        1e-5,  0.5e-5,  2e-5, 0     ];
    
    defaultParam.tEnd = 365;  % 1 year simulation

    % Defaults to zero (if user does not provide)
    defaultParam.birth      = [0, 0, 0, 0];  % default birth rate
    defaultParam.otherDeath = [0, 0, 0, 0];  % default overall mortality

    if nargin < 1
        param = defaultParam;
    else
        % Override values provided by user
        param = defaultParam;
        fieldsUser = fieldnames(userParam);
        for k = 1:length(fieldsUser)
            param.(fieldsUser{k}) = userParam.(fieldsUser{k});
        end
    end

    % Choose continent startCont -> set some initial infected there
    I0_all = [0, 0, 0, 0];
    
    % Check if userParam contains 'userInf'
    if isfield(userParam,'userInf')
        % Use the value given by user
        I0_all(param.startCont) = userParam.userInf;
    else
        % Default 10 thousand
        I0_all(param.startCont) = 1e4;
    end
    % -----------------------------------------------------------------

    % Assemble vector Y0 = [S1; I1; R1; D1; S2; I2; R2; D2; ... S4; I4; R4; D4]
    Y0 = zeros(16,1); 
    idx = 1;
    for i = 1:4
       S_i0 = param.population(i) - I0_all(i);
       I_i0 = I0_all(i);
       R_i0 = 0;
       D_i0 = 0;
       Y0(idx)   = S_i0;   % S_i
       Y0(idx+1) = I_i0;   % I_i
       Y0(idx+2) = R_i0;   % R_i
       Y0(idx+3) = D_i0;   % D_i
       idx = idx + 4;
    end

    % Define time interval
    tSpan = [0, param.tEnd];

    % Call ODE solver
    [tSol, YSol] = ode45(@(t,Y) pandemic_odes4(t, Y, param), tSpan, Y0);

    % Unpack results: 
    %   YSol(:,1) = S1,  (:,2) = I1,  (:,3) = R1,  (:,4) = D1
    %   YSol(:,5) = S2,  (:,6) = I2,  ...
    %   ...
    %   YSol(:,13)= S4, (:,14)= I4,  (:,15)= R4, (:,16)= D4
    S = zeros(length(tSol),4);
    I = zeros(length(tSol),4);
    R = zeros(length(tSol),4);
    D = zeros(length(tSol),4);

    for i = 1:4
        colBase = 4*(i-1);
        S(:,i) = YSol(:, colBase+1);
        I(:,i) = YSol(:, colBase+2);
        R(:,i) = YSol(:, colBase+3);
        D(:,i) = YSol(:, colBase+4);
    end

    t_out = tSol;

    % Plot results: 4 subplots (R, S+I+R, I, D) ===
    figure('Name','Pandemic Simulation - 4 continents','Color','w');
    
    % (A) Recovered (R) plot
    subplot(2,2,1);
    plot(t_out, R(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, R(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, R(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, R(:,4), 'm-', 'LineWidth', 2);
    xlabel('Time [days]');
    ylabel('Recovered R_i(t)');
    legend({'Eurasia','Africa','America','Oceania'}, 'Location','best');
    grid on; 
    title('Recovered (R)');

    % (B) Living population (S + I + R) plot
    subplot(2,2,2);
    plot(t_out, S(:,1)+I(:,1)+R(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, S(:,2)+I(:,2)+R(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, S(:,3)+I(:,3)+R(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, S(:,4)+I(:,4)+R(:,4), 'm-', 'LineWidth', 2);
    xlabel('Time [days]');
    ylabel('Living population');
    legend({'Eurasia','Africa','America','Oceania'}, 'Location','best');
    grid on; 
    title('Number of living (S+I+R)');

    % (C) Infected (I) plot
    subplot(2,2,3);
    plot(t_out, I(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, I(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, I(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, I(:,4), 'm-', 'LineWidth', 2);
    xlabel('Time [days]');
    ylabel('Infected I_i(t)');
    legend({'Eurasia','Africa','America','Oceania'}, 'Location','best');
    grid on; 
    title('Infected (I)');

    % (D) Deaths (D) plot
    subplot(2,2,4);
    plot(t_out, D(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, D(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, D(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, D(:,4), 'm-', 'LineWidth', 2);
    xlabel('Time [days]');
    ylabel('Deaths D_i(t)');
    legend({'Eurasia','Africa','America','Oceania'}, 'Location','best');
    grid on; 
    title('Deaths (D)');

end


% HELPER FUNCTION: ODE equations for 4 continents
function dYdt = pandemic_odes4(t, Y, param)
    % Y = [S1; I1; R1; D1; S2; I2; R2; D2; S3; I3; R3; D3; S4; I4; R4; D4]
    % Unpack:
    S1 = Y(1);  I1 = Y(2);  R1 = Y(3);  D1 = Y(4);
    S2 = Y(5);  I2 = Y(6);  R2 = Y(7);  D2 = Y(8);
    S3 = Y(9);  I3 = Y(10); R3 = Y(11); D3 = Y(12);
    S4 = Y(13); I4 = Y(14); R4 = Y(15); D4 = Y(16);

    % From parameters:
    beta1  = param.beta(1);   beta2  = param.beta(2);
    beta3  = param.beta(3);   beta4  = param.beta(4);
    gamma1 = param.gamma(1);  gamma2 = param.gamma(2);
    gamma3 = param.gamma(3);  gamma4 = param.gamma(4);
    mu1    = param.mu(1);     mu2    = param.mu(2);
    mu3    = param.mu(3);     mu4    = param.mu(4);

    birth1  = param.birth(1);  birth2  = param.birth(2);
    birth3  = param.birth(3);  birth4  = param.birth(4);

    death1  = param.otherDeath(1);  death2  = param.otherDeath(2);
    death3  = param.otherDeath(3);  death4  = param.otherDeath(4);

    alpha = param.alpha; % migration matrix

    % Total population per continent (excluding dead)
    N1 = S1 + I1 + R1;
    N2 = S2 + I2 + R2;
    N3 = S3 + I3 + R3;
    N4 = S4 + I4 + R4;

    % Infection rates (new infections per susceptible)
    inf1 = beta1 * S1 * I1 / N1;
    inf2 = beta2 * S2 * I2 / N2;
    inf3 = beta3 * S3 * I3 / N3;
    inf4 = beta4 * S4 * I4 / N4;

    % Migration flows (susceptible):
    % For each continent i, net migration into i is sum_j(alpha(j,i)*S_j) - sum_j(alpha(i,j)*S_i)
    % Similarly for I and R.
    dS = zeros(4,1);
    dI = zeros(4,1);
    dR = zeros(4,1);
    dD = zeros(4,1);

    S_arr = [S1; S2; S3; S4];
    I_arr = [I1; I2; I3; I4];
    R_arr = [R1; R2; R3; R4];

    for i = 1:4
        % Migration for S:
        mig_S_in = 0;
        mig_S_out = 0;
        mig_I_in = 0;
        mig_I_out = 0;
        mig_R_in = 0;
        mig_R_out = 0;
        for j = 1:4
            if i ~= j
                mig_S_in = mig_S_in + alpha(j,i)*S_arr(j);
                mig_S_out = mig_S_out + alpha(i,j)*S_arr(i);
                mig_I_in = mig_I_in + alpha(j,i)*I_arr(j);
                mig_I_out = mig_I_out + alpha(i,j)*I_arr(i);
                mig_R_in = mig_R_in + alpha(j,i)*R_arr(j);
                mig_R_out = mig_R_out + alpha(i,j)*R_arr(i);
            end
        end

        switch i
            case 1
                % Susceptible change:
                dS(i) = birth1*N1 - inf1 - death1*S1 + mig_S_in - mig_S_out;
                % Infected change:
                dI(i) = inf1 - (gamma1 + mu1 + death1)*I1 + mig_I_in - mig_I_out;
                % Recovered change:
                dR(i) = gamma1*I1 - death1*R1 + mig_R_in - mig_R_out;
                % Death change:
                dD(i) = mu1*I1;
            case 2
                dS(i) = birth2*N2 - inf2 - death2*S2 + mig_S_in - mig_S_out;
                dI(i) = inf2 - (gamma2 + mu2 + death2)*I2 + mig_I_in - mig_I_out;
                dR(i) = gamma2*I2 - death2*R2 + mig_R_in - mig_R_out;
                dD(i) = mu2*I2;
            case 3
                dS(i) = birth3*N3 - inf3 - death3*S3 + mig_S_in - mig_S_out;
                dI(i) = inf3 - (gamma3 + mu3 + death3)*I3 + mig_I_in - mig_I_out;
                dR(i) = gamma3*I3 - death3*R3 + mig_R_in - mig_R_out;
                dD(i) = mu3*I3;
            case 4
                dS(i) = birth4*N4 - inf4 - death4*S4 + mig_S_in - mig_S_out;
                dI(i) = inf4 - (gamma4 + mu4 + death4)*I4 + mig_I_in - mig_I_out;
                dR(i) = gamma4*I4 - death4*R4 + mig_R_in - mig_R_out;
                dD(i) = mu4*I4;
        end
    end

    % Pack derivatives into dYdt vector
    dYdt = zeros(16,1);
    for i = 1:4
        base = 4*(i-1);
        dYdt(base+1) = dS(i);
        dYdt(base+2) = dI(i);
        dYdt(base+3) = dR(i);
        dYdt(base+4) = dD(i);
    end
end
