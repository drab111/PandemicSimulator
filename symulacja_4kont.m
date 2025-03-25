function [t_out, S, I, R, D] = symulacja_4kont(userParam)
% SYMULACJA_4KONT - model SIRD z migracją dla 4 kontynentów
%
% [t_out, S, I, R, D] = symulacja_4kont(userParam)
%
% PARAMETRY (w strukturze userParam):
%   .startCont  = 1,2,3,4 (kontynent, gdzie zaczyna się pandemia)
%   .population = [pop1, pop2, pop3, pop4] (liczba mieszkańców)
%   .beta       = [b1, b2, b3, b4]
%   .gamma      = [g1, g2, g3, g4]
%   .mu         = [m1, m2, m3, m4]
%   .alpha      = macierz 4x4, alpha(i,j) = migracja z i do j
%   .tEnd       = koniec czasu symulacji
%   .birth      = [br1, br2, br3, br4] (szybkość narodzin)
%   .otherDeath = [d1, d2, d3, d4] (śmiertelność z innych przyczyn)
%
% ZWRACA:
%   t_out   - wektor czasu
%   S, I, R, D - macierze [length(t_out) x 4], 
%               w kolumnach kontynenty: 
%               1 - Eurazja, 2 - Afryka, 3 - Ameryka, 4 - Oceania


    % Ustawiamy parametry domyślne
    % Przybliżone populacje:
    % 1 - Eurazja, 2 - Afryka, 3 - Ameryka, 4 - Oceania
    defaultParam.startCont  = 1;  % Domyślnie pandemia startuje w Eurazji
    defaultParam.population = [5.3e9, 1.4e9, 1.0e9, 0.05e9]; 
    defaultParam.beta  = [0.25, 0.30, 0.20, 0.15];  
    defaultParam.gamma = [0.05, 0.06, 0.04, 0.04];  
    defaultParam.mu    = [0.01, 0.015, 0.008, 0.005]; 

    % Macierz migracji alpha (4x4). alpha(i,j) = przepływ z i do j.
    defaultParam.alpha = [ ...
        0,     0.5e-5,  3e-5, 0.3e-5;
        3.5e-5,  0,     1e-5, 0.1e-5;
        2.5e-5,  0.7e-5,  0,    1e-5;
        1e-5,  0.5e-5,  2e-5, 0     ];
    
    defaultParam.tEnd = 365;  % 1 rok symulacji

    % Domyślnie zero (jeśli user nie poda)
    defaultParam.birth      = [0, 0, 0, 0];  % domyślna szybkość narodzin
    defaultParam.otherDeath = [0, 0, 0, 0];  % domyślna śmiertelność ogólna

    if nargin < 1
        param = defaultParam;
    else
        % Nadpisujemy wartości, które użytkownik podał
        param = defaultParam;
        fieldsUser = fieldnames(userParam);
        for k = 1:length(fieldsUser)
            param.(fieldsUser{k}) = userParam.(fieldsUser{k});
        end
    end

    % Wybór kontynentu startCont -> tam dajemy pewną liczbę zarażonych na starcie
    I0_all = [0, 0, 0, 0];
    
    % Sprawdzamy czy userParam zawiera 'userInf'
    if isfield(userParam,'userInf')
        % Używamy wartości podanej przez użytkownika
        I0_all(param.startCont) = userParam.userInf;
    else
        % Domyślnie 10 tys.
        I0_all(param.startCont) = 1e4;
    end
    % -----------------------------------------------------------------

    % Składamy wektor Y0 = [S1; I1; R1; D1; S2; I2; R2; D2; ... S4; I4; R4; D4]
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

    % Definicja przedziału czasu
    tSpan = [0, param.tEnd];

    % Wywołanie solvera ODE
    [tSol, YSol] = ode45(@(t,Y) pandemic_odes4(t, Y, param), tSpan, Y0);

    % Rozpakowanie wyników: 
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

    % Rysowanie wyników: 4 pod-wykresy (R, S+I+R, I, D) ===
    figure('Name','Symulacja pandemii - 4 kontynenty','Color','w');
    
    % (A) Wykres wyleczonych (R)
    subplot(2,2,1);
    plot(t_out, R(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, R(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, R(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, R(:,4), 'm-', 'LineWidth', 2);
    xlabel('Czas [dni]');
    ylabel('Wyleczeni R_i(t)');
    legend({'Eurazja','Afryka','Ameryka','Oceania'}, 'Location','best');
    grid on; 
    title('Wyleczeni (R)');

    % (B) Wykres populacji żyjącej (S + I + R)
    subplot(2,2,2);
    plot(t_out, S(:,1)+I(:,1)+R(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, S(:,2)+I(:,2)+R(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, S(:,3)+I(:,3)+R(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, S(:,4)+I(:,4)+R(:,4), 'm-', 'LineWidth', 2);
    xlabel('Czas [dni]');
    ylabel('Populacja żyjąca');
    legend({'Eurazja','Afryka','Ameryka','Oceania'}, 'Location','best');
    grid on; 
    title('Liczba żywych (S+I+R)');

    % (C) Wykres liczby zarażonych (I)
    subplot(2,2,3);
    plot(t_out, I(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, I(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, I(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, I(:,4), 'm-', 'LineWidth', 2);
    xlabel('Czas [dni]');
    ylabel('Zarażeni I_i(t)');
    legend({'Eurazja','Afryka','Ameryka','Oceania'}, 'Location','best');
    grid on; 
    title('Zarażeni (I)');

    % (D) Wykres liczby zgonów (D)
    subplot(2,2,4);
    plot(t_out, D(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_out, D(:,2), 'r-', 'LineWidth', 2);
    plot(t_out, D(:,3), 'g-', 'LineWidth', 2);
    plot(t_out, D(:,4), 'm-', 'LineWidth', 2);
    xlabel('Czas [dni]');
    ylabel('Zgony D_i(t)');
    legend({'Eurazja','Afryka','Ameryka','Oceania'}, 'Location','best');
    grid on; 
    title('Zgony (D)');

end



% POMOCNICZA FUNKCJA: Równania ODE dla 4 kontynentów
function dYdt = pandemic_odes4(t, Y, param)
    % Y = [S1; I1; R1; D1; S2; I2; R2; D2; S3; I3; R3; D3; S4; I4; R4; D4]
    % Rozpakuj:
    S1 = Y(1);  I1 = Y(2);  R1 = Y(3);  D1 = Y(4);
    S2 = Y(5);  I2 = Y(6);  R2 = Y(7);  D2 = Y(8);
    S3 = Y(9);  I3 = Y(10); R3 = Y(11); D3 = Y(12);
    S4 = Y(13); I4 = Y(14); R4 = Y(15); D4 = Y(16);

    % Z parametrow:
    beta1  = param.beta(1);   beta2  = param.beta(2);
    beta3  = param.beta(3);   beta4  = param.beta(4);

    gamma1 = param.gamma(1);  gamma2 = param.gamma(2);
    gamma3 = param.gamma(3);  gamma4 = param.gamma(4);

    mu1    = param.mu(1);     mu2    = param.mu(2);
    mu3    = param.mu(3);     mu4    = param.mu(4);

    % parametry: b_i (narodziny) i d_i (inne zgony)
    if isfield(param,'birth')
        b1 = param.birth(1); b2 = param.birth(2);
        b3 = param.birth(3); b4 = param.birth(4);
    else
        b1 = 0; b2 = 0; b3 = 0; b4 = 0;
    end

    if isfield(param,'otherDeath')
        d1 = param.otherDeath(1);
        d2 = param.otherDeath(2);
        d3 = param.otherDeath(3);
        d4 = param.otherDeath(4);
    else
        d1 = 0; d2 = 0; d3 = 0; d4 = 0;
    end

    alpha = param.alpha;  % macierz 4x4 migracji

    % Populacje żywe:
    N1 = S1 + I1 + R1;
    N2 = S2 + I2 + R2;
    N3 = S3 + I3 + R3;
    N4 = S4 + I4 + R4;

    % --------------------------
    % Kontynent 1 (Eurazja)
    % --------------------------
    infec_1 = beta1 * (S1*I1 / max(N1,1));  
    recov_1 = gamma1*I1;           
    death_1 = mu1*I1;              

    % Migracje - S:
    migS_in_1  = alpha(2,1)*S2 + alpha(3,1)*S3 + alpha(4,1)*S4;
    migS_out_1 = alpha(1,2)*S1 + alpha(1,3)*S1 + alpha(1,4)*S1;
    % Migracje - I:
    migI_in_1  = alpha(2,1)*I2 + alpha(3,1)*I3 + alpha(4,1)*I4;
    migI_out_1 = alpha(1,2)*I1 + alpha(1,3)*I1 + alpha(1,4)*I1;
    % Migracje - R:
    migR_in_1  = alpha(2,1)*R2 + alpha(3,1)*R3 + alpha(4,1)*R4;
    migR_out_1 = alpha(1,2)*R1 + alpha(1,3)*R1 + alpha(1,4)*R1;

    % Narodziny i inne zgony
    births_1   = b1 * N1;     % nowo narodzeni -> S
    otherD_1_S = d1 * S1;     % inne zgony z S
    otherD_1_I = d1 * I1;     % inne zgony z I
    otherD_1_R = d1 * R1;     % inne zgony z R

    dS1 = + births_1 ...
          - infec_1 ...
          - otherD_1_S ...
          + migS_in_1 - migS_out_1;

    dI1 = + infec_1 ...
          - recov_1 ...
          - death_1 ...
          - otherD_1_I ...
          + migI_in_1 - migI_out_1;

    dR1 = + recov_1 ...
          - otherD_1_R ...
          + migR_in_1 - migR_out_1;

    dD1 = + death_1 ...
          + (otherD_1_S + otherD_1_I + otherD_1_R);

    % --------------------------
    % Kontynent 2 (Afryka)
    % --------------------------
    infec_2 = beta2 * (S2*I2 / max(N2,1));
    recov_2 = gamma2*I2;
    death_2 = mu2*I2;

    migS_in_2  = alpha(1,2)*S1 + alpha(3,2)*S3 + alpha(4,2)*S4;
    migS_out_2 = alpha(2,1)*S2 + alpha(2,3)*S2 + alpha(2,4)*S2;
    migI_in_2  = alpha(1,2)*I1 + alpha(3,2)*I3 + alpha(4,2)*I4;
    migI_out_2 = alpha(2,1)*I2 + alpha(2,3)*I2 + alpha(2,4)*I2;
    migR_in_2  = alpha(1,2)*R1 + alpha(3,2)*R3 + alpha(4,2)*R4;
    migR_out_2 = alpha(2,1)*R2 + alpha(2,3)*R2 + alpha(2,4)*R2;

    births_2   = b2 * N2;
    otherD_2_S = d2 * S2;
    otherD_2_I = d2 * I2;
    otherD_2_R = d2 * R2;

    dS2 = + births_2 ...
          - infec_2 ...
          - otherD_2_S ...
          + migS_in_2 - migS_out_2;

    dI2 = + infec_2 ...
          - recov_2 ...
          - death_2 ...
          - otherD_2_I ...
          + migI_in_2 - migI_out_2;

    dR2 = + recov_2 ...
          - otherD_2_R ...
          + migR_in_2 - migR_out_2;

    dD2 = + death_2 ...
          + (otherD_2_S + otherD_2_I + otherD_2_R);

    % --------------------------
    % Kontynent 3 (Ameryka)
    % --------------------------
    infec_3 = beta3 * (S3*I3 / max(N3,1));
    recov_3 = gamma3*I3;
    death_3 = mu3*I3;

    migS_in_3  = alpha(1,3)*S1 + alpha(2,3)*S2 + alpha(4,3)*S4;
    migS_out_3 = alpha(3,1)*S3 + alpha(3,2)*S3 + alpha(3,4)*S3;
    migI_in_3  = alpha(1,3)*I1 + alpha(2,3)*I2 + alpha(4,3)*I4;
    migI_out_3 = alpha(3,1)*I3 + alpha(3,2)*I3 + alpha(3,4)*I3;
    migR_in_3  = alpha(1,3)*R1 + alpha(2,3)*R2 + alpha(4,3)*R4;
    migR_out_3 = alpha(3,1)*R3 + alpha(3,2)*R3 + alpha(3,4)*R3;

    births_3   = b3 * N3;
    otherD_3_S = d3 * S3;
    otherD_3_I = d3 * I3;
    otherD_3_R = d3 * R3;

    dS3 = + births_3 ...
          - infec_3 ...
          - otherD_3_S ...
          + migS_in_3 - migS_out_3;

    dI3 = + infec_3 ...
          - recov_3 ...
          - death_3 ...
          - otherD_3_I ...
          + migI_in_3 - migI_out_3;

    dR3 = + recov_3 ...
          - otherD_3_R ...
          + migR_in_3 - migR_out_3;

    dD3 = + death_3 ...
          + (otherD_3_S + otherD_3_I + otherD_3_R);

    % --------------------------
    % Kontynent 4 (Oceania)
    % --------------------------
    infec_4 = beta4 * (S4*I4 / max(N4,1));
    recov_4 = gamma4*I4;
    death_4 = mu4*I4;

    migS_in_4  = alpha(1,4)*S1 + alpha(2,4)*S2 + alpha(3,4)*S3;
    migS_out_4 = alpha(4,1)*S4 + alpha(4,2)*S4 + alpha(4,3)*S4;
    migI_in_4  = alpha(1,4)*I1 + alpha(2,4)*I2 + alpha(3,4)*I3;
    migI_out_4 = alpha(4,1)*I4 + alpha(4,2)*I4 + alpha(4,3)*I4;
    migR_in_4  = alpha(1,4)*R1 + alpha(2,4)*R2 + alpha(3,4)*R3;
    migR_out_4 = alpha(4,1)*R4 + alpha(4,2)*R4 + alpha(4,3)*R4;

    births_4   = b4 * N4;
    otherD_4_S = d4 * S4;
    otherD_4_I = d4 * I4;
    otherD_4_R = d4 * R4;

    dS4 = + births_4 ...
          - infec_4 ...
          - otherD_4_S ...
          + migS_in_4 - migS_out_4;

    dI4 = + infec_4 ...
          - recov_4 ...
          - death_4 ...
          - otherD_4_I ...
          + migI_in_4 - migI_out_4;

    dR4 = + recov_4 ...
          - otherD_4_R ...
          + migR_in_4 - migR_out_4;

    dD4 = + death_4 ...
          + (otherD_4_S + otherD_4_I + otherD_4_R);

    % Zwracamy wektor pochodnych
    dYdt = [
        dS1; dI1; dR1; dD1;
        dS2; dI2; dR2; dD2;
        dS3; dI3; dR3; dD3;
        dS4; dI4; dR4; dD4
    ];
end