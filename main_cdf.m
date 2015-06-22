% clear all;
% Reset rand
reset(RandStream.getDefaultStream);
% Pathloss exponent
beta = 6;
% Cell radius
radius = 500;%
% Cutoff parameter
r_cut = 25;
% Cellular setting
Bandwidth = 20*10^(6); N_0 = 1.38*10^(-23)*300;
% Noise variance (-97.8 dBm)
sigma_2 = Bandwidth*2*N_0;
% Carrier frequency in MHz
f_c = 2.4*10^3;
% Attenuation in dB (reference in the paper)
L_dB = -55.9 + 38*log10(r_cut) + (24.5 + f_c/616.67)*log10(f_c);
% Number of users
numOfUsers_vec = 8;
% Ratio K/N
c = 0.5;
% Precoding scheme OLP = 1; MRT = 2; ZF = 3; RZF = 4
precoder_type = 4;
% Rate_max -- Rate_min
rate_max = 1.5; rate_min = 1.5;

% Number of iterations over initial positions
nIter = 1000;
% spatial step in meter
dtheta = 50;
% time step in minute
dtstep = 0.5;
% Time interval T in minutes
max_time = 0.5*(24*60);
% Diffusion coefficient
diffusion_coefficient = dtheta^2/(4*dtstep);
% Number of steps
num_Tsteps = round(max_time/dtstep);

% Energy vector
energy_vec = [];
% Energy mean vector
mean_energy_teo = zeros(length(numOfUsers_vec), 1);
% Energy variance vector
var_energy_teo = zeros(length(numOfUsers_vec), 1);

% Loop over the number of users (specified in numOfUsers_vec)
for ljk = 1:length(numOfUsers_vec)
    
    % Number of users
    numOfUsers = numOfUsers_vec(ljk);
    % Number of antennas
    numOfAnt = numOfUsers/c;
    
    % Rate constraints
    user_rates = rate_min + (rate_max - rate_min)*rand(numOfUsers, 1);
    
    % SNR constraints
    gamma = 2.^user_rates - 1;
    
    % Computing the normalization parameter eta
    switch precoder_type
        
        case 1
            % Optimal linear precoding (OLP)
            eta = 1 - c*mean(gamma./(1 + gamma));
        case 2
            % Maximum ration transmit (MRT)
            eta = 1 - c*mean(gamma);
        case 3
            % Zero-forcing (ZF)
            eta = 1 - c;
        case 4
            % Regularized zero-forcing (RZF)
            eta = 1 - c*(mean(gamma)/( 1 + mean(gamma)));
    end
    
    % Constant value to normalize pathloss
    norm_value = 2*r_cut^beta*(10^(-0.1*L_dB));
    
    % Mean value of pathloss function
    E_x_pathloss = (radius^beta*(2/(beta+2) + r_cut^beta/(radius^beta)))/norm_value;
    % Theoretical energy average
    mean_energy_teo(ljk) = max_time*mean(gamma)*(c*sigma_2/eta)*E_x_pathloss;
    
    % Theoretical energy variance
    var_energy_teo(ljk) = mean(gamma.^2)*(c*sigma_2/eta)^2*...
        var_analytic_2D(max_time, beta, radius, diffusion_coefficient)/norm_value^2;
    % Normalizing variance to number of users
    var_energy_teo(ljk) = var_energy_teo(ljk)/numOfUsers;
    
    tic; counter = 0;
    % Loop over initial positions
    for ii = 1:nIter
        
        % Create user positions over time using random walk.
        r_vec_out = brownian_motion(max_time, dtheta, dtstep, radius, numOfUsers);
        % Average channel gains (pathloss)
        chanGain_vec = norm_value./(r_vec_out.^beta+r_cut^beta);
        
        % Initializing power vector per channel realization
        pow_per_channel_realization = zeros(num_Tsteps,1);
        
        % Loop over channel realizations
        for tt = 1:num_Tsteps
            % Average changain per small-scale fading realization
            chanGain = chanGain_vec(tt,:);
            % Small-scale fading
            W = sqrt(0.5)*(randn(numOfAnt, numOfUsers)+1i*randn(numOfAnt, numOfUsers));
            % Channel matrix
            H = W*sqrt(diag(chanGain));
            % Asymptotic precoder design (OLP, MRT, ZF, RZF)
            [lambda, tx_pow] = precoder_design(precoder_type, chanGain, gamma.', sigma_2, c, eta);
            
            % Designing precoder matrix V
            switch precoder_type
                
                case 1
                    % Optimal linear precoding (OLP)
                    V = (H/(H'*H + numOfAnt*(diag(1./lambda))))*(diag(1./lambda));
                case 2
                    % Maximum ratio transmit (MRT)
                    V = H/numOfAnt;
                case 3
                    % Zero-forcing (ZF)
                    V = H/(H'*H);
                case 4
                    % Optimal regularization parameter
                    rho = 1/mean(gamma) - c/( 1 + mean(gamma));
                    % Regularized zero-forcing (RZF)
                    V = (W*W' + numOfAnt*rho*eye(numOfAnt))\H;
            end
            
            % Power loading
            V = V*sqrt(diag(tx_pow));
            % Power consumption per channel realization
            pow_per_channel_realization(tt) = trace(V*V');
            %
            %             % User of interest
            %             userIndex = 1;
            %             % Channel and Precoding vector user k
            %             h_k = H(:, userIndex); v_k = V(:, userIndex);
            %             % Interference over user k
            %             I_k = real(h_k'*(V*V')*h_k - abs(h_k'*v_k)^2);
            %             % Increasing counter
            %             counter = counter + 1;
            %             % Current SNR of user k
            %             rate_k(ljk, counter) = log2(1 + abs(h_k'*v_k)^2/(I_k + sigma_2));
        end
        
        % Computing the energy consumption at time = dtstep*num_Tsteps
        energy_vec(ii, ljk) = dtstep*sum(pow_per_channel_realization); %#ok<*SAGROW>
        
        if round(ii/100) == ii/100
            toc;
            x = ii %#ok<NOPTS>
        end
        
    end
    % Energy mean
    mean_energy(ljk) = mean(energy_vec(:,ljk));
    % Variance energy
    var_energy(ljk) = var(energy_vec(:,ljk));
end