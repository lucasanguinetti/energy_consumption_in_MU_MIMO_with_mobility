function [lambda, tx_pow] = precoder_design(precoder_type, chanGain, gamma, sigma_2, c, eta);

% Number of users
numOfUsers = length(chanGain);
% Lambda coefficients
lambda = zeros(numOfUsers, 1);
% Transmit power vector
tx_pow = zeros(numOfUsers, 1);
% Parameter A (Average of the inverse of average channel attenuation)
A = mean(gamma./chanGain);
% Average transmit power
avg_pow = c*sigma_2/eta*A;
% Average of SINR constraints
avg_gamma = mean(gamma);

switch precoder_type
    
    case 1 % OLP
        % kth lambda coefficient
        lambda = gamma./(chanGain*(eta));
        % kth transmit power
        tx_pow = gamma./(chanGain*(eta)^2).*...
            (avg_pow + sigma_2./chanGain.*( 1 + gamma).^2);
        
    case 2 % MRT
        % kth transmit power
        tx_pow = gamma./chanGain.*(avg_pow + sigma_2./chanGain);
        
    case 3 % ZF
        % kth transmit power
        tx_pow = gamma.*sigma_2;
    case 4% RZF
        % kth transmit power
        tx_pow = gamma./(chanGain.*avg_gamma^2).*...
            (avg_pow + sigma_2./chanGain*( 1 + avg_gamma)^2);
        
end