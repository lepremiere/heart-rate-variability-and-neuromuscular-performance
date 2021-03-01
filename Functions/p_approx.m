function [a, b, r_crit] = p_approx(r, r0, n, alpha, method, mode)
% This function approximates the probability for a correlation coefficient
% beeing inferior to a certain true population effect.
%
% USAGE
% E.g.:
% [a, b, r_crit] = p_approx(0.77, 0, 42, 0.05, 'Fisher', 'superiority')
%
% INPUT
% r:        Correlation coefficient that should be tested
%           against r0. 
% r0:       The Null-hypothesis to be tested against.
% n:        Number of observation in determining r. 
%           Minimum n = 3.
% alpha:    Confidence level.
% method:   Character vector that determines which method 
%           is used to calculate "a, b, r_crit". Options:
%           'Exact', uses the exact distribution. 
%           'Fisher', uses Fisher's z-transformation approximates.
% mode:     Character vector that determines which test
%           should be performed. Options: 'equivalence', 
%           tests if r is different from r0. 'superiority',
%           tests if r is greater than r0. 'inferiority', 
%           tests if r is smaller than r0.
%
% OUTPUT
% a:        Alpha error probability.
% b:        Beta error probability.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 01.February.2021
% lepremiere
%---------------------------------------------------------------------------------------------------
% Symplifying estimations by using symmetry property
r   = abs(r);   
r0  = abs(r0);

if(r == 1)
    r = 0.99999;
elseif(n < 3)          % Minimum of df = 3 to return valid results
    a = NaN;
    b = NaN;
    return;
end

precision   = 1e-4;    % Determines the precision for the approximation
s           = r > r0;  % "Sideness" for equivalence testing

switch method
    case 'Fisher'
    % Z-transformation and standar error
    r0_z    = atanh(r0);                                                  
    r_z     = atanh(r);
    se      = 1/sqrt(n-3);
    
    switch mode
        case 'equivalence'
            % Critical r, alpha and beta error probability 
            % for two-sided test
            z_crit      = norminv(abs(1-alpha/2));
            r_crit(1,1) = tanh(r0_z - z_crit * se);
            r_crit(2,1) = tanh(r0_z + z_crit * se); 
            a           = abs(s - normcdf((r_z - r0_z)/se));
            b           = sum(normcdf((r_z - atanh(r_crit))/se));
        case 'superiority'
            % Critical r, alpha and beta error probability 
            % for one-sided test
            z_crit      = norminv(abs(1-alpha));
            r_crit      = tanh(r0_z + z_crit * se);
            a           = 1 - normcdf((r_z - r0_z)/se);
            b           = abs(s - normcdf((r_z - atanh(r_crit))/se));
        case 'inferiority'
            % Critical r, alpha and beta error probability 
            % for one-sided test
            z_crit      = norminv(abs(1-alpha));
            r_crit      = tanh(r0_z - z_crit * se);
            a           = normcdf((r_z - r0_z)/se);
            b           = normcdf((r_z - atanh(r_crit))/se);
    end
    
    case 'Exact'
        x      = -1:precision:1;                   % Array for which the probability density function should be calculated
        r0_pdf = corrdist(x, r0, n);               % Exact probability density function for r0
        r0_pdf = r0_pdf./trapz(r0_pdf);            % Normalizing
        r0_cum = cumtrapz(r0_pdf);                 % Cummulative density function for r0

        r_pdf  = corrdist(x, r, n);                % Exact probability density function for r
        r_pdf  = r_pdf./trapz(r_pdf);              % Normalizing
        r_cum  = cumtrapz(r_pdf);                  % Cummulative density function for r

        switch mode
            case 'equivalence'
                urcrit  = find(r0_cum >= 1-alpha/2, 1);     % Upper bound of h0
                lrcrit  = find(r0_cum >= alpha/2, 1);       % Lower bound of h0
                mrcrit  = find(r_cum >= 0.5, 1);
                a       = abs(s - r0_cum(mrcrit));
                b       = trapz(r_pdf(lrcrit:urcrit));      
            case 'inferiority'
                urcrit  = size(x,2);                        % Upper bound of h0
                lrcrit  = find(r0_cum >= alpha, 1);         % Lower bound of h0
                mrcrit  = find(r_cum >= 0.5, 1);
                a       = r0_cum(mrcrit);                   % Probability of cpdf_h1 at r
                b       = 1 - r_cum(lrcrit);   
            case 'superiority'
                urcrit  = find(r0_cum >= 1-alpha, 1);       % Upper bound of h0
                lrcrit  = 1;                                % Lower bound of h0
                mrcrit  = find(r_cum >= 0.5, 1);
                a       = abs(s - r0_cum(mrcrit));          % Probability of cpdf_h1 at r
                b       = r_cum(urcrit);   
        end
        r_crit = x([lrcrit, urcrit]);

end

% Debug plot
verbose = 0; 
if(verbose)
    figure(1);
    clf;
    hold on;

    idx0    = find(r0_cum >= 0.5, 1, 'first');
    idx     = find(r_cum >= alpha/2, 1, 'first');
    idx2    = find(r_cum >= 1-alpha/2, 1, 'first');
    
    h = area(x(lrcrit:urcrit),r_pdf(lrcrit:urcrit), 'LineStyle', 'none');
    h.FaceColor = 'b';
    h.FaceAlpha = 0.3;
    
    h = area(x(1:lrcrit),r0_pdf(1:lrcrit));
    h.FaceColor = 'r';
    h.FaceAlpha = 0.3;
    h = area(x(urcrit:end),r0_pdf(urcrit:end));
    h.FaceColor = 'r';
    h.FaceAlpha = 0.3;
    
    plot(x,r0_pdf,'r', 'LineWidth', 2);
    xline(x(idx0), 'r', 'LineWidth', 2);

    plot(x,r_pdf,'--b', 'LineWidth', 2);
    plot(x(idx:idx2), max(r_pdf)/2.*ones(length(x(idx:idx2)),1), 'b')
    
    xline(x(idx0), 'r', 'LineWidth', 2);
    xline(x(mrcrit),'b', 'LineWidth', 2);
    xline(x(urcrit),'g', 'LineWidth', 2);
    xline(x(lrcrit),'g', 'LineWidth', 2);
    
    title("Crit r: " + num2str(x(lrcrit))+ "/" + num2str(x(urcrit))...
            + " P-Value: " + num2str(a) + " beta: " + num2str(b) + ...
            ' Mean r: ' + num2str(x(mrcrit)) + ' Mean r0: ' + num2str(x(idx0)));          
end
end