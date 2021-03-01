function BW = find_BW(y, Fs)
% This function removes damped oscillations and the linear trend of force data, 
% and returns the mean weight after removal.
%
% USAGE
% E.g.:
% BW = find_BW(y, 1000)
%
% INPUT
% y:        Force data stored in a array (n x 1).
%
% Fs:       Sampling frequency of the force data.
%
% OUTPUT
% BW:       Scalar with the optimized mean weight of the force data.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 19.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

x = [1/Fs:1/Fs:length(y)/Fs]';                                                  % Time
y = y;                                                                          % Force

y_detrended = detrend(y);                                                       % Remove Linear Trend
yu = max(y_detrended);
yl = min(y_detrended);
yr = (yu-yl);                                                                   % Range of ‘y’
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                            % Returns approximate zero-crossing indices of argument vector
zt = x(zci(y_detrended));
per = 2*mean(diff(zt));                                                         % Estimate period

ym = mean(y_detrended);                                                         % Estimate offset
fit = @(b,x)  b(1) .* exp(b(2).*x) .* (sin(2*pi*x./b(3) + 2*pi/b(4))) + b(5);   % Objective function to fit
fcn = @(b) norm(fit(b,x) - y_detrended);                                        % Least-squares cost function
options = optimset('Display', 'off');                                           % Optimizer options
[s,~] = fminsearch(fcn, [yr; -10;  per;  -1;  ym], options);                    % Minimise least-squares   
BW = mean(y-fit(s,x));                                                          % Ouput

% Debug plot
verbose = 0;
if(verbose == 1)
    subplot(3,1,1)
    plot(x,y_detrended,'b', 'LineWidth',1.5)
    hold on
    plot(x,fit(s,x), '--r')
    hold off
    grid
    xlabel('Time')
    ylabel('Amplitude')
    legend('Original Data',  'Fitted Curve')
    text(0.3*max(xlim),0.7*min(ylim),...
        sprintf('$y = %.3f\\cdot e^{%.0f\\cdot x}\\cdot sin(2\\pi\\cdot x\\cdot %.0f%.3f)$',...
        [s(1:2); 1./s(3:4)]), 'Interpreter','latex')

    subplot(3,1,2)
    y_temp = y-fit(s,x);
    plot(x,y_temp)
    yline(mean(y_temp));
end
end