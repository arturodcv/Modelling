% Gene regulatory circuit

max_time = 100;                 % Maximum time for the simulation
gil = gillespie(max_time);      % Gillespie algorithm
q = qssa(max_time);             % QSSA algorithm
leaping = tau_leap(max_time);   % Tau-leap algorithm

% Kolmogorov-Smirnov test for the results of gillespie and qssa algorithms
significant_alpha = 0.05;       % Significant level for the test

[h,p_value] = kstest2(gil(200:10000),q(200:10000),'alpha',significant_alpha);

if h == 1
    fprintf(['The null hypothesis H0: "The data from gillespie and qssa algorithms ' ...
        'follows the same continuous distribution" has been rejected under\n the ' ...
        'Kolmogorov-Smirnov test with a %f %% sigficance level and a p-value of %f.' ], ...
        significant_alpha, p_value)
else 
    fprintf(['The null hypothesis H0: "The data from gillespie and qssa algorithms' ...
        'follows the same continuous distribution" has been accepted under the \n' ...
        'Kolmogorov-Smirnov test with a %f %% sigficance level and a p-value of %f.' ], ...
        significant_alpha, p_value)
end

% Transition rates
function [value] = transition_rates(x,R_hat,k1,k2,E,u11,b11)
    value = [R_hat + k1*x(2); k2*x(1);b11*x(1)*(x(1)-1)*(E-x(2));u11*x(2)];
end

% Select the process
function [jota] = which_process_occur (z2,w)
    for j = 1:4 
        upper = sum(w(1:j));
        if z2*sum(w) < upper
            jota = j;
            break
        end
    end    
end

% Gillespie algorithm
function [x1] = gillespie(T)
    rng(sum(1000000*clock)); % Initialisation of the random number generator
    x1 = zeros(1,2000);      % Variable we will study the statistics on
    w_time = zeros(1,2000);  % Waiting time
    
    % Parameters 
    S = 500; E = 5; a = 3; r = 0.4; Kd = 10; kdeg = 2; b11 = kdeg/(E*S);
    R = r/(kdeg*sqrt(Kd)); v11 = 1; kappa1 = a/(kdeg*sqrt(Kd)); kappa2 = 1;
    
    % Parameters
    k1 = kappa1*b11*S^2; k2 = kappa2*kdeg; u11 = v11*b11*S^2; R_hat = R*kdeg*S;
    
    % X1 is number of transcription factor molecules
    % X2 is number of bound promoter sites in the gene promoter region
    x = [0,0];                 % Initial values for variables
    t = 0;                     % Initial time
    r = [1,0;-1,0;-2,1;2,-1];  % Stochastic dynamics of each process
    
    iter = 0;
    while t<T
        z1 = rand; z2 = rand;
        w = transition_rates(x,R_hat,k1,k2,E,u11,b11);
        w0 = sum(w);
        tau = log(1/z1)/w0;
        j = which_process_occur(z2,w);
        t = t+tau;
        x = x + r(j,:); 
        iter = iter + 1;
        x1(iter) = x(1); 
        w_time(iter) = tau;
    end
    figure(1)
    plot(x1(200:iter));
    figure(2)
    histogram(w_time(200:iter));
    figure(3)
    boxplot(x1(200:iter))
end

% QSSA algorithm
function [x1] = qssa(T)
    rng(sum(1000000*clock)); % Initialisation of the random number generator
    x1 = zeros(1,2000);      % Variable we will study the statistics on
    w_time = zeros(1,2000);  % Waiting time
    
    % Parameters 
    S = 500; E = 5; a = 3; r = 0.4; Kd = 10; kdeg = 2; b11 = kdeg/(E*S);
    R = r/(kdeg*sqrt(Kd)); v11 = 1; kappa1 = a/(kdeg*sqrt(Kd)); kappa2 = 1;
    
    % Parameters
    k1 = kappa1*b11*S^2; k2 = kappa2*kdeg; u11 = v11*b11*S^2; R_hat = R*kdeg*S;
    
    % X1 is number of transcription factor molecules
    % X2 is number of bound promoter sites in the gene promoter region
    x = [0,0];                 % Initial values for variables
    t = 0;                     % Initial time
    r = [1,0;-1,0;-2,1;2,-1];  % Stochastic dynamics of each process
    
    iter = 0;
    while t<T
        z1 = rand; z2 = rand;
        w = transition_rates(x,R_hat,k1,k2,E,u11,b11); 
        w0 = sum(w); 
        tau = log(1/z1)/w0; 
        j = which_process_occur(z2,w);
        t = t+tau;
        x(1) = x(1) + r(j,1); 
        p = (x(1)/S)^2/((x(1)/S)^2+u11);
        x(2) = binornd(E,p);
        iter = iter + 1;
        x1(iter) = x(1);
        w_time(iter) = tau;
    end
    figure(4)
    plot(x1(200:iter));
    figure(5)
    histogram(w_time(200:iter));
    figure(6)
    boxplot(x1(200:iter))
end


% Derivative of the transition rates with respect to x
function [value] = w_derivate(x,k1,k2,E,u11,b11)
    value = [0, k2, b11*(2*x(1)-1)*(E-x(2)), 0; k1, 0, -b11*x(1)*(x(1)-1),u11];
end

% Tau-leap algorithm
function [x1] = tau_leap(T)
    rng(sum(1000000*clock)); % Initialisation of the random number generator
    x1 = zeros(1,2000);      % Variable we will study the statistics on
    w_time = zeros(1,2000);  % Waiting time
    
    % Parameters 
    S = 500; E = 5; a = 3; r = 0.4; Kd = 10; kdeg = 2; b11 = kdeg/(E*S);
    R = r/(kdeg*sqrt(Kd)); v11 = 1; kappa1 = a/(kdeg*sqrt(Kd)); kappa2 = 1;
    epsilon = E/S;
    
    % Parameters
    k1 = kappa1*b11*S^2; k2 = kappa2*kdeg; u11 = v11*b11*S^2; R_hat = R*kdeg*S;
    
    % X1 is number of transcription factor molecules
    % X2 is number of bound promoter sites in the gene promoter region
    x = [0,0];                 % Initial values for variables
    t = 0;                     % Initial time
    r = [1,0;-1,0;-2,1;2,-1];  % Stochastic dynamics of each process
    
    iter = 0;
    while t<T
        w = transition_rates(x,R_hat,k1,k2,E,u11,b11);
        w0 = sum(w); 
        xi = sum(w.*r);
        dw = w_derivate(x,k1,k2,E,u11,b11);
        tau = min(min(epsilon*w0./abs(xi*dw)));
        deltax = sum(r.*w)*tau;
        x = x + deltax;
        t = t+tau;
        iter = iter + 1; 
        x1(iter) = x(1);
        w_time(iter) = tau;
    end
    figure(7)
    plot(x1(250:iter));
    figure(8)
    histogram(w_time(200:iter));
    figure(9)
    boxplot(x1(200:iter))
end


    
