% Code to calculate the maximum and minimum complexity curves

D = 6; % Number of points to form the symbol
N = factorial(D); % Number of possible symbols in the PDF
num_pontos = 200; % Number of points to form the curves

PDF_eq = ones(1,N)* (1/N); % Equiprobable PDF
PDF_deg = zeros(1,N); % PDF with only one symbol/state of maximum probability
PDF_deg(1) = 1; % Assigning maximum probability to the first symbol of PDF_deg
S_max = log2(N); % Maximum entropy
% Normalization constant of complexity
Q_0 = -2*(((N+1)/N)*log2(N+1) - 2*log2(2*N) + log2(N))^-1;

% Definition of the canonical basis (a vector with 1 at the i-th position and 0 elsewhere)
e = zeros(N,N);
for i = 1:N
    e(i,i) = 1;
end

%%%%%%%%% Calculation of Minimum Complexity

H_min = zeros(num_pontos, 1); % Create matrix to store entropy
C_min = zeros(num_pontos, 1); % Create matrix to store minimum complexity
id = 1;  % Index for matrices H and C

% Loop over 'rho' from 0 to 1/N, with 'num_pontos' equally spaced values
for rho = linspace(0, 1/N, num_pontos)
    % Equation that defines 'P' - the PDF that minimizes complexity
    P = rho * sum(e(1:N-1, :),1) + (1 - rho*(N-1))*e(N,:);
    % Calculates normalized entropy of distribution 'P'
    H_min(id) = S(P)/S_max;
    % Calculates Jensen-Shannon divergence
    D_js = S((P + PDF_eq) / 2) - (1/2) * (S(PDF_eq) + S(P));
    % Calculates statistical complexity of 'P' and stores the values
    C_min(id) = Q_0 * H_min(id) * D_js;
    % Advance index to next point in vectors H and C
    id = id + 1;
end

%%%%%%%%% Calculation of Maximum Complexity

H_max = zeros(num_pontos, 1); % Create matrix to store entropy
C_max = zeros(num_pontos, 1); % Create matrix to store maximum complexity
id = 1;  % Index for matrices H and C

% Loop over L_max paths
for i = 1:N-1
    % Loop over 'rho' from 0 to 1/(N-i+1), with 'num_pontos' equally spaced values
    for rho = fliplr(linspace(0, 1/(N-i+1), round(num_pontos/4)))
        % Equation that defines 'P' - the PDF that maximizes complexity
        P = rho*e(i,:) + ((1.0-rho)/(N-i))*sum(e(i+1:N, :),1);
        % Calculates normalized entropy of distribution 'P'
        H_max(id) = S(P)/S_max;
        % Calculates Jensen-Shannon divergence
        D_js = S((P + PDF_eq) / 2) - (1/2) * (S(PDF_eq) + S(P));
        % Calculates statistical complexity of 'P' and stores the values
        C_max(id) = Q_0 * H_max(id) * D_js;
        % Advance index to next point in vectors H and C
        id = id + 1;
    end
end

% Plot the curves
figure(1)
clf
hold on
plot(H_max, C_max, '-k', 'LineWidth',3);
plot(H_min, C_min, '-k', 'LineWidth',3);
ylabel('Complexity','Interpreter','latex', 'FontSize', 30);
xlabel('Entropy','Interpreter','latex', 'FontSize', 30);

% Function to calculate Shannon entropy
function [out] = S(P)
out = 0;
for i = 1:length(P)
    if P(i) == 0
        continue
    end
    out = out - P(i).*log2(P(i));
end
end
