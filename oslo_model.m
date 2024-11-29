% Oslo model

set(groot, 'defaultTextInterpreter' , 'latex')
set(groot, 'defaultLegendInterpreter' , 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')

L = [100 200 400];
[H, Z_T, Z] = oslo_model_simulation(L);

function [H, Z_T, Z] = oslo_model_simulation(L)

    S_all = [];

    
    for l = 1:length(L)
    H = zeros(1, L(l)+1);
    Z_T = randi([1 2], 1, L(l));
    num_iterations = L(l)^2;

    for t = 1:num_iterations
        
        H(1) = H(1) + 1;
        Z = calculate_zi(H);
        unstable = true;
        while unstable
            unstable = false;
            for i = 1:L(l)
                if Z(i) > Z_T(i)
                    unstable = true;

                    H(i) = H(i) - 1;
                    H(i+1) = H(i+1) + 1;
                    H(end) = 0;
                    Z_T(i) = randi([1 2]);
                    Z = calculate_zi(H);
                end
            end
        end
    end

    % Plot results
    figure;
    xaxis = 1:L(l);
    label = "$L=" + L(l) + "$";
    plot(xaxis, H(1, 1:end-1), 'DisplayName', label)
    xlim([1 length(H)-1])
    xlabel('n')
    ylabel('h')
    legend
    grid on
    title('h(n)')

    % Size of avalanche measurement

    num_iterations = 1000;
    S = zeros(1, num_iterations);
    S_max = length(H)-1;

    for t = 1:num_iterations
        H(1) = H(1) + 1;
        Z = calculate_zi(H);
        % Relax
        unstable = any(Z > Z_T);
        while unstable
            unstable = false;
            for i = 1:L(l)
                if Z(i) > Z_T(i)
                    unstable = true;
                    S(t) = S(t) + 1;
                    H(i) = H(i) - 1;
                    H(i+1) = H(i+1) + 1;
                    Z_T(i) = randi([1 2]);
                    Z = calculate_zi(H);
                end
            end
        end
    end

    % Plot S/Smax to t
    S_all = [S_all, S];
    t = 1:num_iterations;
    figure;
    label = "$L=" + L(l) + "$";
    stem(t, S/S_max, 'Marker','none', 'DisplayName', label)
    xlabel('t')
    ylabel('$\frac{S}{S_{max}}$')
    title('$\frac{S}{S_{max}}(t)$')
    legend
    end

    % Probablity P(S,L)

    function [counts, unique_sizes] = count_unique_sizes(S)
    unique_sizes = unique(S);
    counts = histcounts(S, [unique_sizes, max(unique_sizes)+1]);
    end

    [n_of_occurence_S, unique_S] = count_unique_sizes(S_all);
    P = n_of_occurence_S/sum(n_of_occurence_S);

    figure;
    loglog(unique_S, P)
    xlabel('S')
    ylabel('P(S,L)')
    ylim([min(P) max(P)])
    title('P(S,L)(S)')
    grid on;

    figure;
    loglog(unique_S(1:94), P(1:94))
    xlabel('S')
    ylabel('P(S,L)')
    ylim([min(P) max(P)])
    title('P(S,L)(S)')
    grid on;

    figure;
    logfit(unique_S(1:94), P(1:94), 'loglog', 'markertype', '.', 'markersize', 5.5, 'linewidth', 1.5)
    xlabel('S')
    ylabel('P(S,L)')
    title('P(S,L)(S)')
    grid on;
end

function Z = calculate_zi(H)
    Z = H(1:end-1) - H(2:end);
end