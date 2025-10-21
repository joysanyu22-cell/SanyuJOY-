clc; clear;

% Application to real-world situations
m = 80; g = 9.81; c = 0.25;
f_terminal = @(v) m*g - c*v.^2;
df_terminal = @(v) -2*c*v;
true_v = sqrt(m*g/c);

% Secant Method for terminal velocity
function [v2, secant_terminal_vals] = secant_terminal_recursive(f, v0, v1, n, secant_terminal_vals)
    if n == 0
        v2 = v1;
    else
        v2 = v1 - f(v1)*(v1 - v0)/(f(v1) - f(v0));
        secant_terminal_vals(end+1) = v2;
        if abs(v2 - v1) < 1e-6
            fprintf('\nTerminal Velocity Secant≈ %.6f (to 6 d.p.)\n', v2);
        else
            [v2, secant_terminal_vals] = secant_terminal_recursive(f, v1, v2, n-1, secant_terminal_vals);
        end
    end
end

% Newton-Raphson Method for terminal velocity
function [v_new, newton_terminal_vals] = newton_terminal_recursive(f, df, v, n, newton_terminal_vals)
    if n == 0
        v_new = v;
    else
        v_new = v - f(v)/df(v);
        newton_terminal_vals(end+1) = v_new;
        if abs(v_new - v) < 1e-6
            fprintf('\nTerminal Velocity NRM≈ %.6f (to 6 d.p.)\n', v_new);
        else
            [v_new, newton_terminal_vals] = newton_terminal_recursive(f, df, v_new, n-1, newton_terminal_vals);
        end
    end
end

% Measure computation time for terminal velocity
tic;
[v2, secant_terminal_vals] = secant_terminal_recursive(f_terminal, 30, 40, 20, [30, 40]);
secant_terminal_time = toc;

tic;
[v_new, newton_terminal_vals] = newton_terminal_recursive(f_terminal, df_terminal, 30, 20, 30);
newton_terminal_time = toc;

% Results
fprintf('True Terminal Velocity: %.4f m/s\n', true_v);