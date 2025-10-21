clear; clc;
% Knapsack Problem

% Recursive approach
function v = knapsackR(W, wt, v, n)
    if n == 0 || W == 0
        v = 0;
    elseif wt(n) > W
        v = knapsackR(W, wt, v, n-1);
    else
        v = max(v(n) + knapsackR(W-wt(n), wt, v, n-1), knapsackR(W, wt, v, n-1));
    end
end
% Dynamic programming approach
function V = knapsackDP(W, wt, v)
    n = length(v);
    V = zeros(n+1, W+1);
    for i = 2:n+1
        for w = 2:W+1
            if wt(i-1) < w
                V(i, w) = max(V(i-1, w), v(i-1) + V(i-1, w-wt(i-1)));
            else
                V(i, w) = V(i-1, w);
            end
        end
    end
    V = V(n+1, W+1);
end

% Fibonacci Sequence
% Recursive approach
function f = fibR(n)
    if n <= 1
        f = n;
    else
        f = fibR(n-1) + fibR(n-2);
    end
end
% Dynamic programming approach
function f = fibDP(n)
    f = zeros(1, n+1);
    f(1) = 0; f(2) = 1;
    for i = 3:n+1
        f(i) = f(i-1) + f(i-2);
    end
    f = f(n+1);
end

% Test and plot
n_values = 10:5:30;
recursive_times_knapsack = zeros(size(n_values));
dp_times_knapsack = zeros(size(n_values));
recursive_times_fib = zeros(size(n_values));
dp_times_fib = zeros(size(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    W = n *10;
    wt = [10 ,20 ,30];
    v = [6,10,12];
    
    tic; knapsackR(W, wt, v, length(v)); recursive_times_knapsack(i) = toc;
    tic; knapsackDP(W, wt, v); dp_times_knapsack(i) = toc;
    tic; fibR(n); recursive_times_fib(i) = toc;
    tic; fibDP(n); dp_times_fib(i) = toc;
end
figure
plot(n_values, recursive_times_knapsack, 'b-o', n_values, dp_times_knapsack, 'g-s');
xlabel('Fibonacci number (n)'); ylabel('Time (s)'); title('Knapsack Problem'); legend('Recursive', 'Dynamic programming');
figure
plot(n_values, recursive_times_fib, 'r-o', n_values, dp_times_fib, 'b-s');
xlabel('Fibonacci number (n)'); ylabel('Time (s)'); title('Fibonacci Sequence'); legend('Recursive', 'Dynamic programming');