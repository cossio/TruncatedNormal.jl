addpath([pwd, '/mvrandn'], [pwd, '/trandn'])
format long

X = mvrandn([100, 200], [120, 230], [1 0.9; 0.9 1], 1000000);

mean(X, 2)  % <x1> and <x2>

mean(X(1,:) .* X(1,:))  % <x1^2>
mean(X(2,:) .* X(2,:))  % <x2^2>
mean(X(1,:) .* X(2,:))  % <x1*x2>

