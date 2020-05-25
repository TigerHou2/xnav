close all hidden
clear;clc

profile on

disp('Internal syms loop:')
tic
for i = 1:50
    x_val = external_syms();
end
toc

disp(' ')

disp('External syms loop:')
tic
syms x
for i = 1:50
    x_val = external_syms(x);
end
toc

profile viewer


%% Functions

function x_val = internal_syms()
    syms x
    eqn = x == 1;
    soln = solve(eqn,x);
    x_val = double(soln);
end

function x_val = external_syms(x)
    if nargin == 0
        syms x
    end
    eqn = x == 1;
    soln = solve(eqn,x);
    x_val = double(soln);
end