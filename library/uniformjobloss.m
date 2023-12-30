function [g_new] = uniformjobloss(g, L_new)
[I, ~] = size(g);
g_new = zeros(I,2);
L_old = sum(g(:,2));
g_new(:,2) = L_new/L_old*g(:,2);
g_new(:,1) = g(:,1)+g(:,2) - g_new(:,2);



