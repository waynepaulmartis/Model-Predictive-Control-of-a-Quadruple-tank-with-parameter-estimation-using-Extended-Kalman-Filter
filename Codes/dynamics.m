%% 
% *DYNAMICS OF THE SYSTEM*

syms x3dot x4dot x3 x4 x1 x2 b1 b2 b3 b4 b5 g1 g2 u1 u2 c1
P = [b1+b2*cos(x2) b3+b4*cos(x2); b3+b4*cos(x2) b5]
Q = -c1*sin(x2)*[x4 x3+x4; -x3 0]
R = [g1*cos(x1) + g2*cos(x1+x2); g2*cos(x1+x2)]

eqn=P*[x3dot; x4dot] == -Q*[x3; x4] - R + [u1; u2]
%%
S = solve(eqn, [x3dot; x4dot])
%%
S.x3dot
S.x4dot