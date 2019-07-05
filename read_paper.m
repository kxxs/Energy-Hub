F_CHP = sdpvar;
Q_WARG = sdpvar;
Q_CHP = sdpvar;
W_CHP = sdpvar;
R_WARG = sdpvar;
F_in = sdpvar;
R_d = sdpvar;
Q_d = sdpvar;
W_d = sdpvar;

eta_Q = 0.25;
eta_W = 0.5;
eta_R = 0.8;
H_CHP = [eta_Q 1 0; % F Q W
         eta_W 0 1];
H_WARG = [eta_R 1];
A_CHP = [1 0 0 0 0;   % F
         0 -1 -1 0 0; % Q
         0 0 0 -1 0;]; % W
A_WARG = [0 1 0 0 0;
          0 0 0 0 -1];
Z_CHP = H_CHP*A_CHP;
Z_WARG = H_WARG*A_WARG;
Z = [Z_CHP; Z_WARG];
X = [1 0 0 0 0];
Y = [0 0 0 0 1;
    0 0 1 0 0;
    0 0 0 1 0];
Cost = 100*F_in;
Cons = [[X;Y;Z]*[F_CHP;Q_WARG;Q_CHP;W_CHP;R_WARG] == [F_in; R_d; Q_d; W_d; 0; 0; 0];
        R_d >= 50;
        Q_d >= 30;
        W_d >= 100;
        F_CHP >=0; Q_WARG>=0; Q_CHP>=0;W_CHP>=0;R_WARG>=0];
ops = sdpsettings('solver','gurobi','verbose',1);
double(F_in)
double(Cost)
solvesdp(Cons,Cost,ops);
result.cost = double(Cost);
result.F_CHP = double(F_CHP);
result.Q_WARG = double(Q_WARG);
result.Q_CHP = double(Q_CHP);