Demand = [100;50;30]; % Demand of W,R,Q
Vin_W1 = sdpvar;
Vin_W2 = sdpvar;
Vin_F1 = sdpvar;
Vin_F2 = sdpvar;
Vout_W1 = sdpvar;
Vout_W2 = sdpvar;
Vout_R1 = sdpvar;
Vout_R2 = sdpvar;
Vout_Q1 = sdpvar;
Vout_Q2 = sdpvar;
Vout_Q3 = sdpvar;
V = sdpvar(11,1);

Cons = [
  Vin_W1 == Vout_W1; 
  Vout_W1 == V(1,1);
  Vin_W2 == V(2,1);
  Vin_F1 == V(3,1);
  Vin_F2 == V(4,1);
  V(5,1) + V(6,1) + V(7,1) == V(3,1)*0.5;
  V(8,1) + V(9,1) == 0.25 * V(3,1);
  V(10,1) + V(11,1) == 0.8 * V(4,1);
  Vout_R1 == V(2,1)*1.2 + V(5,1)*1.2;
  Vout_R2 == V(9,1) * 0.9;
  Vout_Q3 == V(7,1)*0.3;
  Vout_W2 == V(6,1);
  Vout_Q1 == V(8,1);
  Vout_Q2 == V(11,1); V>=0;
  Vout_W1 + Vout_W2 >= Demand(1);
  Vout_R1 + Vout_R2 >= Demand(2);
  Vout_Q1 + Vout_Q2 + Vout_Q3 >= Demand(3);
];
Cost = 220*(Vin_W1+Vin_W2) + 2*(Vin_W1+Vin_W2)^2 + 100*(Vin_F1 + Vin_F2) + (Vin_F1 + Vin_F2)^2;

ops = sdpsettings('solver','gurobi','verbose',1);
solvesdp(Cons,Cost,ops);

result.Cost = double(Cost)
result.V = double(V)
% 5.0416667 * 10^4