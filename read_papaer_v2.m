Time = 1:24;
          %1    2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
          %17   18  19  20  21  22  23  24
Solar  = [ 0,   0,  0,  0,  0,  3,  5,  16, 26, 34, 38, 40, 42, 40, 36, 26,...
           12,  4,  1,  0,  0,  0,  0,  0;
          ]*2.5;
Demand = [45,   47, 46, 44, 45, 53, 56, 77, 80, 92, 95, 94, 85, 83, 84, 82,...
          89,   101,110,115,110,90, 70, 50;   % W (electricity)
          17,   16, 16, 16, 16, 18, 20, 25, 25, 27, 30, 37, 35, 34, 33, 32....
          30,   29, 29, 29, 27, 23, 23, 24;   % R (cold)
          30,   30, 30, 30, 32, 32, 36, 40, 42, 35, 32, 28, 16, 17, 21, 23,...
          28,   32, 25, 26, 27, 30, 32, 35;   % Q (head)
          ]; 
Gas_In = zeros(1,24);
Elec_In = zeros(1,24);

for hour = 1:24
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
      Vout_W1 + Vout_W2 >= Demand(1,hour);
      Vout_R1 + Vout_R2 >= Demand(2,hour);
      Vout_Q1 + Vout_Q2 + Vout_Q3 >= Demand(3,hour);
    ];
    Cost = 220*(Vin_W1+Vin_W2) + 2*(Vin_W1+Vin_W2)^2 + 100*(Vin_F1 + Vin_F2) + (Vin_F1 + Vin_F2)^2;

    ops = sdpsettings('solver','gurobi','verbose',1);
    solvesdp(Cons,Cost,ops);

    Elec_In(hour) = double(Vin_W1+Vin_W2);
    Gas_In(hour) = double(Vin_F1+Vin_F2);
    result(hour).Cost = double(Cost);
    result(hour).V = double(V);
end
figure()
plot(Time,Elec_In)
hold on 
plot(Gas_In)