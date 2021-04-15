// Standard Normal variate using Box-Muller transform.
function randn_bm() {
    var u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
}


console.log("Welcome");

//Simulate true signal and measurement
//Define system

//number of cycles to run
N = 1000;

//time step, delta t
dt = 0.001;

//time vector containing all time stamps in seconds
t = []
for (i = 0; i < N; i++) {
    t[i] = dt * (i + 1);
    //console.log(t[i]);
}


//define matrices

//system matrix - state. 
F = math.matrix([[1, dt], //Corresponds to coefficients for h_{k-1} and hdot_{k-1} in position kinematic eqn
    [0, 1]]);             //Corresponds to coefficients for h_{k-1} and hdot_{k-1} in vel kinematic eqn

//system matrix - input    
G = math.matrix([[-.5 * dt * dt], //Corresponds to coefficient for g in pos eqn
    [-dt]]);                    //Corresponds to coefficient for g in vel eqn

//observation matrix  
H = math.matrix([[1, 0]]);

//process noise covariance
Q = math.matrix([[0, 0], [0, 0]]);

//input = accel due to gravity (in m/s^2)
u = 9.80665;

//2x2 identity matrix
I = math.identity(2);
console.log("Done");


//define init position and vel
y0 = 100; //m
v0 = 0; //m/s

//initialize state vector (true state) by taking blank 2xN matrix
xt = math.zeros(2, N); //true state vec

//set 0th column of xt to y0, v0
xt = math.subset(xt, math.index(math.range(0, 2), 0), [y0, v0]) //true initial state

console.log("xt first", xt);

aDueToGravTerm = math.multiply(G, u);

console.log(aDueToGravTerm);
console.log(G);
console.log(u);

test = math.subset(xt, math.index(math.range(0, 2), 0));

console.log("xt subset is", math.subset(xt, math.index(math.range(0, 2), 0)));
console.log("f is", F);
console.log("size of f is", math.size(F), "and size of subset is", math.size(test));


math.multiply(F, math.subset(xt, math.index(math.range(0, 2), i - 1)));

console.log("size of G is", math.size(G));


//go through all columns of the matrix, replacing with the calculated state
for (i = 1; i < N; i++) {
    //propogate the states through the prediction equation for our entire time range
    xt = math.subset(xt, 
        //select current column (current state slot)
        math.index(math.range(0, 2), i), 
        //what to replace the selected subset with: the predicted state for that time, using equation x_{k} = F_{k-1}x_{k-1} + G_{k-1}u_{k-1}
        math.add(math.multiply(F, math.subset(xt, math.index(math.range(0, 2), i - 1))), 
                aDueToGravTerm)); 
}

console.log("Final xt is", xt);


//Generate noisy measurement matrix from true state

//Measurement noise covariance mat reduces to scalar val, since measurement system has stdev of 2m (variance 4 m^2), and there's only one term in output vec
R = 4; //in m^2

//create 1xN zeros matrix
v = math.zeros(1, N);

//replace all zeros in top row of v with random numbers
for (i = 0; i < N; i++) {
    v = math.subset(v, math.index(math.range(0, 1), i), randn_bm());
}

//matrix of the measurement noises at each timestamp
v = math.multiply(v, Math.sqrt(R)); //measurement noises

//multiply the system matrix H by the state matrix, which just selects the position component of each predicted state. Then add the measurement noise
z = math.add(math.multiply(H, xt), v); //the calculated (simulated) noisy measurement

//our fake measurements for the experiment are now in z



//Perform Kalman filter estimation
//Init the state vector (est state)
x = math.zeros(2, N); //estimated state vec

//guess for initial state (height of 105m), while the actual initial position is 100m
x = math.subset(x, math.index(math.range(0, 2), 0), [105, 0]);

//Initialize covariance matrix: error of 10 m^2 for initial position, assume error of 0.01 m^2/s^2 for initial vel (obj at rest). No cross terms. 
P = math.matrix([[10, 0], [0, 0.01]]); //Covariance mat for initial state error

//Loop thru and perform Kalman filter equations recursively
for (i = 1; i < N; i++) {
    //KALMAN FILTER ALGORITHM uses a prediction followed by a correction to determine states of the filter
    //Using info about dynamics of the state, the filter will project forward and predict next state
    //The correction/update part involves comparing a measurement w/what we predict that measurement should have been based on predicted states
    //Predictor-corrector format is applied recursively at each time step

    //PREDICT

    //Predict state vector, same as we did above, using state dynamic equation
    x = math.subset(xt, 
        //select current column (current state slot)
        math.index(math.range(0, 2), i), 
        //what to replace the selected subset with: the predicted state for that time, using state dynamic equation x_{k} = F_{k-1}x_{k-1} + G_{k-1}u_{k-1}
        math.add(math.multiply(F, math.subset(xt, math.index(math.range(0, 2), i - 1))), 
                aDueToGravTerm));


    //Predict state error covariance matrix. We multiply the system state mat (F) by the previous state's measurement error covariance matrix by transpose of F, 
    //then add process noise covar mat (0 in this case)

    //Equation: P_{k|k-1} = F_{k-1}P_{k-1}F_{k-1}^T + Q_{k-1}
    P = math.add(math.multiply(math.multiply(F, P), math.transpose(F)), Q);

    //Calculate Kalman gain matrix
    //H is a matrix necessary to define the output equation
    //R is measurement noise covariance
    K = math.multiply(math.multiply(P, math.transpose(H)), math.inv(math.add(math.multiply(math.multiply(H, P), math.transpose(H)), R)));

    //console.log(K);

    //CORRECT
    
    //Update state vector by scaling the "innovation" (diff between measurement of output, z_{k}, and predicted output, H_{k}xhat_{k|k-1}) by calculated Kalman gain mat,
    //in order to correct prediction by appropriate amount: xhat_{k} = xhat_{k|k-1} + K_{k} * (z_{k} - H_{k}xhat_{k|k-1})

    //get fake measurement for this timestamp
    measuredOutput = math.subset(z, math.index(math.range(0, 1), i)); //where z is a 1xN matrix of measurements
    console.log("Measured output is", measuredOutput);

    //get the predicted state vector for this timestamp, as predicted above
    currStateCol = math.subset(x, math.index(math.range(0, 2), i));

    //calculate predicted output (predicted position)
    predictedOutput = math.multiply(H, currStateCol);
    
    //calculate how much we should add to the predicted state vector to correct it
    toAdd = math.multiply(K, math.subtract(measuredOutput, predictedOutput));

    x = math.subset(x, math.index(math.range(0, 2), i), math.add(currStateCol, toAdd));

    //Update state error covariance: P_{k} = (I - K_{k}H_{k})P_{k|k-1}
    P = math.multiply(math.subtract(I, math.multiply(K, H)), P);
}



/*
//Plot results
//Plot states
figure(1);
subplot(211);
plot(t,z,'g-',t,x(1,:),'b--','LineWidth',2);
hold on; plot(t,xt(1,:),'r:','LineWidth',1.5)
xlabel('t (s)'); ylabel('x_1 = h (m)'); grid on;
legend('Measured', 'Estimated', 'True');
subplot(212);
plot(t, x(2,:), 'b--', 'LineWidth', 2);
hold on; plot(t, xt(2,:), 'r:', 'LineWidth', 1.5)
xlabel('t (s)'); ylabel('x_2 = v (m/s)'); grid on;
legend('Estimated', 'True');


//Plot estimation errors
figure(2);
subplot(211);
plot(t, x(1,:)-xt(1,:),'m','LineWidth', 2)
xlabel('t (s)'); ylabel('\Deltax_1 (m)');grid on;
subplot(212);
plot(t, x(2,:)-xt(2,:), 'm', 'LineWidth', 2)
xlabel('t (s)'); ylabel('\Deltax_2 (m/s)'); grid on;*/