var canvDim = 600;
var estCtrXInit = canvDim / 2 + canvDim / 4;
var actualCtrXInit = canvDim / 2 - canvDim / 4;
var initialY = 70;
var rangeFinderBoxWidth = 100;
var rangeFinderBoxCtrX = actualCtrXInit;
var rangeFinderBoxLeftX = rangeFinderBoxCtrX - rangeFinderBoxWidth / 2;
var rangeFinderBoxTop = 575;
var ballRad = 15;
var zeroLineYVal = 550;
var g = 9.80665;
var time = 0;
var yPosActual = initialY;
var desiredSec = 10;
var predictedOutputs = [];

//time step, delta t
dt = 0.01;

var xArray = [];
var yArrayActual = [];

//plotly actual position data
var actualPosData = {
    x: xArray,
    y: yArrayActual,
    mode: "lines",
    type: "scatter",
    line: {
        color: 'rgb(55, 128, 191)',
        width: 1
    },
    name: "Actual position"
};

var yArrayEstimated = [];

//plotly estimated position data
var estimatedPosData = {
    x: xArray,
    y: yArrayEstimated,
    mode: "lines",
    type: "scatter",
    line: {
        color: 'rgb(255, 128, 191)',
        width: 1
    },
    name: "Estimated position"
};

var yArrayMeasured = [];
var yArrayMeasuredLive = [];

//plotly simulated measurement data
var measuredPosData = {
    x: xArray,
    y: yArrayMeasuredLive,
    mode: "lines",
    type:"scatter",
    line: {
        color: 'rgb(23, 255, 191)',
        width: 1
    },
    name: "\'Measured\' position"
};

var allData = [measuredPosData, actualPosData, estimatedPosData];

//plotly graph layout
var layout = {
    autosize: false,
    width: 750,
    height: 300,
    margin: {
        l: 50,
        r: 50,
        b: 50,
        t: 50,
        pad: 4
    },
    xaxis: {range: [0, 10], title: "Time (s)"},
    yaxis: {range: [50, 550], title: "Position (pixels)"},
    title: "Position vs. Time Plot"
};

/**
 * The position of the object as it undergoes free fall is estimated despite uncertainty about its initial position and uncertainty
 * in measurements of its position taken with a "laser rangefinder."
 */

//Standard Normal variate using Box-Muller transform.
function randn_bm() {
    var u = 0, v = 0;
    while (u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while (v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u )) * Math.cos(2.0 * Math.PI * v);
}

function setLineDash(list) {
    drawingContext.setLineDash(list);
}

function cleanCanvas() {
    let cnv = createCanvas(canvDim, canvDim);
    background(230, 230, 230);
    cnv.position(800, 100);
}

var setup = function() {
    cleanCanvas();

    //run Kalman simulation for all discrete timepoints once at beginning
    runKalmanSim();

    button = createButton('Click to restart simulation');
    button.position(800, 80);
    button.mouseClicked(onButtonClicked);
}

var draw = function() {
    stroke(0, 0, 0);
    background(230, 230, 230);
    textSize(14);
    text('Actual Newtonian motion', 150, 40);

    textSize(10);
    textAlign(CENTER);
    text('Newtonian motion estimated by applying Kalman filter to simulated rangefinder measurement data', 350, 10, 200, 100);
    rect(rangeFinderBoxLeftX, rangeFinderBoxTop, 100, 100);

    time = frameCount * dt;
    textSize(15);
    text('Rangefinder', 150, 590);
    text('Time elapsed: ' + round(time, 3) + ' s', 520, 590);

    //calculate y pos of actual falling ball, using kinematic equation
    yPosActual = initialY + 0 * time + 0.5 * g * time * time;
    yPosActual = constrain(yPosActual, initialY, zeroLineYVal - ballRad);

    xArray.push(time);
    yArrayActual.push(yPosActual);

    indexer = frameCount - 1;
    indexer = constrain(indexer, 0, predictedOutputs.length - 1);

    yPosEstimated = constrain(predictedOutputs[indexer], initialY, zeroLineYVal - ballRad);

    yArrayEstimated.push(yPosEstimated);

    yPosMeasured = constrain(yArrayMeasured[indexer], initialY, zeroLineYVal - ballRad);

    yArrayMeasuredLive.push(yPosMeasured);

    push();
    fill(55, 128, 191);
    ellipse(actualCtrXInit, yPosActual, 2 * ballRad);
    push();
    fill(255, 128, 191);
    setLineDash([4, 4]);
    ellipse(estCtrXInit, yPosEstimated, 2 * ballRad);
    pop();
    pop();

    line(0, zeroLineYVal, canvDim, zeroLineYVal);

    stroke(255, 0, 0);
    line(rangeFinderBoxCtrX, rangeFinderBoxTop, rangeFinderBoxCtrX, yPosEstimated + ballRad);

    if (yPosActual < zeroLineYVal - ballRad) {
        Plotly.newPlot("posvtimeplot", allData, layout);
    }
}

function resetDataArrays() {
    xArray = [];
    yArrayActual = [];

    //plotly actual position data
    actualPosData = {
        x: xArray,
        y: yArrayActual,
        mode: "lines",
        type: "scatter",
        line: {
            color: 'rgb(55, 128, 191)',
            width: 1
        },
        name: "Actual position"
    };

    yArrayEstimated = [];

    //plotly estimated position data
    estimatedPosData = {
        x: xArray,
        y: yArrayEstimated,
        mode: "lines",
        type: "scatter",
        line: {
            color: 'rgb(255, 128, 191)',
            width: 1
        },
        name: "Estimated position"
    };

    yArrayMeasuredLive = [];

    //plotly simulated measurement data
    var measuredPosData = {
        x: xArray,
        y: yArrayMeasuredLive,
        mode: "lines",
        type:"scatter",
        line: {
            color: 'rgb(23, 255, 191)',
            width: 1
        },
        name: "\'Measured\' position"
    };

    allData = [measuredPosData, actualPosData, estimatedPosData];
}

var onButtonClicked = function() {
    console.log("Button cicked!");
    frameCount = 0;
    cleanCanvas();
    resetDataArrays();
}

var mouseClicked = function() {
}

function runKalmanSim() {
    console.log("Welcome to Kalman Sim fxn");

    //number of cycles to run
    N = desiredSec / dt;

    //time vector containing all discrete time stamps for simulation (in seconds)
    t = []
    for (i = 0; i < N; i++) {
        t[i] = dt * (i + 1); //first timestamp is (dt)sec, last is (dt * N)sec
    }

    //define matrices

    //system matrix - state. This matrix gets multiplied to the previous state's matrix
    F = math.matrix([[1, dt], //Corresponds to coefficients for h_{k-1} and hdot_{k-1} in position kinematic eqn
        [0, 1]]);             //Corresponds to coefficients for h_{k-1} and hdot_{k-1} in vel kinematic eqn

    //system matrix - input. This matrix gets multiplied to g, the graviational constant 
    G = math.matrix([[.5 * dt * dt], //Corresponds to coefficient for g in pos eqn
        [dt]]);                    //Corresponds to coefficient for g in vel eqn

    //observation matrix - basically just used to select position from state matrix
    H = math.matrix([[1, 0]]);

    //process noise covariance matrix - can be set to zero since assume there's minimal inherent process error in the kinematic equations
    Q = math.matrix([[0, 0], [0, 0]]);

    //input vector = accel due to gravity (in m/s^2) -- a constant in this example
    u = g;

    //2x2 identity matrix
    I = math.identity(2);

    //define ACTUAL init position and vel
    y0 = initialY; //m
    v0 = 0; //m/s

    //initialize ACTUAL (assuming perfect system w/no noise) state vector array by taking blank 2xN matrix
    xt = math.zeros(2, N); 

    //set 0th column of xt to [y0 v0], which is the ACTUAL true initial state of the system
    xt = math.subset(xt, math.index(math.range(0, 2), 0), [y0, v0]);

    //console.log("xt first", xt);

    //multiply the G matrix by u to get the constant G*u term for the Kalman state estimation equation
    aDueToGravTerm = math.multiply(G, u);


    /*------DBUG
    console.log(aDueToGravTerm);
    console.log(G);
    console.log(u);

    test = math.subset(xt, math.index(math.range(0, 2), 0));

    console.log("xt subset is", math.subset(xt, math.index(math.range(0, 2), 0)));
    console.log("f is", F);
    console.log("size of f is", math.size(F), "and size of subset is", math.size(test));

    math.multiply(F, math.subset(xt, math.index(math.range(0, 2), i - 1)));

    console.log("size of G is", math.size(G));
    -----END DBUG*/


    //go through all columns of the matrix, replacing with the calculated state
    for (i = 1; i < N; i++) {
        //propogate the states through the prediction equation for our entire time range
        xt = math.subset(xt, 
            //select current column (current state slot)
            math.index(math.range(0, 2), i), 
            //what to replace the selected subset with: the predicted state for that time, using equation x_{k} = F_{k-1}x_{k-1} + G_{k-1}u_{k-1}
            math.add(math.multiply(F, math.subset(xt, math.index(math.range(0, 2), i - 1))), aDueToGravTerm)); 
    }

    //print out final array
    //console.log("Final xt is", xt);

    //Generate noisy measurement matrix from true state

    //Measurement noise covariance mat reduces to scalar val, since measurement system has stdev of 2m (variance 4 m^2), and there's only one term in output vec
    R = 4; //in m^2

    //create 1xN zeros matrix for v, the measurement noise vector
    v = math.zeros(1, N);

    //replace all entries of v with random numbers to simulate noise
    for (i = 0; i < N; i++) {
        v = math.subset(v, math.index(math.range(0, 1), i), randn_bm());
    }

    //matrix of the measurement noises at each timestamp
    v = math.multiply(v, Math.sqrt(R)); //final measurement noises are now in v

    //multiply the system matrix H by the state matrix, which just selects the position component of each predicted state. Then add the measurement noise
    z = math.add(math.multiply(H, xt), v); //the calculated (simulated) noisy measurement

    /**
     * Our fake laser rangefinder measurements for the experiment are now in z. All we did was compute the array of actual values of the object's position
     * at each timestep according to Newtonian physics, then added our measurement noise to it.
    */



    //NEXT: perform the Kalman filter estimation of the object's position at each timestep
    //initialize ESTIMATED state vector array by taking blank 2xN matrix
    x = math.zeros(2, N); 

    //simulate guess for estimated initial state, while the ACTUAL initial position is initialY
    x = math.subset(x, math.index(math.range(0, 2), 0), [initialY + (0.05 * initialY), 0]);

    //initialize covariance matrix: error of 10 m^2 for initial position (since it was a guess
    //(obj is p much guaranteed to be rest). No cross terms. 
    P = math.matrix([[10, 0], [0, 0.01]]); //covar mat for initial state error

    //Loop thru and perform Kalman filter equations recursively
    for (i = 1; i < N; i++) {
        //KALMAN FILTER ALGORITHM uses a prediction followed by a correction to determine states of the filter
        //Using info about dynamics of the state, the filter will project forward and predict next state
        //The correction/update part involves comparing a measurement w/what we predict that measurement should have been based on predicted states
        //Predictor-corrector format is applied recursively at each time step

        //-------PREDICT-------

        //Predict state vector, same as we did above, using state dynamic equation  *****CHECK xt vs. x
        x = math.subset(x, 
            //select current column (current state slot)
            math.index(math.range(0, 2), i), 
            //what to replace the selected subset with: the predicted state for that time, using state dynamic equation x_{k} = F_{k-1}x_{k-1} + G_{k-1}u_{k-1}
            math.add(math.multiply(F, math.subset(x, math.index(math.range(0, 2), i - 1))), aDueToGravTerm));


        //Predict state error covariance matrix. We multiply the system state mat (F) by the previous state's measurement error covariance matrix by transpose of F, 
        //then add process noise covar mat (0 in this case)

        //Equation: P_{k|k-1} = F_{k-1}P_{k-1}F_{k-1}^T + Q_{k-1}
        P = math.add(math.multiply(math.multiply(F, P), math.transpose(F)), Q);

        //Calculate Kalman gain matrix
        //H is a matrix necessary to define the output equation. It's called the observation matrix. In this case it's just used to select position from state matrix
        //R is measurement noise covariance, which is 4 m^2 in this example
        K = math.multiply(math.multiply(P, math.transpose(H)), math.inv(math.add(math.multiply(math.multiply(H, P), math.transpose(H)), R)));

        //console.log(K);

        //--------CORRECT-------
        
        //Update state vector by scaling the "innovation" (diff between measurement of output, z_{k}, and predicted output, H_{k}xhat_{k|k-1}) by calculated Kalman gain mat,
        //in order to correct prediction by appropriate amount: xhat_{k} = xhat_{k|k-1} + K_{k} * (z_{k} - H_{k}xhat_{k|k-1})

        //get fake measurement for this timestamp
        measuredOutput = math.subset(z, math.index(math.range(0, 1), i)); //where z is a 1xN matrix of fake measurements
        //console.log("Fake measurement is", measuredOutput);
        yArrayMeasured.push(measuredOutput);

        //get the predicted state vector for this timestamp (assuming no noise in system), as predicted in for loop above
        currStateCol = math.subset(x, math.index(math.range(0, 2), i));

        //calculate predicted output (predicted position at this timestamp)
        predictedOutput = math.multiply(H, currStateCol);
        
        console.log(predictedOutput.get([0, 0]));
        predictedOutputs[i - 1] = predictedOutput.get([0, 0]);

        //find error between measured position of obj and the predicted value
        error = math.subtract(measuredOutput, predictedOutput);
        
        //calculate how much we should add to the predicted state vector to correct it, by multiplying Kalman gain matrix by error
        toAdd = math.multiply(K, error);

        //apply the correction to the state predictions
        x = math.subset(x, math.index(math.range(0, 2), i), math.add(currStateCol, toAdd));

        //update state error covariance: P_{k} = (I - K_{k}H_{k})P_{k|k-1}
        P = math.multiply(math.subtract(I, math.multiply(K, H)), P);
    }
}