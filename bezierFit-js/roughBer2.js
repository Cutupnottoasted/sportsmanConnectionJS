// var d = []; // an array of objects containing x, y coordinates
// var deg = 3; // the degree of the bezier curve
// var stop = -4; // stoping criteria which the algorithm becomes 10^(stop)

// var p = []; // an array of objects containing x, y control coordinates
// var t = []; // a vector of nodes based on spread and bed of data
//             // vaar t = aff_angle(d); 

// i = size (d, 1); used to assign the number of rows in matrix d to the var i

// page 14 ex. figure 4 on page 18
// degree = 2
// d = ( {1, 1}, {2, 3}, {4, 2}, {5, 5} )

// aff_angle function
// returns affine invariant vector of nodes from ordered data d
// assumes that ordered data d is in the form of x, y coordinates

function aff_angle(d) // aff_angle(x)
{
    var n = d.length;
    var means = [];
    for (var i = 0; i < d[0].length; i++) // d[0].length gives number of columns
    {
        var sum = 0;
        for (var j = 0; j < n; j++)
            sum += d[j][i];
    
        var mean = sum / n;
        means.push(mean);
        console.log("column " + (i + 1) + " mean: " + mean);
    }
    
    
    var onesMatrix = [];
    for (var i = 0; i < n; i++) 
    {
        var row = [];
        for (var j = 0; j < d[0].length; j++)
            row.push(1);
        
        onesMatrix.push(row);
    }
    
    console.log("matrix of ones:");
    console.log(onesMatrix);
    
    var diagMatrix = []; 
    for (var i = 0; i < d[0].length; i++)
    {
        var row = [];
        for (var j = 0; j < d[0].length; j++)
        {
            if (i === j)
                row.push(means[i]);
            else
                row.push(0);
        }
        diagMatrix.push(row);
    }
    
    console.log("diagonal matrix with mean values on diagonal: ");
    console.log(diagMatrix);
    
    var Xbar = [];
    for (var i = 0; i < n; i++)
    {
        var row = [];
        for (var j = 0; j < d[0].length; j++)
            row.push(d[i][j] - onesMatrix[i][j] * diagMatrix[j][j]);
        
        Xbar.push(row);
    }
    
    console.log("The matrix representing the deviation of each element in d from the column-wise mean");
    console.log(Xbar); 
    return Xbar;   
}

function mxbern2(t, deg)
{
    // [n, m] = size (t);
}



// setting matrix of ordered points
var d = [
    [1, 1],
    [2, 3],
    [4, 2],
    [5, 5]
];


console.log("matrix of data points: ");
console.log(d);
console.log("length of d: " + d.length);
console.log("column lenght of d: " + d[0].length);
var deg = 3; // configured for real world
var stop = -4; // configured for real world

// step 1

// skip hold_was_off

var i = d.length;  // number of data points
var j = deg + 1; // number of control points
var t = aff_angle(d);

console.log("returned deviation matrix as var t: ");
console.log(t);

    
// Xbar = X - ones(size(X)) * diag(mean(X)); 
