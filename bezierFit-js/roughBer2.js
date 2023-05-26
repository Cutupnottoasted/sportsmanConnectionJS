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

// require will not work in adobe illustrator, a translation function must be made 
// translation, matrix product, and inverse needed
const math = require('mathjs'); 

// function aff_angle(d) // aff_angle(x)
// {
//     // n = size(X, 1);
//     var n = d.length;
//     // Xbar = X - ones(size(X)) * diag(mean(X));
//     var means = [];
//     for (var i = 0; i < d[0].length; i++) // d[0].length gives number of columns
//     {
//         var sum = 0;
//         for (var j = 0; j < n; j++)
//             sum += d[j][i];
    
//         var mean = sum / n;
//         means.push(mean);
//         console.log("column " + (i + 1) + " mean: " + mean);
//     }
    
    
//     var onesMatrix = [];
//     for (var i = 0; i < n; i++) 
//     {
//         var row = [];
//         for (var j = 0; j < d[0].length; j++)
//             row.push(1);
        
//         onesMatrix.push(row);
//     }
    
//     console.log("matrix of ones:");
//     console.log(onesMatrix);
    
//     var diagMatrix = []; 
//     for (var i = 0; i < d[0].length; i++)
//     {
//         var row = [];
//         for (var j = 0; j < d[0].length; j++)
//         {
//             if (i === j)
//                 row.push(means[i]);
//             else
//                 row.push(0);
//         }
//         diagMatrix.push(row);
//     }
    
//     console.log("diagonal matrix with mean values on diagonal: ");
//     console.log(diagMatrix);
    
//     var Xbar = [];
//     for (var i = 0; i < n; i++)
//     {
//         var row = [];
//         for (var j = 0; j < d[0].length; j++)
//             row.push(d[i][j] - onesMatrix[i][j] * diagMatrix[j][j]);
        
//         Xbar.push(row);
//     }
    
//     console.log("The matrix representing the deviation of each element in d from the column-wise mean");
//     console.log(Xbar); 

//     // Xcov = Xbar' * Xbar/n;
//     var XbarT = math.transpose(Xbar);
//     console.log("Xbar matrix translated:");
//     console.log(XbarT);

//     n = Xbar.length;
//     var Xcov = math.multiply(XbarT, Xbar).map(row => math.divide(row, n));
//     console.log("Xcov matrix:");
//     console.log(Xcov);

//     // A = inv(Xcov);
//     var A = math.inv(Xcov); 
//     console.log("A inverse matrix:");
//     console.log(A); 
    
//     // V = X(2 : n, :) - X(1 : n - 1, :);
//     // X(2:n, :) selects submatrix of X starting from second row up to last
//     // including all the columns ':' AKA extracts all rows of X except the first
//     //
//     // X(l:n - 1, :) selects a submatrix of X stating from l-th row up to n - 1 row
//     // including al lthe columns. Excluding the first l row and the last l row
//     //
//     // x(2:n, :) - X(l:n-1, :) subtracts the two submatrices calculated above, by element
//     // the subtraction is performed for each corresponding element in the matrices

//     n = A.length;
//     var l = 1; 

//     var V = d.slice(1, n).map((row, i) =>
//         row.map((value, j) => value - d[i + l][j]));
        
//     console.log("The submatricies subtracted element-wise");
//     console.log(V);

//     // t = diag(V * A * V').^ (1/2);
//     // diag extracts the diagonal elements of the resulting matrix
//     // .^ element wise exponential operator
//     var product = math.multiply(V, A, math.transpose(V));
//     var diagonal = math.diag(product);
//     t = math.map(diagonal, value => math.pow(value, 1/2)); 
//     console.log("matrix multiplication, and extracted diagonal elements (^ 1/2):")
//     console.log(t);

//     // M^2[X](X_(i - 1), X(i + 1))
//     // obtain these values

//     // V_2 = X(3:n, :) - X(1:n-2, :);
//     n = d.length; 
//     var V_2 = d.slice(2, n).map((row, i) => 
//         row.map((value, j) => value - d[i][j]));
//     console.log("Slice and subtracing from 3rd row to last row from d excluding the first n-2 rows");
//     console.log(V_2);

//     // t_2 = diag(V_2 * A * V_2')
//     // diag extracts the diagonal elements of resulting matrix
//     product = math.multiply(V_2, A, math.transpose(V_2));
//     var t_2 = math.diag(product); 
//     console.log("resulting matrix product after diagonal element are extracted");
//     console.log(t_2);

//     // find the theta_i values (takes into account the bending of data)
//     // theta = zeros(n - 1, 1); 
    
// }

function aff_angle(d)
{
    var n = d.length;
    const bar = find_bar(d, n);
    console.log("Find bar matrix");
    console.log(bar);
    const A = find_A(bar, n);
    console.log("A matrix: ");
    console.log(A);
    // V = X(2: n, :) - X(1: n - 1, :)
    const V = m_difference(d, n, 1);
    console.log("V matrix: ");
    console.log(V);
    // t = diag(V * A * V') .^ (1/2)
    var temp = m_multiply(V, A);
    temp = m_multiply(temp, transpose(V));
    temp = m_diag(temp);
    const t = m_sqroot(temp);
    console.log("t matrix: ")
    console.log(t);
    // V_2 = X(3: n, :) - X(1:n-2,:);
    const V_2 = m_difference(d, n, 2);
    console.log("V_2 matrix: ");
    console.log(V_2);
    // t_2 = diag(V_2 * A * V_2'); 
    temp = m_multiply(V_2, A);
    temp = m_multiply(temp, transpose(V_2)); 
    const t_2 = m_diag(temp);
    console.log("t_2 matrix: ");
    console.log(t_2);
    
    // theta = zeros(n - 1, 1);
    var theta = Array(n - 1).fill(0);
    // for j = 2:n - 1...
    for (var j = 1; j <= n - 2; j++) {
        theta[j] = Math.min(Math.PI - Math.acos((t[j - 1] - 2 + t[j] - 2 - t_2[j - 1]) 
                                                    / (2 * t[j] * t[j - 1])), Math.PI / 2);
    }
    console.log("theta matrix: ");
    console.log(theta);

    // h = zeros(n - 1, 1);
    var h = Array(n - 1).fill(0);
    // h(1) = t(1)*( 1 + (1.5*theta(2)*t(2))/(t(1) + t(2)) );
    var h_1 = t[0] * (1 + (1.5 * theta[1] * t[1]) / (t[0] + t[1]));
    console.log(h_1);
    // for j = 2:n-2
    var hj;
    for (let j = 1; j <= n - 3; j++) {
        h[j] =
          t[j] *
          (1 +
            (1.5 * theta[j] * t[j - 1]) / (t[j - 1] + t[j]) +
            (1.5 * theta[j + 1] * t[j + 1]) / (t[j] + t[j + 1]));
        if (h[j] != 0)
            h_j = h[j];
      }
    console.log(h_j);

    // h(n-1) = t(n-1) * ( 1 + (1.5*theta(n-1)*t(n-2))/(t(n-2)+t(n-1)) ); 
    var h_3 = t[n - 2] * (1 + (1.5 * theta[n - 2] * t[n - 3]) / (t[n - 3] + t[n - 2]));
    console.log(h_3);

    // missing just the last 3 lines of code

}



// bar = d - ones(size(d)) * diag(mean(d))
function find_bar(d, n)
{    
    // calculate the mean of each column
    var col_means = find_col_means(d);
    console.log("Find column means: "); 
    console.log(col_means);
    // create a diagonal matrix from col_means
    var diagMatrix = create_diag_matrix(col_means);
    console.log("find diagonal matrix: ");
    console.log(diagMatrix);
    
    // create a matrix of ones based on the original matrix
    var onesMatrix = Array.from({ length: n }, () => Array(d[0].length).fill(1));
    // find the product of the diagonal means matrix and the matrix of ones
    productMatrix = matrixProduct(onesMatrix, diagMatrix);
    // subtract the initial with the productMatrix element wise
    return bar = d.map((row, i) => row.map((value, j) => value - productMatrix[i][j]));
}

function find_col_means(d)
{
    var rows = d.length;
    var columns = d[0].length;
    var means = Array(columns).fill(0);

    for (var i = 0; i < columns; i++) 
    {
        var sum = 0;
        for (var j = 0; j < rows; j++)
        {
            sum += d[j][i];
        }
        means[i] = sum / rows; 
    }
    return means;
}

function create_diag_matrix(m)
{
    var size = m.length;
    diagMatrix = Array(size).fill(0).map(() => Array(size).fill(0));
    for (var i = 0; i < size; i++) 
    {
        diagMatrix[i][i] = m[i];
    }
    return diagMatrix;
}

function matrixProduct(m1, m2)
{
    var result = [];
    var numRows1 = m1.length,
        numCols1 = m1[0].length,
        numRows2 = m2.length,
        numCols2 = m2[0].length;
    
    if (numCols1 !== numRows2)
        throw new Error("Number of columns in m1 must be equal to number of rows in m2");
    
    for (var i = 0; i < numRows1; i++) {
        result[i] = [];
        for (var j = 0; j < numCols2; j++) {
            let sum = 0;
            for (let k = 0; k < numCols1; k++) {
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

// Xcov = Xbar' * Xbar/n
// A = inv(Xcov)
// find the covariance of bar matrix and returns the inverse
function find_A(bar, n)
{
    var barT = transpose(bar);
    var barDiv = m_divide(bar, n);
    var cov= m_multiply(barT, barDiv);
    var inverse = inv(cov);
    return inverse;
}

// transposes and returns matrix m
function transpose(m)
{
    var result = m[0].map((_, i) => m.map(row => row[i]));
    return result;
}
// finds the inverse of matrix m
function inv(m)
{
    var det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    var inverse = [
        [m[1][1] / det, -m[0][1] / det],
        [-m[1][0] / det, m[0][0] / det]
    ];
    return inverse;
}
// divides matrix m by n
function m_divide(m, n)
{
    var result = m.map(row => row.map(value => value / n));
    return result;
}
// finds matrix product of m1 and m2
function m_multiply(m1, m2)
{
    var result = [];
    for (var i = 0; i < m1.length; i++) {
        result[i] = [];
        for (var j = 0; j < m2[0].length; j++) {
            result[i][j] = 0;
            for (var k = 0; k < m2.length; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}
// finds the difference from the same matrix
// where i is = n - i
function m_difference(m, n, i)
{
    var result = [];
    for (var j = i; j < n; j++) {
        result.push([
            m[j][0] - m[j - i][0],
            m[j][1] - m[j - i][1]
        ]);
    }
    return result;
}
// returns the diagonal of a matrix
function m_diag(m)
{
    var result = [];
    for (var i = 0; i < m.length; i++)
        result.push(m[i][i]);
    return result;
}

function m_sqroot(m)
{
    var result = []
    result = m.map(element => element ** (1/2));
    return result;
}
// main

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

// console.log("returned deviation matrix as var t: ");
// console.log(t);

// var bez_mat = mxbern2(t, deg);

    
// Xbar = X - ones(size(X)) * diag(mean(X)); 
