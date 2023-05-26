main();
function main()
{
    if (documents.length < 1)
        return; 
    var myDoc = app.activeDocument;
    var myLine = myDoc.pathItems.add();

    //set stroked to true so we can see the path
    myLine.stroked = true;
    var sel = activeDocument.selection;
    if (!(sel instanceof Array) || sel.length < 2) return;
    
    selectPaths(sel, 1, paths_array);
    if (paths_array.length < 1)
        return;

    for(var i = 0; i < paths_array.length; i++)
    {
        draw(path_array[i], myLine); 
    }

    
    
    
}

function draw(path_item, myLine)
{
    var path_point = path_item.pathPoints;
    
    for(var i = 0; i < path_point.length; i++)
    {
        var newPoint = myLine.pathPoints.add();
        newPoint.anchor = path_point[i].anchor;
        newPoint.leftDirection = path_point[i].leftDirection;
        newPoint.rightDirection = path_point[i].rightDirection;
        newPoint.pointType = path_point[i].pointType; 
    }
}

(function () { 'use strict';

// to suit your point format, run search/replace for '.x' and '.y';
// for 3D version, see 3d branch (configurability would draw significant performance overhead)

// square distance between 2 points
function getSqDist(p1, p2) {

    var dx = p1.x - p2.x,
        dy = p1.y - p2.y;

    return dx * dx + dy * dy;
}

// square distance from a point to a segment
function getSqSegDist(p, p1, p2) {

    var x = p1.x,
        y = p1.y,
        dx = p2.x - x,
        dy = p2.y - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((p.x - x) * dx + (p.y - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = p2.x;
            y = p2.y;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = p.x - x;
    dy = p.y - y;

    return dx * dx + dy * dy;
}
// rest of the code doesn't care about point format

// basic distance-based simplification
function simplifyRadialDist(points, sqTolerance) {

    var prevPoint = points[0],
        newPoints = [prevPoint],
        point;

    for (var i = 1, len = points.length; i < len; i++) {
        point = points[i];

        if (getSqDist(point, prevPoint) > sqTolerance) {
            newPoints.push(point);
            prevPoint = point;
        }
    }

    if (prevPoint !== point) newPoints.push(point);

    return newPoints;
}

function simplifyDPStep(points, first, last, sqTolerance, simplified) {
    var maxSqDist = sqTolerance,
        index;

    for (var i = first + 1; i < last; i++) {
        var sqDist = getSqSegDist(points[i], points[first], points[last]);

        if (sqDist > maxSqDist) {
            index = i;
            maxSqDist = sqDist;
        }
    }

    if (maxSqDist > sqTolerance) {
        if (index - first > 1) simplifyDPStep(points, first, index, sqTolerance, simplified);
        simplified.push(points[index]);
        if (last - index > 1) simplifyDPStep(points, index, last, sqTolerance, simplified);
    }
}

// simplification using Ramer-Douglas-Peucker algorithm
function simplifyDouglasPeucker(points, sqTolerance) {
    var last = points.length - 1;

    var simplified = [points[0]];
    simplifyDPStep(points, 0, last, sqTolerance, simplified);
    simplified.push(points[last]);

    return simplified;
}

// both algorithms combined for awesome performance
function simplify(points, tolerance, highestQuality) {

    if (points.length <= 2) return points;

    var sqTolerance = tolerance !== undefined ? tolerance * tolerance : 1;

    points = highestQuality ? points : simplifyRadialDist(points, sqTolerance);
    points = simplifyDouglasPeucker(points, sqTolerance);

    return points;
}


})();


// n = floor limit of pathPoints in pathObject
function getPathItems(n , paths_array)
{
    if (documents.length < 1) return;
    var selection = activeDocument.selection;

    // if pathObject does not extend Array or pathPoints are < 1 return
    if (!(selection instanceof Array) || selection.length < 1) return;
    // else extract all path items from selection
    selectPaths(selection, n, paths_array);
}

// selects only pathItems in the selection
// search for pathItems in groupItems & compoundPaths
// *unless grouped items are compoundPaths, then ignore
// point_floor = floor limit for pathPoints length
function selectPaths(selection, point_floor, paths_array)
{ 
    for (var i = 0; i < selection.length; i++)
    {
        if (selection[i].locked || selection[i].hidden)
            continue;
        else if (selection[i].typename == "PathItem")
        {
            if((point_floor && selection[i].pathPoints.length <= point_floor)
                || selection[i].closed || selection[i].guides || selection[i].clipping)
                continue;

            paths_array.push(selection[i]);
        }
        else if (selection[i].typename == "GroupItem")
            selectPaths(selection[i].pageItems, point_floor, paths_array);
        else if (selection[i].typename == "CompoundPathItem")
            selectPaths(selection[i].pageItems, point_floor, paths_array);
        
    }
}