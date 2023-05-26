
main();

function main()
{
    if (documents.length < 1)
        return; 

    var sel = activeDocument.selection;
    if (!(sel instanceof Array) || sel.length < 2) return;

    var target_path = sel[0];
    var target_point = []; 

    var i;

    if (target_path.typename == "PathItem")
    {
        var path_points = target_path.pathPoints;

        for (i = 0; i < path_points.length; i++)
        {
            if (isAnchor(path_points[i])) 
            {
                if (target_path.length < 1)
                {
                    target_point = path_points[i].anchor;
                }
                else // more than 1 path point
                {
                    target_point = [];
                    break; 
                }
            }
        }
    }
    // if two or more points are selected find the center
    if (target_point.length < 1)
    {
        var vb = target_path.visibleBounds;
        target_point = [(vb[0] + vb[2]) / 2, (vb[1] + vb[3]) / 2];
    }

    var paths = [];
    selectPaths(sel.slice(1), 0, paths);
    var j;

    for (i = 0; i < paths.length; i++)
    {
        path_points = paths[i].pathPoints;
        
        for (var j = 0; j < path_points.length; j++)
        {
            if (isAnchor(path_points[j]))
            {
                target_path.duplicate().translate(path_points[j].anchor[0] - target_point[0],
                                                    path_points.anchor[1] - target_point[1]);
            }
        }
    }
}

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

