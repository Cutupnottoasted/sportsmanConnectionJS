main();
function main()
{
    if (documents.length < 1)
        return;
    
    var sel = activeDocument.selection; 

    if (!(sel instanceof Array) || sel.length < 1)
        return;
    
    var paths_array = [];
    selectPaths(sel, 1, paths_array);
    if (paths_array.length < 1)
        return;
    
    var num_anchors = 0;

    for (var i = 0; i < paths_array.length; i++)
    {
        path_points = paths_array[i].pathPoints;
        num_anchors += path_points.length; 
    }
    alert(num_anchors);
}

function selectPaths(sel, point_floor, paths_array) 
{
    for (var i = 0; i < sel.length; i++)
    {
        if (sel[i].locked || sel[i].hidden)
        {
            continue;
        }
        else if (sel[i].typename == "PathItem")
        {
            if ((point_floor && sel[i].pathPoints.length <= point_floor)
                || sel[i].guides || sel[i].clipping) 
            {
                continue; 
            }
            paths_array.push(sel[i]);
        }
        else if (sel[i].typename == "GroupItem")
        {
            selectPaths(sel[i].pageItems, point_floor, paths_array);
        }
        else if (sel[i].typename == "CompoundPathItem")
        {
            selectPaths(sel[i].pageItems, point_floor, paths_array);
        }
    }
}