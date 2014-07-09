%/* ************************************************** */
function [val, funcval] = tree_do_refine(qtree, func, maxErrPerNode, maxLevel, resPerNode, t)
[err, funcval] = tree_compute_error(qtree, func, resPerNode,t);

if qtree.level == maxLevel || err <=  maxErrPerNode 
val = false; 
return;
end

val = true;
end