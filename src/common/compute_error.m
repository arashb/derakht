function [ err ] = compute_error( c, canalytical, error_type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch error_type
    case 0
        norm = max(max(max(canalytical)));
        diff = c - canalytical;
        err = max(max(max(abs(diff))))/norm;
    case 1
        norm = sqrt(sum(sum(sum(canalytical.*canalytical))));
        diff = c - canalytical;
        err = sqrt(sum(sum(sum(diff.*diff))))/norm;
end
end
