% [average,stdev] = average(data)
%                      
%   data        ->  multidimensional array (points x elecs x epoch)              
%
%   average     <-  matrix of average data (elecs in columns)
%   stdev       <-  matrix of standard deviation of the average
%
%   -- Note: Works only with Scan 4.1+ data files
function [a,s] = average(data)

s = std(data,0,3);
a = mean(data,3);