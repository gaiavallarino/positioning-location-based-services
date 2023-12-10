clear all
close all
clc
%import RTK points
trajectory = importdata( 'RTK.csv' );
%import smartphone points
smartPhone = importdata( 'smartphone.csv' );


%plot trajectories

%cycle on smartphones points
for i=1:numpoints
	%find the nearest RTK points
	for j=1:numRTK
	
	end
	%fill vector of errors
end

%plot and calculate statistics of the errors

