%{

Time is in elapsed days from time 0.
Here are the headers for each column:
1) River Mile
2) Elevation of the Top of the Top Layer (ft NAVD88) (i.e. water surface elevation)
3) Elevation of the Top of the Bottom Layer (ft NAVD88)
4) Bed Elevation (ft NAVD88)
5) Salinity of the Top Layer (ppt)
6) Salinity of the Bottom Layer (ppt) 
7) Velocity of the Top Layer (ft/sec) (positive is downstream)
8) Velocity of the Bottom Layer (ft/sec) (positive is downstream)


%}
clearvars -except data

%%-----------------READ IN OUTPUT FILE
[fileID, path] = uigetfile(".txt","Select input file");

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ["\t", " "];

% Specify column names and types
opts.VariableNames = ["Time", "River Mile", "Top Layer Elev", "Bottom Layer Elev", "Bed Elev", "Salanity Top", "Salinity Bottom", "Vel. Top", "Vel. Bottom"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Specify variable properties
opts = setvaropts(opts, "Time", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Time", "EmptyFieldRule", "auto");

% Import the data
data = readtable(fileID, opts);


%% Clear temporary variables
clear opts

%%-----------------CONVERT DATA INTO MATRIX (Location, Variable, Time)
szA = height(data);
szB = width(data);

start = 1;
dstart = start;

for i = 2:szA  
    if data{i,1} == "Time"
        break
    end
end
stop = i - 1;
dstop = stop;
intervals = szA/stop;

for i = 1:intervals

    times(i) = data{dstart,3}
    temp(start+1:stop,1:8) = data{dstart+1:dstop,2:9};
    dstart = dstart + stop;
    dstop = dstop + stop;
    dataOut(:,:,i) = temp;
end



%%--------------- PlotToe Position

timeLength = size(dataOut,3);


for i = 1:timeLength
    for j = size(dataOut,1):-1:1
        if dataOut(j,6,i) > 8.9
            temp2(i,1) = times(i);
            temp2(i,2) = dataOut(j,1);
            break
        end
    end
end

plot(temp2(:,1),temp2(:,2))
xlabel('Time elapsed (days)')
ylabel('Location of toe (River Mile)')


Rise = temp2(timeLength,2)-temp2(1,2);
Run = temp2(timeLength,1)-temp2(1,1);
MeanVelocity = Rise/Run;
predict(1,1) = temp2(timeLength,1);
predict(2,1) = temp2(timeLength,1) + 10;
predict(1,2) = temp2(timeLength,2);
predict(2,2) = temp2(timeLength,2) + 10*MeanVelocity;
hold on
plot(predict(:,1),predict(:,2))





    





