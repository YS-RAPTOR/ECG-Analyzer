sampleRate = 360;

filesHeader = 'Sample_';
fileMaxNo = 8;
data = {};

close;

for i = 1:fileMaxNo
    loadedData = load(['Data/', filesHeader, num2str(i), '.mat']);
    variableName = fieldnames(loadedData);
    data{i} = eval(['loadedData.', variableName{1}]);
end

%Find QRS Complexesy
for sampleNo = 3
    close;
    %PreProcessing
    data{sampleNo} = FilterData(data{sampleNo}, 3, sampleRate);

    %Find QRS Complexes
    [Q, R, S] = FindQRS(data{sampleNo}, sampleRate);

    %Data for Part A and B
    RRInterval = mean(R(2:end) - R(1:(end - 1)));
    QRSInterval = mean(S - Q);
    

    [Pb, Te] = FindPT(data{sampleNo}, Q, R, S, QRSInterval);

    %Calculations for Part A and B
    RRInterval = RRInterval/sampleRate;
    bpm = 60/RRInterval;
    QRSInterval = QRSInterval/sampleRate;
end

function output = AnalyseData(data, name, sampleRate)
    noOfSamples = length(data);
    time = noOfSamples/sampleRate;

    figure('name',name);
    subplot(2, 1, 1)
    plot(data);
    
    subplot(2, 1, 2)
    Y = fft(data);
    Y = Y(1:noOfSamples/2);
    Magnitude = real(Y);

    X = 0:(noOfSamples/2 - 1);
    X = X/time;

    plot(X, Magnitude, 'r');
    ylabel('Magnitude');
    xlabel('Frequency');

    output = NaN;

end
