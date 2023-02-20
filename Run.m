clear;
clc;
warning('off');
close;

sampleRate = input('Enter the sample rate of the data: ');
fileName = input('Enter the header of the files to be analysed: ', 's');
field = input('Enter the field name of the data: ', 's');

data = {};

loadedData = load([fileName]);
data{1} = loadedData.(field);


disp('-----------------------------------')
disp('ECG Analyser')
disp('-----------------------------------')

samplesToAnalyse = [1];

for index = 1:length(samplesToAnalyse)
    sampleNo = samplesToAnalyse(index);

    %PreProcessing
    dataToAnalyse = FilterData(data{sampleNo}, 3, sampleRate);
    X = (1:length(dataToAnalyse))/sampleRate;

    %Plot data
    Figure = figure('name', ['Sample Number ', num2str(sampleNo)]);
    hold on
    plot(X, dataToAnalyse, 'k');

    %Find QRS Complexes
    [Q, R, S] = FindQRS(dataToAnalyse, sampleRate);

    %Plot QRS Complexes
    scatter(Q/sampleRate, dataToAnalyse(Q), 'rx');
    scatter(R/sampleRate, dataToAnalyse(R), 'gv');
    scatter(S/sampleRate, dataToAnalyse(S), 'b+');

    %Find QRS Interval
    QRSInterval = mean(S - Q)/sampleRate;

    %Find P wave beginning and T wave Ends
    [Pb, Te] = FindPT(dataToAnalyse, Q, R, S, QRSInterval);

    %Plot P and T waves
    PbNoNan = Pb(~isnan(Pb));
    TeNoNan = Te(~isnan(Te));
    scatter(PbNoNan/sampleRate, dataToAnalyse(PbNoNan), 'ms');
    scatter(TeNoNan/sampleRate, dataToAnalyse(TeNoNan), 'cd');

    %Find PR and QT Intervals
    PRIntervals = 1:length(Pb);
    QTIntervals = 1:length(Te);

    for i = 1:length(Pb)
        %Check if no P wave is present
        if ~isnan(Pb(i))
            %Check to see if the PR Interval can be calculated
            if Pb(i) < S(end)
                PRIntervals(i) = Q(i) - Pb(i);
            else
                PRIntervals(i) = NaN;
            end
        end

        %Check if no P wave is present
        if ~isnan(Te(i))
            %Check to see if the PR Interval can be calculated
            if Te(i) > Q(1)
                QTIntervals(i) = Te(i) - Q(i - 1);
            else
                QTIntervals(i) = NaN;
            end
        end
    end

    PRInterval = mean(PRIntervals/sampleRate, 'omitnan');
    QTInterval = mean(QTIntervals/sampleRate, 'omitnan');

    

    %Check if results are normal
    if PRInterval > 0.12 && PRInterval < 0.2
        PRNormal = 'Normal';
    else
        PRNormal = 'Abnormal';
    end

    if QTInterval < 0.38
        QTNormal = 'Normal';
    else
        QTNormal = 'Abnormal';
    end

    if QRSInterval < 0.1
        QRSNormal = 'Normal';
    else
        QRSNormal = 'Abnormal';
    end
        

    %Display Results
    disp(['Sample Number ', num2str(sampleNo)])
    fprintf('PR Interval: %.2fs (%s)\n', PRInterval, PRNormal);
    fprintf('QT Interval: %.2fs (%s)\n', QTInterval, QTNormal);
    fprintf('QRS Duration: %.2fs (%s)\n', QRSInterval, QRSNormal);
    disp('-----------------------------------');
    
    legend('Filtered Signal', 'Q', 'R', 'S', 'P', 'T')
    hold off;
end