function [Q, R, S] = FindQRS(data, sampleRate)

    %constants Pass 1
    heightThresholdP1 = 0.25;
    
    %constants Pass 2
    heightThresholdP2 = 0.5;

    %Finds Peaks and Troughs of Dataset
    [~, TurningPointIndex] = FindTurningPoints(data);
    
    %Finds Max Value and Addumes to be the most likely QRS Complex
    [MostLikelyRValue, MostLikelyR] = max(data);
    
    % Pass 1
    minHeight = MostLikelyRValue * heightThresholdP1;

    %Minimum rr interval is 160ms.
    minDistance = 0.16 * sampleRate;

    [~, confirmedRpeaks, width] = findpeaks(data, 'MinPeakHeight', minHeight, 'MinPeakDistance', minDistance);
    
    %Pass 2
    minHeight = median(data(confirmedRpeaks)) * heightThresholdP2;
    
    [~, confirmedRpeaks] = findpeaks(data, 'MinPeakHeight', minHeight, 'MinPeakDistance', minDistance);
    
    %Iterates through present R Values to find QRS Complexes
    R = 1:length(confirmedRpeaks);
    Q = 1:length(confirmedRpeaks);
    S = 1:length(confirmedRpeaks);

    TurningPointIndexCopy = TurningPointIndex;

    for i = 1:length(confirmedRpeaks)
        Rindex = find(TurningPointIndex == confirmedRpeaks(i));
        [Q(i), R(i), S(i), TurningPointIndexCopy] = GetQRSfromPeaks(TurningPointIndexCopy, Rindex, 1, length(data));
    end

    [~, ~, Q] = find(Q);
    [~, ~, R] = find(R);
    [~, ~, S] = find(S);

end

function [Q, R, S, Array] = GetQRSfromPeaks(Array, Index, Replace, MaxVal)

    R = Array(Index);
    if Replace
        Array(Index) = 0;
    end

    if(Index ~= 1)
        Q = Array(Index - 1);
        if Replace
            Array(Index - 1) = 0;
        end
    else
        Q = 1;
    end

    if(Index ~= length(Array))
        S = Array(Index + 1);
        if Replace
            Array(Index + 1) = 0;
        end
    else
        S = MaxVal;
    end
    
end




