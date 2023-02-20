function [Pb, Te] = FindPT(data, Q, R, S, minSegmentSize)
    %Logic to segment the data for seperate analysis of data.

    noOfSegments = length(R) + 1;
    startIndex = 1;
    BME = 1;

    Pb = (1:noOfSegments) * NaN;
    Te = (1:noOfSegments) * NaN;

    for i = 1:noOfSegments
        if(i <= length(Q))
            endIndex = Q(i);
        else
            endIndex = length(data);
            BME = 3;
        end
        
        %Preprocess Segment
        X = startIndex:endIndex;
        segment = detrend(data(X));
        
        %Ignore if QRS Complex is at very end or beginning
        if (endIndex - startIndex) < minSegmentSize
            BME = 2;
            continue
        end
        
        [Pb(i), Te(i)] = DetectPT(segment, BME);
        Pb(i) = Pb(i) + startIndex - 1;
        Te(i) = Te(i) + startIndex - 1;

        BME = 2;

        if(i <= length(Q))
            startIndex = S(i);
        end

    end
end

function [Pb, Te] = DetectPT(segment, BME)

    BEthreshold = 0.45;
    distanceScorePercentage = 0.63;
    
    PP = NaN;
    TP = NaN;

    X = 1:length(segment);

    %Split segment to two parts
    halfSegment = round(length(segment)/2);
    
    Left = segment(1:halfSegment);
    Right = segment(halfSegment:end);

    %ordered in height
    [~, PossiblePWaves] = findpeaks(Right, 'SortStr', 'descend');
    PossiblePWaves = PossiblePWaves + halfSegment - 1;
    [~, PossibleTWaves] = findpeaks(Left, 'SortStr', 'descend');

    %The beginning segment can have, only p wave.
    if BME == 1
        if isempty(PossiblePWaves) || isempty(PossibleTWaves)
            joint = [PossiblePWaves, PossibleTWaves];
            PP = joint(1);
        else
            Pval = abs(segment(PossiblePWaves(1)));
            Tval = abs(segment(PossibleTWaves(1)));

            if(Tval < Pval * BEthreshold)
                PP = PossiblePWaves(1);
            end

        end
    %The ending segment can have, only t wave.
    elseif BME == 3
        if isempty(PossiblePWaves) || isempty(PossibleTWaves)
            joint = [PossiblePWaves, PossibleTWaves];
            TP = joint(1);
        else
            Pval = abs(segment(PossiblePWaves(1)));
            Tval = abs(segment(PossibleTWaves(1)));

            if(Pval < Tval * BEthreshold)
                TP = PossibleTWaves(1);
            end
        end
    %If the middle segment cannot find both twave and pwave.
    else
        if isempty(PossiblePWaves) || isempty(PossibleTWaves)
            joint = [PossiblePWaves, PossibleTWaves];
            peak = joint(1);

            segmentAfterPeak = [peak,length(segment)];
            segmentBeforePeak = [1,peak];
            
            %Use changes in gradient to find bumps in the signal to classify as either t wave or p wave.
            dydx = gradient(segment(:))./gradient(X(:));

            [TurningPointValues, TurningPointIndex] = FindTurningPoints(dydx);
            
            inBetween = find(TurningPointIndex > segmentAfterPeak(1) & TurningPointIndex < segmentAfterPeak(2)); 
            [~, maxL] = max(TurningPointValues(inBetween));

            if TurningPointIndex(inBetween(maxL)) < halfSegment
                Te = TurningPointIndex(inBetween(maxL) + 1);
            else
                Pb =  TurningPointIndex(inBetween(maxL) - 1);
            end
            
            inBetween = find(TurningPointIndex > segmentBeforePeak(1) & TurningPointIndex < segmentBeforePeak(2));
            [~, maxL] = max(TurningPointValues(inBetween));

            if TurningPointIndex(inBetween(maxL)) < halfSegment
                Te = TurningPointIndex(inBetween(maxL) + 1);
            else
                Pb =  TurningPointIndex(inBetween(maxL) - 1);
            end
            
            return

            
        end
    end

    if isnan(TP) && isnan(PP)
        pWaveScore = 1:length(PossiblePWaves);
        for i = 1:length(PossiblePWaves)
            heightScore = (segment(PossiblePWaves(i))/segment(PossiblePWaves(1))) * (1 - distanceScorePercentage);
            distanceScore = (abs(halfSegment - PossiblePWaves(i))/halfSegment) * distanceScorePercentage;
            pWaveScore(i) = heightScore + distanceScore;
        end

        tWaveScore = 1:length(PossibleTWaves);
        for i = 1:length(PossibleTWaves)
            heightScore = (segment(PossibleTWaves(i))/segment(PossibleTWaves(1))) * (1 - distanceScorePercentage);
            distanceScore = (abs(halfSegment - PossibleTWaves(i))/halfSegment) * distanceScorePercentage;
            tWaveScore(i) = heightScore + distanceScore;
        end

        [~, PP] = max(pWaveScore);
        PP = PossiblePWaves(PP);

        [~, TP] = max(tWaveScore);
        TP = PossibleTWaves(TP);
    end
    
    %Shift the P wave to the start of the p wave and shift t wave to end of t wave.
    [~, TurningPointIndex, width] = findpeaks(segment);
    halfSegment = round((PP + TP)/2);
    
    %If a p wave cannot be found in the ending segment
    if isnan(PP)
        Pb = NaN;
    else
        Pi = find(TurningPointIndex == PP);
        Pb = ceil(PP - width(Pi));

        if Pb < 1
            Pb = 1;
        elseif Pb < halfSegment && ~isnan(TP)
            Pb = halfSegment;
        end
    end

    %If a t wave cannot be found in the beginning segment
    if isnan(TP)
        Te = NaN;
    else
        Ti = find(TurningPointIndex == TP);
        Te = floor(TP + width(Ti));

        if Te > length(segment)
            Te = length(segment);
        elseif Te > halfSegment && ~isnan(PP)
            Te = halfSegment;
        end
    end

    %Move away if on top of eachother
    if Pb == Te
        Pb = Pb + ceil((PP - Pb) * 0.10);
        Te = Te - floor((Te - TP) * 0.10);
    end

end
