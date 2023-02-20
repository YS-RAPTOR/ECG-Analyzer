function [TurningPointValues, TurningPointIndex, PeakValues, PeakIndices, TroughValues, TroughIndices] = FindTurningPoints(data)
    [PeakValues, PeakIndices] = findpeaks(data);
    [TroughValues, TroughIndices] = findpeaks(-data);
    
    TurningPointValueIndexTable = [PeakValues, PeakIndices; -TroughValues, TroughIndices];
    TurningPointValueIndexTable = sortrows(TurningPointValueIndexTable, 2);

    TurningPointValues = TurningPointValueIndexTable(:,1);
    TurningPointIndex = TurningPointValueIndexTable(:,2);
    
end