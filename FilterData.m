function FilteredData = FilterData(data, order, sampleRate)
    Fn = sampleRate/2;

    cutoff = [0.5,10];
  
    [b , a] = butter(order, cutoff/Fn, 'bandpass');
    FilteredData = filtfilt(b, a, data); 
end