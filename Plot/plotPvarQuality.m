function plotPvarQuality(filename)
    [data8bit,~,~]=load8BitBinary(filename,'dontsave');
    
    % Fetch date
    currentDirectory = pwd;
    [~,foldername] = fileparts(currentDirectory);
    
    % Plot PointwiseVariance
    plotPointwiseVariance(data8bit(1:800,:))
    xlabel('Samples@5GS/s')
    ylabel('PointwiseVariance')
    title([foldername(1:10) ':' filename])

    % Saving PointwiseVariance
    outputFilename = ['Pvar-' filename '.jpg'];
    outputFiletype = '-djpeg';
    print(outputFilename,outputFiletype);
    
    % Plot StackedWaveforms
    plotStackedWaveforms(data8bit(1:800,:))
    xlabel('Samples@5GS/s')
    ylabel('StackedWaveforms')
    title([foldername(1:10) ':' filename])
    
    % Saving StackedWaveforms
    outputFilename = ['Stacked-' filename '.jpg'];
    outputFiletype = '-djpeg';
    print(outputFilename,outputFiletype);
end