%Path of all Scripts = C:\Users\vgodava1\Desktop\MWDB\DataR
%User Input for Path of scripts
Path = input('Enter Absolute Path of folder with videos : ','s');
Resolution = input('Enter R value to split frames : ');
Bins = input('Enter the value of bins(n) : '); 
Output = fopen('C:\Users\vgodava1\Desktop\MWDB\DataR\outputfile1.chst','w');
VideoFiles = dir(fullfile(Path,'*mp4'));
for File = 1:size(VideoFiles)   
    %disp (VideoFiles(File).name)
    Video = VideoReader(fullfile(Path,VideoFiles(File).name));
    FrameNo = 1;
    while hasFrame(Video)
        FrameInRGB = readFrame(Video);
        %Converting Frame to Grayscale
        Frame = rgb2gray(FrameInRGB);
        %whos Frame
        [FrameHeight,FrameWidth] = size(Frame);
        HeightSplits = FrameHeight/Resolution;
        FrameHeightSplits = zeros(Resolution+1,1); % to initialize array 
        for counter = 2:Resolution + 1  
            FrameHeightSplits(counter) = ((counter - 1) * HeightSplits);      
        end
        WidthSplits = FrameWidth/Resolution;
        FrameWidthSplits = zeros(Resolution+1,1); 
        for counter = 2:Resolution + 1  
            FrameWidthSplits(counter) = ((counter - 1) * WidthSplits);      
        end
        SmallFrameCount = 1;
        for NewHeight = 1:Resolution
            for NewWidth = 1:Resolution
                SmallFrame = Frame(FrameHeightSplits(NewHeight)+1:FrameHeightSplits(NewHeight+1), FrameWidthSplits(NewWidth)+1:FrameWidthSplits(NewWidth+1));
                [PixelFreq,Range] = imhist(SmallFrame,Bins);
                PixelFreqOutput = '';
                for Value = 1:size(PixelFreq)
                    if isempty(PixelFreqOutput)
                        PixelFreqOutput = strcat(PixelFreqOutput,';',int2str(PixelFreq(Value)));
                    else
                        PixelFreqOutput = strcat(PixelFreqOutput,',',int2str(PixelFreq(Value)));
                    end    
                end    
                fprintf(Output,'%s',int2str(File),';',int2str(FrameNo),';',int2str(SmallFrameCount));
                fprintf(Output,'%s',PixelFreqOutput);
                fprintf(Output,'\n');
                SmallFrameCount = SmallFrameCount + 1;
            end 
        end       
    FrameNo = FrameNo + 1;
    end
end
fclose(Output);







