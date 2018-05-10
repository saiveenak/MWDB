%Path of all Scripts = C:\Users\vgodava1\Desktop\MWDB\DataR
%User Input for Path of scripts
Path = input('Enter Absloute Path of folder with videos : ','s');
R = input('Enter R value to split frames : ');

%File Handler for Writing Output to File
Output = fopen('C:\Users\vgodava1\Desktop\MWDB\DataR\output3.mvect','w');



VideoFiles = dir(fullfile(Path,'*mp4'));
for File = 1:size(VideoFiles)
    VideoName = fullfile(Path,VideoFiles(File).name);
    Video = VideoReader(VideoName);
    FrameInRGB = readFrame(Video);
    Frame = rgb2gray(FrameInRGB);
    [FrameHeight,FrameWidth] = size(Frame);
    
	filepath = fullfile(Path,VideoFiles(File).name);
    ffilepath = strcat('C:\Users\vgodava1\Desktop\MWDB\DataR\Siftapp\Outputfolder\Siftapp.exe',{' '},filepath,'>','C:\Users\vgodava1\Desktop\MWDB\DataR\motion.txt');
    system(char(ffilepath));
    answer = fileread('C:\Users\vgodava1\Desktop\MWDB\DataR\motion.txt');
    eachvector = strsplit(answer,';');
    eachvector = eachvector';
    length = size(eachvector) - 1;
    FinalMotionV = zeros(length(1),8);
    for vector = 1:length(1)
        abc = strsplit(char(eachvector(vector)),',');
        length2 = size(abc);
        for matrix = 1:length2(2)
           value = abc{matrix};
           FinalMotionV(vector,matrix) = str2double(value);
        end
    end
    FinalMotionV = FinalMotionV';
    NumHeight = FrameHeight/R;
    NumSplits = FrameWidth/R;    
    
    Cellloc = zeros(1,length(1));
    FinalOutputMatrix =  [Cellloc;FinalMotionV];   
    TotalVectors = size(FinalOutputMatrix);
	
    for motionvectorframes = 1:length(1)
        xcord = FinalOutputMatrix(8,motionvectorframes);
        ycord = FinalOutputMatrix(9,motionvectorframes);
		
        ycelly = floor(ycord/NumHeight);
        if (ycell >= R)
            ycell = ycell - 1;
        end
        xcell = floor(xcord/NumWidth);
        if (xcell >= R)
            xcell = xcell - 1;
        end
        cell = (ycell * R) + (xcell + 1);
		
        FinalOutputMatrix(1,motionvectorframes) = cell;
    end
    
    FinalOutputMatrix([1 2],:) = FinalOutputMatrix([2 1 ],:);
    FinalOutputMatrix = FinalOutputMatrix';
    FinalOutputMatrix = sortrows(FinalOutputMatrix,[1 2]);
    FinalOutputMatrix = FinalOutputMatrix';
    for vectors = 1:length(1)
        fprintf(Output,'%s',int2str(File));
            MotionVectorValues = '';
            for eachvector = 1:9
                if (eachvector == 1 || eachvector == 2)
                    MotionVectorValues = strcat(MotionVectorValues,';',num2str(FinalOutputMatrix(eachvector,vectors)));
                else
                    MotionVectorValues = strcat(MotionVectorValues,',',num2str(FinalOutputMatrix(eachvector,vectors)));
                end     
            end
             fprintf(Output,'%s',MotionVectorValues);
             fprintf(Output,'\n');         
    end  
end
fclose(Output);