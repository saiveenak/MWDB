%Path of all Scripts = C:\Users\vgodava1\Desktop\MWDB\DataR
%User Input for Path of scripts
Path = input('Enter Absolute Path of folder with videos : ','s');
Resolution = input('Enter R value to split frames : ');

Output = fopen('C:\Users\vgodava1\Desktop\MWDB\DataR\outfile2.sift','w');
VideoFiles = dir(fullfile(Path,'*mp4'));
for File = 1:size(VideoFiles)   
    Video = VideoReader(fullfile(Path,VideoFiles(File).name));
    FrameNo = 1;
    while hasFrame(Video)
        FrameInRGB = readFrame(Video);
        Frame = rgb2gray(FrameInRGB);
        [FrameHeight,FrameWidth] = size(Frame);
        
        %Splitting Resolution R into cells
        NumHeight = FrameHeight/Resolution;
        NumWidth = FrameWidth/Resolution;
        
        FrameDouble = im2double(Frame);
        [frames,descriptors,gss,dogss] = sift(FrameDouble);
        Outputmatrix = [frames;descriptors];
        TotalVectors = size(Outputmatrix);
        CellNumber = zeros(1,TotalVectors(2));
        FinalOutputMatrix =  [CellNumber;Outputmatrix];   
        
        for siftvectorsframes = 1:TotalVectors(2)
            xcord = FinalOutputMatrix(2,siftvectorsframes);
            ycord = FinalOutputMatrix(3,siftvectorsframes);
			
            
            ycell = floor(ycord/NumHeight);
            if (ycell >= R)
                ycell = ycell - 1;
            end
			
            
            xcell = floor(xcord/NumWidth);
            if (xcell >= R)
                xcell = xcell - 1;
            end
            
            cell = (ycell * R) + (xcell + 1);
            FinalOutputMatrix(1,siftvectorsframes) = cell;
        end
        Finalmatrixsize = size(FinalOutputMatrix);
        for siftvectorsframes = 1:Finalmatrixsize(2)
			fprintf(Output,'%s',int2str(File),';',int2str(FrameNo));
            SiftVectorValues = '';
            for eachvector = 1:133
                if (eachvector == 1 || eachvector == 2)
                    SiftVectorValues = strcat(SiftVectorValues,';',num2str(FinalOutputMatrix(eachvector,siftvectorsframes)));
                else
                    SiftVectorValues = strcat(SiftVectorValues,',',num2str(FinalOutputMatrix(eachvector,siftvectorsframes)));
                end     
            end
             fprintf(Output,'%s',SiftVectorValues);
             fprintf(Output,'\n');
        end       
    FrameNo = FrameNo + 1;
    end
end
fclose(Output);