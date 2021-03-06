%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                      Code to read from file and put into matrix                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Video = input('Enter index of Video of input subsequence (V):- ','s');                
Frames = input('Enter Frame Range in the format - FirstFrame,LastFrame:- ','s');
FrameRange = str2double(strsplit(Frames,','));
k = input('Enter number of output Frame sequences required (k)','s'); 

InputDir = 'C:\Users\Joel\Desktop\ASU_Stuff\Sem1\MWD\Phase2\Phase2_output\';                  %PATH TO BE CHANGED
InputPath = fullfile(InputDir,'in_file_d.mpca');             

Input = fileread(InputPath);                 
AllVector = strsplit(Input,'\n');
AllVector = AllVector';
lenVectorCount = size(AllVector)-1;
InputMatrix = zeros(lenVectorCount(1),1);
lenEachVector = 0;
for eachvector = 1:lenVectorCount(1)
    eachVectorList = strsplit(char(AllVector(eachvector)),',');
    lenEachVector = size(eachVectorList);
    for matrix = 1:lenEachVector(2)
        value = eachVectorList{matrix};
        InputMatrix(eachvector,matrix) = str2double(value);
    end
end

InputMatrixSize = size(InputMatrix);
startV = 0;endV =0 ;
for vector = 1:InputMatrixSize(1)
    if (InputMatrix(vector,1) == str2double(Video)) && ((InputMatrix(vector,2) == FrameRange(1,1)))
        if(startV == 0)
            startV = vector;
        end    
    end
    if (InputMatrix(vector,1) == str2double(Video)) && ((InputMatrix(vector,2) == (FrameRange(1,2) + 1)))
        if(endV == 0)
            endV = vector - 1;
        end    
    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                      Code to read from file ends here                           %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix1 = InputMatrix(startV:endV,1:end);
length=size(Matrix1);
i=(FrameRange(2)-FrameRange(1))*(Matrix1(length(1),3));
Mat1=zeros(i,length(2));
frame_it=1;
frame_no=Matrix1(frame_it,2);
it=1;
cell_it=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%           Aggregating n Motion Vectors to 1 per cell in Video 1                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (Matrix1(it,2)==frame_no)
    Mat1(cell_it,1)=Matrix1(it,1);
    Mat1(cell_it,2)=Matrix1(it,2);
    Mat1(cell_it,3)=Matrix1(it,3);
    cell_num=Matrix1(it,3);
    num_vec=0;
    while (Matrix1(it,3)==cell_num)
        num_vec=num_vec+1;
        for it2=4:length(2)            
            Mat1(cell_it,it2)=Mat1(cell_it,it2)+Matrix1(it,it2);
        end
        if (it+1 > length(1))
            break;
        end
        it=it+1;
    end
    for it2=4:length(2)
        x=Mat1(cell_it,it2);
        Mat1(cell_it,it2)=x/num_vec;
    end
    if (it+1 >length(1))
            break;
    end
    cell_it=cell_it+1;
    frame_no=Matrix1(it,2);
end
Matrix1=Mat1;
TotalVideo = InputMatrix(InputMatrixSize(1),1);

ksim = zeros(TotalVideo,3);
for eachvideo = 1:TotalVideo
    startV = 0;countV = 0 ;
    for vector = 1:InputMatrixSize(1)
        if (InputMatrix(vector,1) == eachvideo)
           if (startV == 0)
                startV = vector;
           end 
           countV = countV + 1;
        end
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%           Aggregating n Motion Vectors to 1 per cell in Video 2                 %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Matrix2 = InputMatrix(startV:(startV + countV - 1),1:end);
    length2=size(Matrix2);
    Mat2=zeros(Matrix2(length(1),3)*Matrix2(length(1),3),length2(2));
    frame_it=1;
    frame_no=Matrix2(frame_it,2);
    it=1;
    cell_it=1;
    while (Matrix2(it,2)==frame_no)
        Mat2(cell_it,1)=Matrix2(it,1);
        Mat2(cell_it,2)=Matrix2(it,2);
        Mat2(cell_it,3)=Matrix2(it,3);
        cell_num=Matrix2(it,3);
        num_vec=0;
        while (Matrix2(it,3)==cell_num)
            num_vec=num_vec+1;
            for it2=4:length2(2)           
                Mat2(cell_it,it2)=Mat2(cell_it,it2)+Matrix2(it,it2);
            end
            if (it+1 > length2(1))
                break;
            end
            it=it+1;
        end
        for it2=4:length2(2)
            x=Mat2(cell_it,it2);
            Mat2(cell_it,it2)=x/num_vec;
        end
        if (it+1 >length2(1))
                break;
        end
        cell_it=cell_it+1;
        frame_no=Matrix2(it,2);
    end
      
    lenSubSequence = size(Matrix1);
    lenVideo=size(Matrix2);
    SimilarityMatrix = zeros(FrameRange(2)-FrameRange(1)+ 1,Matrix2(lenVideo(1),2)-1);
	Matrix2=Mat2;
	lenVideo=size(Matrix2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%           Generate Similarity Matrix for the subsequence and each Video         %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    start = 1;i = 1;j = 1;
    while (start < lenSubSequence(1))
        for eachcellM1 = 1:Matrix1(lenSubSequence(1),3)
            for eachcellM2 = 1:lenVideo
                if (Matrix2(eachcellM2,3) == eachcellM1)
                    vectorV1 = Matrix1(start+eachcellM1-1,4:end);
                    vectorV2 = Matrix2(eachcellM2,4:end);
                    sim=0;
                    len=size(vectorV1);
                    for it=1:len(2)
                        diff=vectorV1(it)-vectorV2(it);
                        sim=sim+abs(diff);
                    end
                    SimilarityMatrix(i,j) = SimilarityMatrix(i,j) + sim;
                    j = j + 1;
                end    
            end
            j = 1;
        end    
        start = start + Matrix1(lenSubSequence(1),3);
        i = i + 1;
    end
    append_matrix=zeros(i-1,1);
    append_matrix=append_matrix+9999;
    Sim_mat=horzcat(append_matrix,SimilarityMatrix);
    SimilarityMatrix=Sim_mat/Matrix1(lenSubSequence(1),3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%           Computation of Least Computation Subsequence for a given video        %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sizeSim = size(SimilarityMatrix);
    frameMatrix = zeros(1,sizeSim(1));
    min = 99999;runningsum = 0;
    for i = 1:sizeSim(1)
        min = 99999;
        for j = 1:sizeSim(2)
            if(any(frameMatrix == j))
               continue;
            else   
                if(i == 1) && (SimilarityMatrix(i,j) < min)
                    min = SimilarityMatrix(i,j);
                    frameref = j; 
                end
                sum = 0;
                if(i ~=1)
                    fms = size(frameMatrix);
                    s = 1;
                    for q = 1:i-1
                        sum = sum + SimilarityMatrix(s,frameMatrix(1,q));
                        s = s + 1;
                    end
                end
                sum = sum/(i-1);
                if((i ~= 1) && (abs(SimilarityMatrix(i,j) - sum) < min))
                    min = abs(SimilarityMatrix(i,j) - sum);
                    frameref = j; 
                end
            end    
        end
        runningsum = runningsum + SimilarityMatrix(i,frameref);
        frameMatrix(1,i) = frameref;
    end
    
   first = 99999;last = 1;
   lenframeMatrix = size(frameMatrix);
   for i = 1:lenframeMatrix(2)
      if(frameMatrix(1,i) > last)
          last = frameMatrix(1,i);
      end
      if(frameMatrix(1,i) < first)
          first = frameMatrix(1,i);
      end
   end 
   
   ksim(eachvideo,1) = runningsum/(FrameRange(2) - FrameRange(1));
   ksim(eachvideo,2) = first;
   ksim(eachvideo,3) = last;
   ksim(eachvideo,4) = eachvideo; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%               Code to return k most similar frame sequences                     %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ksim = sortrows(ksim,1);
disp('Most Matching Sequence in Descending order -');
for i = 1:str2double(k)
    eachvideo = ksim(i,4);
    first = ksim(i,2);
    last = ksim(i,3);
    dist = ksim(i,1);
    result = ['Video ',num2str(eachvideo),' with Frames ranging from :- ',num2str(first),' - ',num2str(last)];
    disp(result);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                          Code to visualise output                               %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VideoPath = 'C:\Users\Joel\Desktop\ASU_Stuff\Sem1\MWD\Phase2\InputVideo\';
VideoFiles = dir(fullfile(VideoPath,'*mp4'));
for i = 1:str2double(k)
    eachvideo = ksim(i,4);
    first = ksim(i,2);
    last = ksim(i,3);
    Video = VideoReader(fullfile(VideoPath,VideoFiles(eachvideo).name));
    FrameNo = 0;
    while hasFrame(Video)     
        FrameNo = FrameNo+1;
        FrameInRGB = readFrame(Video);
        Frame = rgb2gray(FrameInRGB);
        if(first <= FrameNo)&&(FrameNo <= last)
            imshow(Frame);
        end
    end
end