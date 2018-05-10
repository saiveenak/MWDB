Matrix1 = dlmread('C:\Users\vgodava1\Desktop\MWDB\Phase3\Phase3_output\filename_d_k.gspc',',');
MatrixSize = size(Matrix1);
m = input(' Enter the top M frames to be visualized: ');
s1 = input(' Enter the seed video,frame value: ','s');
s2 = input(' Enter the seed video,frame value: ','s');
s3 = input(' Enter the seed video,frame value: ','s');

%%%%%%%% Calculating  K value %%%%%%%%%
vcount = 1;
fcount = 1;
v = 1;
f = 1;
k = 0;
while ( vcount == 1 && fcount == 1 && f < MatrixSize(1) )
    if ( Matrix1(v,1) == vcount && Matrix1(f + 1,2) == fcount )
        k = k + 1;
    end  
    v = v + 1;
    f = f + 1;
end 
%%%%%%%% K value computation ends %%%%%%%%%%%%%%%%
%%%%%%%% Video Matrix acual representation %%%%%%%%%%
VideoMatrix = zeros(63,2);
for iter = 1:63
    VideoMatrix(iter,1) = iter;
end
%VideoMatrix(1:totalvideos,2) = [1,2,3,4,5];
%VideoMatrix(1:totalvideos,2) = [1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,28,29,3,30,31,32,33,34,35,36,37,38,39,4,40,41,42,43,44,45,46,47,48,49,5,50,51,52,53,54,55,56,57,58,59,6,60,61,62,63];
AA = [1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,24,25,26,27,28,29,3,30,31,32,33,34,35,36,37,38,39,4,40,41,42,43,44,45,46,47,48,49,5,50,51,52,53,54,55,56,57,58,59,6,60,61,62,63,7,8,9];
AA = AA';
VideoMatrix(1:63,2) = AA;
%%%%%%%% Video Matrix representation ends %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Videos and Frames counts %%%%%%%%%
totalvideos = Matrix1(f,1);
CountsMatrix = zeros(totalvideos,2);
framecount = 1;
videocount = 1;
i = 1;
while ( videocount <= totalvideos )
    while ( Matrix1(i,1) == videocount && i < MatrixSize(1) )
        i = i + 1;
    end
    if ( videocount < totalvideos )     
        framecount = ceil((i-1)/k);
    else 
        framecount = ceil(i/k);
    end    
    CountsMatrix(videocount,1) = videocount;
    if (videocount == 1)
        CountsMatrix(videocount,2) = framecount;
    else
        CountsMatrix(videocount,2) = (framecount - sum(CountsMatrix(1:(videocount-1),2)));
    end
    videocount = videocount + 1;
end 
%%%%%%%%%%% Video and Frame count ends %%%%%%%%
 
%%%%%%%%%%% Ordering of video,frames %%%%%%%%%%% 
fframecount = sum(CountsMatrix(1:end,2));
TempMatrix = zeros(framecount,3);
index1 = 0;
a = 1;
iter = 1;
val = 1;
while ( iter <= MatrixSize(1) )
        if ( Matrix1(iter,1) == 1)
            index1 = Matrix1(iter,1)*Matrix1(iter,2);
        else
            index1 = 0;
            a = 1;
            while ( a < Matrix1(iter,1) )
                index1 = index1 + CountsMatrix(a,2);
                a = a + 1;
            end
            index1 = index1 + Matrix1(iter,2);
        end    
    TempMatrix(val,3) = index1;
    TempMatrix(val,1) = a;
    TempMatrix(val,2) = Matrix1(iter,2);
    val = val + 1;
    iter = iter + k;
end
FTempMatrix = zeros(fframecount,3);
FTempMatrix(1,1:end) = TempMatrix(1,1:end);
fiter = 2;
for titer = 2:framecount
    if ( TempMatrix(titer,3) ~= TempMatrix(titer-1,3) )
        FTempMatrix(fiter,1:end) = TempMatrix(titer,1:end);
        fiter = fiter + 1;
    end    
end  
FTempMatrix(all(FTempMatrix==0,2),:)=[];
%%%%%%%% End of ordering of video,frames  %%%%%%%%%%

%%%%%%%% Similarity Matrix creation %%%%%%%%%%%%
newsize = size(FTempMatrix);
fframecount = newsize(1);
SimilarityMatrix = zeros(fframecount,fframecount);
for fiter = 1:MatrixSize(1)
    for tcount = 1:fframecount
        if ( Matrix1(fiter,1) == FTempMatrix(tcount,1) && Matrix1(fiter,2) == FTempMatrix(tcount,2) )
            for kcount = 1:k
                for ttcount = 1:fframecount
                    if ( Matrix1(fiter,3) == FTempMatrix(ttcount,1) && Matrix1(fiter,4) == FTempMatrix(ttcount,2) )
                        SimilarityMatrix(tcount,ttcount) = Matrix1(fiter,5);
                    end    
                end    
            end    
        end    
    end    
end
%%%%%%%% End of Similarity Matrix creation %%%%%%%%%

%%%%%%% Seed Nodes matrix creation %%%%%%%% 
s1 = strsplit(s1,',');
s2 = strsplit(s2,',');
s3 = strsplit(s3,',');
s1 = str2double(s1);
s2 = str2double(s2);
s3 = str2double(s3);
IndexMatrix = zeros(3,1);
sm = [s1;s2;s3];
for seed = 1:3
    for iter = 1:fframecount
        if ( FTempMatrix(iter,1) == sm(seed,1) && FTempMatrix(iter,2) == sm(seed,2) )
            IndexMatrix(seed,1) = (FTempMatrix(iter,3));
        end    
    end    
end   
%%%%%%% Seed Nodes matrix %%%%%%%%%%%%

%%%%%%%% Page Rank Computation %%%%%%%%%%%%%%%%%%%%%
NewSimMatrix = SimilarityMatrix';
csum = 0;
for citer = 1:fframecount
    csum = 0;
    for riter = 1:fframecount
        csum = csum + NewSimMatrix(riter,citer);
    end 
    NewSimMatrix(1:fframecount,citer) = NewSimMatrix(1:fframecount,citer)/csum;
end

%%%%%%%% RPR-2 Formula Implementation %%%%%%%%%%%%%%%%%%%%%%
TeleportationMatrix = NewSimMatrix;
b = input('Enter the beta value(between 0 and 1): ');
for iter = 1:seed
    SM = zeros(fframecount,1);
    SM( (IndexMatrix(iter)),1 ) = 1;
    Prod1 = (ones(fframecount,fframecount) - ((1-b)*TeleportationMatrix));
    InvProd1 = pinv(Prod1);
    Product2 = b*SM;
    if (iter ==1)
        S1 = InvProd1*Product2; %%% PIE MATRIX %%%%
    end
    if (iter ==2)
        S2 = InvProd1*Product2;
    end
    if (iter ==3)
        S3 = InvProd1*Product2;
    end    
end
%%%%% RPR-1 logic %%%%%%%%%
P1 = sum(S1((1:fframecount),1));
P2 = sum(S2((1:fframecount),1));
P3 = sum(S3((1:fframecount),1));
PowerP = [P1,P2,P3];
seedcount = zeros(3,2);
Scrit = max(PowerP);
counts = 0;
for newiter = 1:seed
    if( abs((Scrit-PowerP(newiter))) <= 0.00000001 )
        seedcount(newiter,2) = 1; 
        seedcount(newiter,1) = newiter;
        counts = counts + 1;
    end  
end
values = 0;
if (seedcount(1,2) == 1)
    values = values + S1;
end 
if (seedcount(2,2) == 1)
    values = values + S2;
end 
if (seedcount(3,2) == 1)
    values = values + S3;
end 
if( counts == 1)
    finval = values;
else
    finval = (values/3);
end    
RankMatrix = zeros(fframecount,3);
RankMatrix(1:fframecount,1) = FTempMatrix(1:fframecount,1);
RankMatrix(1:fframecount,2) = FTempMatrix(1:fframecount,2);
RankMatrix(1:fframecount,3) = finval;   
NewRank = sortrows(RankMatrix,3);
NewRank = flipud(NewRank);
for mj = 1:m
    disp(NewRank(mj,1:2));
end

%%%%%%%% Code to display video,frames %%%%%%%%%
VideoPath = 'C:\Users\vgodava1\Desktop\MWDB\Phase3\InputVideo\';
VideoFiles = dir(fullfile(VideoPath,'*mp4'));
for i = 1:m
    eachvideo = NewRank(i,1);
    first = NewRank(i,2);
    Video = VideoReader(fullfile(VideoPath,VideoFiles(eachvideo).name));
    FrameNo = 1;
    while hasFrame(Video)     
        FrameNo = FrameNo+1;
        FrameInRGB = readFrame(Video);
        Frame = rgb2gray(FrameInRGB);
        if(first == (FrameNo-1) )
            %disp(NewRank(i,1:2));
            imshow(Frame);
        end
    end
end