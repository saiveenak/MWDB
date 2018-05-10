Matrix1 = dlmread('C:\Users\vgodava1\Desktop\MWDB\Phase3\Phase3_output\filename_d.lsh',',');
VideoPath = 'C:\Users\vgodava1\Desktop\MWDB\Phase3\InputVideo\';
VideoFiles = dir(fullfile(VideoPath,'*mp4'));
vid=input('Enter the vid num-','s');
vid1=strcat(vid,'.mp4');
for i=1:size(VideoFiles,1)
   if (strcmp((VideoFiles(i).name),vid1)) 
       break;
   end
end
vid_num=i;
frame=input('Enter the frame no- ','s');
frame_num=str2num(frame);
Coor1 = input('Enter coordinates1 in the format - X,Y:- ','s');
Coord1 = str2double(strsplit(Coor1,','));
Coor2 = input('Enter coordinates2 in the format - X,Y:- ','s');
Coord2 = str2double(strsplit(Coor2,','));
frames=input('Enter the number of frames to be output- ','s');

length=size(Matrix1);
layers=[];
Matrix2=[];
for i=1:length(1)
    if (Matrix1(i,3)==vid_num)
        if (Matrix1(i,4)==frame_num)
            x=Matrix1(i,6);
            y=Matrix1(i,7);
            if ((Coord1(1)<=x) && (x<=Coord2(1)))
                if ((Coord1(2)<=y) && (y<=Coord2(2)))
                    bucket=[Matrix1(i,1);Matrix1(i,2)];
                    layers=[layers,bucket];
                end
            end
        end
     end
end

l=Matrix1(i-1,1);
disp('no of unique sift vectors in the given region are');
unique_no=size(layers,2)/l;
disp(unique_no);
if(unique_no==0)
    disp('number of sift vectors in the given region are 0. So nothing to compare with');
    exit;
end
buckets=unique(layers','rows');
ans = [];
i=1;flag=0;start=0;bytes=1;val=0;

    for j = 1:size(Matrix1,1)
        if(val~=Matrix1(j,2))
            bytes=bytes+1;
            val=Matrix1(j,2);
        end
        if((buckets(i,1)==Matrix1(j,1))&&(buckets(i,2)==Matrix1(j,2)))
            flag=1;
            if (start==0)
                start=j;
            end
         else
            if (flag==1)
                ans=[ans;Matrix1(start:j-1,1:7)];
                flag=0;
                j=j-1;
                start=0;
                if(i<size(buckets,1))
                    i=i+1;
                else 
                    break;
                end
            end
        end
    end
    
for i=1:size(ans,1)
   if (ans(i,1)~=1) 
       break;
   end
end
count_mat=ans(1:i-1,:);
count_mat=count_mat(:,3);
[unique_rows1,~,ind1] = unique(count_mat,'rows');
unique_rows1(:,end+1) = histc(ind1,unique(ind1));
sum=0;sum1=0;
for i=1:size(unique_rows1,1)
   if(unique_rows1(i,1)==vid_num)
       continue;
   else
       sum=sum+unique_rows1(i,2);
   end
   sum1 = sum1 + unique_rows1(i,2);
end
disp('number of overall sift vectors considered are');
disp(sum);
disp('number of bytes of data from index accessed are');
disp((sum1*7*8));
Ans2 = Matrix1(:,3:4);
[unique_3,~,ind3] = unique(Ans2,'rows');
unique_3(:,end+1) = histc(ind3,unique(ind3));

Ans1 = ans(:,3:4);
[unique_rows,~,ind] = unique(Ans1,'rows');
unique_rows(:,end+1) = histc(ind,unique(ind));

j=1;unique_count=[];
for i=1:size(unique_3,1)
   if ((unique_3(i,1)==unique_rows(j,1))&&(unique_3(i,2)==unique_rows(j,2)))
       val=unique_rows(j,3)/unique_3(i,3);
       val=val*100;
       unique_count=[unique_count;val];
       j=j+1;
       continue;
   end
end
unique_rows=[unique_rows,unique_count];
sorted = sortrows(unique_rows,-4);


final_videos=[];
for k=1:str2double(frames)
    eachvideo = sorted(k,1);
    if (eachvideo==vid_num)
        k=k-1;
        continue;
    end
    Video = VideoReader(fullfile(VideoPath,VideoFiles(eachvideo).name));
    Frameno = sorted(k,2);
    A=strtok(VideoFiles(eachvideo).name,'.');
    final_videos=[final_videos,str2num(A)];
    n = 0;
    while (hasFrame(Video))
        n=n+1;
        FrameInRGB = readFrame(Video);
        Frame = rgb2gray(FrameInRGB);
        if(n==Frameno)
            figure;
            imshow(Frame);
        end
    end
end
