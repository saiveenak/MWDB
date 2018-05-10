%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                      Code to read from file and put into matrix                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                                                                                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = input('Enter The no of hash functions per layer:- ','s');                 
layer=input('Enter the number of layers- ','s');
layers=str2num(layer);
k=str2num(ks);
Matrix1 = dlmread('C:\Users\vgodava1\Desktop\MWDB\Phase3\Phase3_output\filename_d.spc',',');
length=size(Matrix1);
dims=length(2);
Min_val=zeros(1,dims-5);
Max_val=zeros(1,dims-5);
for i=1:dims-5
    Min_val(i)=min(Matrix1(:,i+5));
    Max_val(i)=max(Matrix1(:,i+5));
end
random=zeros(1,dims-5);
Mat_final=zeros(1,7);

for j=1:layers
    disp(j);
    it=1;
    bucket_vals=zeros(length(1),k);
    Mat_fin=zeros(length(1),7);
    for l=1:k
        for i=1:dims-5
            random(i)=Min_val(i)+(Max_val(i)-Min_val(i))*rand(1,1);
        end
        for i=1:length(1)
            vec=Matrix1(i,6:end);
            sum=0;
            ran_val=0;
            vec_val=0;
            for m=1:dims-5
                sum=sum+(random(m)*vec(m));
                ran_val=ran_val+(random(m)*random(m));
                vec_val=vec_val+(vec(m)*vec(m));
            end
            d=sum/(sqrt(ran_val)*sqrt(vec_val));
            if d<0.85
                bucket_vals(i,l)=0;
            else
                bucket_vals(i,l)=1;
            end
        end
    end
    for i=1:length(1)
        Mat_fin(i,1)=j;
        vec=bucket_vals(i,:);
        Mat_fin(i,2)=bi2de(vec);
        for m=1:5
            Mat_fin(i,m+2)=Matrix1(i,m);
        end
    end
    
    Mat_fin1=sortrows(Mat_fin,2);
    Mat_final=[Mat_final;Mat_fin1];
end
Mat_final=Mat_final(2:end,:);
length=size(Mat_final);

OutputPath = 'C:\Users\vgodava1\Desktop\MWDB\Phase3\Phase3_output';
OutputFile2 = 'filename_d.lsh';
OutputPath2 = fullfile(OutputPath,OutputFile2);
Output2 = fopen(OutputPath2,'w');
dlmwrite('C:\Users\vgodava1\Desktop\MWDB\Phase3\Phase3_output\filename_d.lsh',Mat_final);

