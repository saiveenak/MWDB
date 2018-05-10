%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                      Code to read from file and put into matrix                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                                                                                 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = input('Enter Target Dimensionality:- ','s');                 

InputDir = 'C:\Users\Joel\Desktop\ASU_Stuff\Sem1\MWD\Phase2\Phase1_output\';                    
InputPath = fullfile(InputDir,'in_file.mvect');             

Input = fileread(InputPath);                 
AllVector = strsplit(Input,'\n');
AllVector = AllVector';
lenVectorCount = size(AllVector)-1;
Matrix1 = zeros(lenVectorCount(1),1);
lenEachVector = 0;
for eachvector = 1:lenVectorCount(1)
    eachVectorList = strsplit(char(AllVector(eachvector)),',');
    lenEachVector = size(eachVectorList);
    for matrix = 1:lenEachVector(2)
        value = eachVectorList{matrix};
        Matrix1(eachvector,matrix) = str2double(value);
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%                      Code to read from file ends here                           %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutputPath = 'C:\Users\Joel\Desktop\ASU_Stuff\Sem1\MWD\Phase2\Phase2_output\';
OutputFile1 = 'in_file_d.mpca';                                               
OutputPath1 = fullfile(OutputPath,OutputFile1);

OutputFile2 = 'in_file_d2.mpca';
OutputPath2 = fullfile(OutputPath,OutputFile2);
Output2 = fopen(OutputPath2,'w');

InfoMatrix = Matrix1(1:end,1:3);
VectorMatrix = Matrix1(1:end,4:end);

CovarianceMatrix = cov(VectorMatrix);
[BasisMatrix,DiagonalMatrix] = eig(CovarianceMatrix);
[DiagonalMatrix,index] = sort(diag(DiagonalMatrix),'descend');
DiagonalMatrix = diag(DiagonalMatrix);
BasisMatrix = BasisMatrix(:,index);

BasisMatrix = BasisMatrix(1:end,1:str2double(dim));
PCAMatrixFull = VectorMatrix * BasisMatrix;
PCAMatrix = PCAMatrixFull(1:end,1:str2double(dim));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                     Code for <original_index,score>                                         %%%%%%%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BasisMatrixSize = size(BasisMatrix);
IndexMatrix = zeros(BasisMatrixSize(1),1);
for i = 1:BasisMatrixSize(1)
    IndexMatrix(i,1) = i;
end    

ScoreMatrix = horzcat(IndexMatrix,BasisMatrix);

for j = 2:(BasisMatrixSize(2)+1)
    AscPrintMatrix = sortrows(ScoreMatrix,j);
    k = BasisMatrixSize(1);
    indexVal = '';index = 0;score = 0;
    while(k > 0)
        score = AscPrintMatrix(k,j);
        index = AscPrintMatrix(k,1);
        indexVal = strcat(indexVal,'<',num2str(index),',',num2str(score),'>/' );
        k = k - 1;
    end
    fprintf(Output2,'%s',indexVal);
    fprintf(Output2,'\n'); 
end    
fclose(Output2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%                                     Code for <original_index,score> ends here                               %%%%%%%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MvectPCAMatrix = horzcat(InfoMatrix,PCAMatrix);
mvectsize = size(MvectPCAMatrix);
Output1 = fopen(OutputPath1,'w');
for printlen = 1:mvectsize(1)
    MvectVectorValues = '';
    for eachvector = 1:(str2double(dim)+3)
        if (isempty(MvectVectorValues))
            MvectVectorValues = strcat(num2str(MvectPCAMatrix(printlen,eachvector)));
        else
            MvectVectorValues = strcat(MvectVectorValues,',',num2str(MvectPCAMatrix(printlen,eachvector)));
        end
    end    
    fprintf(Output1,'%s',MvectVectorValues);
    fprintf(Output1,'\n');
end 
fclose(Output1);