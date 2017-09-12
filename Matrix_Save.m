A = full(Lp)'; % saving the graph Laplacian to .txt file
A = A + (1/N)*ones(N,N);
fileID = fopen('Lp.txt','w');
fprintf(fileID,'%f\n',A);
fclose(fileID);

B = [m;N]; % saving the graph sizes to text file
fileID = fopen('MN.txt','w');
fprintf(fileID,'%d\n',B);
fclose(fileID);

% saving the edge matrix to file
[I1, J1] = ind2sub(size(Ec),find(Ec == 1));
[I2, J2] = ind2sub(size(Ec),find(Ec == -1));
I = [I1, I2]' -1;
fileID = fopen('E.txt','w');
fprintf(fileID,'%d\n',I);
fclose(fileID);

% saving the state weight matrix
C = Qp';
fileID = fopen('Qp.txt','w');
fprintf(fileID,'%f\n',C);
fclose(fileID);

% saving the regularizers weight
fileID = fopen('gamma.txt','w');
fprintf(fileID,'%f\n',gamma);
fclose(fileID);
