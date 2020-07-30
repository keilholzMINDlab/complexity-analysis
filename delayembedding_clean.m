%%%reads in a time course and creates delay embedding
%%max embedding dimension is 5 by default
%%updated 3/5/20 to clean up for paper

%%%read in preprocessed zscored hcp data from behnaz or have it preloaded
%close all; clear all;
%dat=load ('W:/keilholz-lab/Behnaz/myMatlab/QPP_HCP817vSep18/BR_GWCR_0318.mat')

%%%for reading other files
% fid=fopen(fname,'r');
% formatSpec = '%f';
% ts=fscanf(fid, formatSpec);
% fclose(fid);
%fprintf('Read in data OK! \n');
%ts1=reshape(ts,4800,200);

cdimmat=zeros(4,817,360);
elagmat=zeros(4,817,360);
edimmat=zeros(4,817,360);
le=zeros(4,817,360);
apen=zeros(4,817,360);
for scannum=1:1
    for ind=1:817
        ts=B(ind,scannum);
        ts2=ts{1,:};
        for p=1:360
            ts1=ts2(p,:)';

            %%each column is treated as a time course
           [psrecon, elag, edim]=phaseSpaceReconstruction(ts1);
            %phaseSpaceReconstruction(ts1, elag, edim);
            corDim = correlationDimension(ts1,elag, edim);
            cdimmat(scannum,ind,p)=corDim;
            elagmat(scannum,ind,p)=elag;
            edimmat(scannum,ind,p)=edim;
            lyapExp = lyapunovExponent(ts1,1.39,elagmat(scannum,ind,p),edimmat(scannum,ind,p));
            le(scannum,ind,p)=lyapExp;
            apen(scannum,ind,p)=approximateEntropy(ts1,elagmat(scannum,ind,p));
        end
    end
end

save('delayembeddingclean_nogsr.mat', 'le', 'apen', 'edimmat', 'cdimmat', 'elagmat', '-v7.3');

%script to look at effects of filtering

% for i=1:100
%     test=randn(1,1200);
%     testa(i,:)=bandpass(test, [0.01 0.1], 1.4);
%     testb(i,:)=bandpass(test, [0.1 0.2], 1.4);
%     testc(i,:)=bandpass(test, [0.01 0.2], 1.4);
%     cdtesta(i)=correlationDimension(testa(i,:));
%     cdtestb(i)=correlationDimension(testb(i,:));
%     cdtestc(i)=correlationDimension(testc(i,:));
%     apa(i,:)=pwelch(testa(i,:),[],[],[],1/0.72);
%     apb(i,:)=pwelch(testb(i,:),[],[],[],1/0.72);
%     apc(i,:)=pwelch(testc(i,:),[],[],[],1/0.72);
%     letest(i)=lyapunovExponent(testa(i,:));
%     letestb(i)=lyapunovExponent(testb(i,:));
%     letestc(i)=lyapunovExponent(testc(i,:));
%     aptesta(i)=approximateEntropy(testa(i,:));
%     aptestb(i)=approximateEntropy(testb(i,:));
%     aptestc(i)=approximateEntropy(testc(i,:));
% 
% end

% effects of noise
for i=1:100
    test=randn(1,1200);
    testa(i,:)=bandpass(test, [0.01 0.1], 1.4);
    testb(i,:)= testa(i,:)+ 0.2*randn(1,1200);
    testc(i,:)=testa(i,:)+0.5*randn(1,1200);
    cdtesta(i)=correlationDimension(testa(i,:));
    cdtestb(i)=correlationDimension(testb(i,:));
    cdtestc(i)=correlationDimension(testc(i,:));
    apa(i,:)=pwelch(testa(i,:),[],[],[],1/0.72);
    apb(i,:)=pwelch(testb(i,:),[],[],[],1/0.72);
    apc(i,:)=pwelch(testc(i,:),[],[],[],1/0.72);
    letest(i)=lyapunovExponent(testa(i,:));
    letestb(i)=lyapunovExponent(testb(i,:));
    letestc(i)=lyapunovExponent(testc(i,:));
    aptesta(i)=approximateEntropy(testa(i,:));
    aptestb(i)=approximateEntropy(testb(i,:));
    aptestc(i)=approximateEntropy(testc(i,:));

end


