%calculates center frequency within low freq range along with cc and
%complexity measures
%cleaned up 3.4.20 to check figures in paper and document methods

%must first load 'W:/keilholz-lab/Behnaz/myMatlab/QPP_HCP817vSep18/BR_GWCR_0318.mat'
%817 individuals, gm, wm, csf regressed, glasser parcellation

for scannum=1:1 %four scans, 1 and 2 on day 1, 3 and 4 on day 2

    for ind=1:817  %loop over subjects
        ts=B(ind,scannum);
        ts2=ts{1,:};
        for p=1:360 %loop over parcels
            ts1=ts2(p,:)';
            [pts,f] = pwelch(ts1,[],[],[],1/0.72); %tr=0.72
            ptsall(scannum,ind,p,:)=pts;
        end
    end
end
fprintf('Calculated psd OK! \n');

%calculates amplitude time frequency per roi per subject 
area=zeros(4,817,360);
for scannum=1:1
    for ind=1:817   
        for p=1:360
            for k=1:size(f)
                area(scannum,ind,p)=ptsall(scannum,ind,p,k).*f(k)+area(scannum,ind,p);
            end

        end
    end
end
fprintf('Calculated area OK! \n');

for scannum=1:1
    for ind=1:817 %finds weighted average freq. is discretized.
        for p=1:360
            c=0;
            for k=1:size(f)
                c=ptsall(scannum, ind,p,k).*f(k)+c;
                if (c>(area(scannum,ind,p)/2))
                    avgfreq(scannum,ind,p)=f(k);
                    break           
                end
            end
        end
    end
end
fprintf('Calculated center freq OK! \n');

%calculate correlation between parcels
ccmat=zeros(4,360,360,817);     
for scannum=1:1
    for ind=1:817
        ts=B(ind,scannum);
        ts2=ts{1,:};
        for p1=1:360
            for p2=p1:360
                ts1=ts2(p1,:);
                ts3=ts2(p2,:);
                cctemp=corrcoef(ts1,ts3);
                ccmat(scannum,p1,p2,ind)=cctemp(2,1);
            end
        end
    end
end
fprintf('Calculated correlation matrix OK! \n');

%calculate average center frequency for two parcels to compare to cc
avgfreqmat=zeros(4,360,360,817);   
difffreqmat=zeros(4,360,360,817);
for scannum=1:1
    for ind=1:817
        for p1=1:360
            for p2=p1:360
                avgfreqmat(scannum,p1,p2,ind)=(avgfreq(scannum,ind,p1)+avgfreq(scannum,ind,p2))/2;
                difffreqmat(scannum, p1,p2,ind)=avgfreq(scannum,ind,p1)-avgfreq(scannum,ind,p2);
            end
        end
    end
end
fprintf('Calculated avg center freq and difference OK!\n');

save('freqccall_nogsr.mat', 'ccmat', 'f', 'ptsall', 'avgfreq', 'area', '-v7.3');

% these are for subsequent analyses
freq_ind=mean(avgfreq,3);
freq_roi=mean(avgfreq,2); 
 

% 
% 
% 
    
% for i = 1:817
%     stdind(i)=std(avgfreq(i,:));
% end
% for i=1:360
%     stdroi(i)=std(avgfreq(:,i));
% end
%  mean(stdind)
%  mean(stdroi)

count=0;
for ind=1:817
    for i = 1:360
        for j=(i+1):360
            count=count+1;
            cctemp(count)=ccmat(1,i, j, ind);
            ftemp(count)=avgfreqmat(1,i,j,ind);
        end
    end
end

temp=reshape(avgfreq,4*817,360);
for i=1:360
    [n,edges]=histcounts(temp(:,i), f(1:38));
    hist_parcel(i,:)=n;
end

 