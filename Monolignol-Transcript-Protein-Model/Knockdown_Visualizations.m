function []=Knockdown_Visualizations(g,TargConst,Yact_table,Xtable,Ypred,legend_flag)

if ~isfield(Ypred,'prevhalf')
    Ypred.prevhalf=zeros(size(Ypred.prevfull));
    prevhalf_flag=0;
else
    prevhalf_flag=1;
end

colors=[180, 180, 180; % Gray 1
    204, 0, 0; % Wolfpack Red
    253, 215, 38; % Hunt Yellow
    209, 73, 5; % Pyroman Flame
    66, 126, 147; % Innovation Blue
    65, 86, 161; % Bioindigo
    230 230 230]./255; % Gray 2


Experiments=Yact_table.Properties.RowNames;
GeneLabels=Yact_table.Properties.VariableNames;
Yact=Yact_table{:,:}';
[M,~]=size(Yact);

%%% Group experiments into their constructs
[Gc, Const]=findgroups(cellfun(@(x) x(6:8),Experiments,'UniformOutput',false));
[Gcl, ConstLines]=findgroups(cellfun(@(x) x(6:end-2),Experiments,'UniformOutput',false));

%%

i=find(strcmp(GeneLabels,['t' g])); %index of transcript for goi
ip=i+M/2; %index of protein for goi
TargConst_ind=find(cellfun(@(x) strcmp(x,TargConst),Const)); %indices of the coi

TargConstLine_ind=find(cellfun(@(x) contains(x,TargConst),ConstLines)); %indices for each line for the coi

WT_ind=find(cellfun(@(x) contains(x,'WT_'),Experiments));
Ywt=Yact_table{WT_ind,1:M}';
Ywt_avg=mean(Ywt,2);


TargGene_inds=find(Gc==TargConst_ind); % index for target gene
it=find(Xtable{TargGene_inds(1),:});
it_prot=it+M/2;
TargGeneNames=GeneLabels(it);

for n=1:length(TargConstLine_ind)
    jinds{n}=find(Gcl==TargConstLine_ind(n)); % indices coi
end


Ymeans=zeros(3,length(TargConstLine_ind)); % mean transcript of interest
Ymeans_p=zeros(4,length(TargConstLine_ind)); % mean protein of interest
Ymeans_t=zeros(length(it),length(TargConstLine_ind)); % mean targeted transcript(s)
Ymeans_tprot=zeros(length(it),length(TargConstLine_ind));% mean targeted protein(s)
for n=1:length(TargConstLine_ind)
    if numel(jinds{n})==1
        
        Ymeans_t(:,n)=Yact(it,jinds{n});
        Ymeans_tprot(:,n)=Yact(it_prot,jinds{n});
        Ymeans(:,n)=[Yact(i,jinds{n}) Ypred.full(i,jinds{n}) Ypred.prevfull(i,jinds{n})];
        Ymeans_p(:,n)=[Yact(ip,jinds{n}) Ypred.full(ip,jinds{n}) Ypred.prevfull(ip,jinds{n}) Ypred.prevhalf(ip,jinds{n})];
        
    else
        
        Ymeans_t(:,n)=mean(Yact(it,jinds{n})');
        Ymeans_tprot(:,n)=mean(Yact(it_prot,jinds{n})');
        Ymeans(:,n)=mean([Yact(i,jinds{n})' Ypred.full(i,jinds{n})' Ypred.prevfull(i,jinds{n})']);
        Ymeans_p(:,n)=mean([Yact(ip,jinds{n})' Ypred.full(ip,jinds{n})' Ypred.prevfull(ip,jinds{n})' Ypred.prevhalf(ip,jinds{n})']);
    end
end

% Calculate WTs
Ymeans_twt=mean([Yact(i,WT_ind)' Ypred.full(i,WT_ind)' Ypred.prevfull(i,WT_ind)']);
Ymeans_pwt=mean([Yact(ip,WT_ind)' Ypred.full(ip,WT_ind)' Ypred.prevfull(ip,WT_ind)' Ypred.prevhalf(ip,WT_ind)']);

% Combine Y_t and Y and Y_t and Y_tp
Ymean_trans=100*[Ymeans_twt; Ymeans']./Ywt_avg(i);
Ymean_prot=100*[Ymeans_pwt; Ymeans_p']./Ywt_avg(ip);

% sort lines by targeted knockdown level
[~, sort_it]=sort(sum(Ymeans_t./Ywt_avg(it),1),'descend');
sort_it=[1 1+sort_it];
Ymean_trans=Ymean_trans(sort_it,:);
Ymean_prot=Ymean_prot(sort_it,:);

wtjinds=[WT_ind jinds];
wtjinds=wtjinds(sort_it);

% plot transcript estimates
figure;
hold on

if prevhalf_flag
    Ymean_trans=[Ymean_trans Ymean_trans(:,1)];
    b=bar(Ymean_trans,'LineWidth',1);
    
    b(1).FaceColor=[.9 .9 .9];
    b(2).FaceColor=[245 194 190]./255;
    b(3).FaceColor=[245 228 190]./255;
    b(4).FaceColor=[255 204 153]./255;
    
    b(1).EdgeColor=colors(1,:);
    b(2).EdgeColor=colors(2,:);
    b(3).EdgeColor=colors(3,:);
    b(4).EdgeColor=colors(4,:);
    
else
    b=bar(Ymean_trans,'LineWidth',1);
    
    b(1).FaceColor=[.9 .9 .9];
    b(2).FaceColor=[245 194 190]./255;
    b(3).FaceColor=[245 228 190]./255;
    
    b(1).EdgeColor=colors(1,:);
    b(2).EdgeColor=colors(2,:);
    b(3).EdgeColor=colors(3,:);
end

pause(1)

for n=1:size(Ymean_trans,1)
    plot(b(1).XData(n)+b(1).XOffset,100*Yact(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(1,:),'MarkerSize',8)
    plot(b(2).XData(n)+b(2).XOffset,100*Ypred.full(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(2,:),'MarkerSize',8)
    plot(b(3).XData(n)+b(3).XOffset,100*Ypred.prevfull(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(3,:),'MarkerSize',8)
    if prevhalf_flag
        plot(b(4).XData(n)+b(4).XOffset,100*Ypred.prevhalf(i,wtjinds{n})./Ywt_avg(i),'ko','MarkerFaceColor',colors(4,:),'MarkerSize',8)
    end
end

wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT'; ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title(GeneLabels(i),'FontSize',22)
if legend_flag
    leg_text={'Exper.','New Model Est.', 'Old Model-S1 Est.'};
    if prevhalf_flag
        leg_text={'Exper.','New Model Est.', 'Old Model-S1 Est.','Old Model-S2 Est.'};
    end
    legend(leg_text,'Location','Best')
end

% Plot Protein Estimates
figure;
hold on
if prevhalf_flag
    b=bar(Ymean_prot,'LineWidth',1);
    
    b(1).FaceColor=[.9 .9 .9];
    b(2).FaceColor=[245 194 190]./255;
    b(3).FaceColor=[245 228 190]./255;
    b(4).FaceColor=[255 204 153]./255;
    
    b(1).EdgeColor=colors(1,:);
    b(2).EdgeColor=colors(2,:);
    b(3).EdgeColor=colors(3,:);
    b(4).EdgeColor=colors(4,:);
    
else
    b=bar(Ymean_prot(:,1:(end-1)),'LineWidth',1);
    
    b(1).FaceColor=[.9 .9 .9];
    b(2).FaceColor=[245 194 190]./255;
    b(3).FaceColor=[245 228 190]./255;
    
    b(1).EdgeColor=colors(1,:);
    b(2).EdgeColor=colors(2,:);
    b(3).EdgeColor=colors(3,:);
end

pause(1)

for n=1:size(Ymean_trans,1)
    plot(b(1).XData(n)+b(1).XOffset,100*Yact(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(1,:),'MarkerSize',8)
    plot(b(2).XData(n)+b(2).XOffset,100*Ypred.full(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(2,:),'MarkerSize',8)
    plot(b(3).XData(n)+b(3).XOffset,100*Ypred.prevfull(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(3,:),'MarkerSize',8)
    if prevhalf_flag
        plot(b(4).XData(n)+b(4).XOffset,100*Ypred.prevhalf(ip,wtjinds{n})./Ywt_avg(ip),'ko','MarkerFaceColor',colors(4,:),'MarkerSize',8)
    end
end

wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT'; ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title(GeneLabels(ip),'FontSize',22)
if legend_flag
    leg_text={'Exper.','New Model Est.', 'Old Model-S1 Est.'};
    if prevhalf_flag
        leg_text={'Exper.','New Model Est.', 'Old Model-S1 Est.','Old Model-S2 Est.'};
    end
    legend(leg_text,'Location','Best')
end

%% Plot Targeted Bar plots
colors_gray=gray(3*length(it)+2);

Ymeans_targ=100*[Ywt_avg(it)'; Ymeans_t']./Ywt_avg(it)';
Ymeans_targprot=100*[Ywt_avg(it_prot)'; Ymeans_tprot']./Ywt_avg(it_prot)';

Ymeans_targ=Ymeans_targ(sort_it,:);
Ymeans_targprot=Ymeans_targprot(sort_it,:);

% Plot Transcript Estimates
figure;
hold on
b=bar(Ymeans_targ,'LineWidth',1);
for nt=1:length(it)
    b(nt).FaceColor=colors_gray(2*length(it)+nt+1,:);
    b(nt).EdgeColor=colors_gray(3*nt+1,:);
end

pause(1)
for nt=1:length(it)
    for n=1:size(Ymean_trans,1)
        plot(b(nt).XData(n)+b(nt).XOffset,100*Yact(it(nt),wtjinds{n})./Ywt_avg(it(nt)),'ko','MarkerFaceColor',colors_gray(3*nt+1,:),'MarkerSize',8)
    end
end
wt_line=refline(0,100);
wt_line.LineStyle='--';
wt_line.LineWidth=0.75;
wt_line.Color='k';
ylabel('Abundance (% of WT)','FontSize',18)
set(gca,'XTick',1:length(TargConstLine_ind)+1,'XTickLabel',['WT';ConstLines(TargConstLine_ind(sort_it(2:end)-1))],'FontSize',18)
set(gcf,'Color','w')
ylimits=ylim;
if ylimits(2)<120
    set(gca,'YLim',[0 120])
end
title('Transcripts of Targeted Genes','FontSize',22)
if legend_flag
    legend(TargGeneNames,'Location','Best')
end
