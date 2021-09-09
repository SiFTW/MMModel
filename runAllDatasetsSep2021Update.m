%%
%run all datasets

drugsAll={'SGC-CBP30','JQ1','OTX015','ISOX-DUAL','JQ1andSGC-CBP30'};
drugsSingle={'SGC-CBP30','JQ1','OTX015'};
drugsDual={'ISOX-DUAL','JQ1andSGC-CBP30'};
linesAll={'H929','SKMM-1','KMS'};
experimentArray={{drugsAll,linesAll},{drugsSingle,linesAll},{drugsDual,linesAll}};




%general simulation parameters
maxTime=48; %hours
maxTime=maxTime*60;%minutes
maxValOfDrug=5;%fold-inhibition as a result of drug
TimeOfMidDrug=30;%minutes
slopeOfDrug=0.1;%how switch-like the drug turn on is
HL=1; %should drug decay
absThresholdForSS=0.0001;%how accurate steady state needs to be
timeStepForSS=10000;%how long to simulate for one steady state loop. Increase if looping too many times to find SS.
colorArray=[0.7,0,0;0.7,0,0.7;0,0.6,0;0,0.2,0.7;0.5,0.5,0.5;0,0,0];%plot colours.
capWidth=25;%cap width on error bars
plotAbsValus=0;%normalise to timepoint zero
plotIndex=1;

for i=1:length(experimentArray)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% first we model the drug targetting cMyc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 short hl';
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEcMycTarget(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');

    titleString="cMyc targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');
    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
        
    titleString="cMyc targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short MSE";
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now we model the drug targetting IRF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 short hl';
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEIRF4Target(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');
    
    titleString="IRF4 targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');  
    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
  
    titleString="IRF4 targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short- MSE";
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');  
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now we model the drug targetting both
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 short hl';
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEDualTarget(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');
    
    titleString="dual targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');  

    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
   
    titleString="Dual targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL short- MSE";
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% we now repeat everything but with long half life IRF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% first we model the drug targetting cMyc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 long hl'
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEcMycTarget(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');

    

    titleString="cMyc targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL Long";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');  
    
    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
    titleString="cMyc targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- MSE- IRF4 HL long";
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now we model the drug targetting IRF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 long hl';
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEIRF4Target(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');
    
    titleString="IRF4 targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL long";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');   
    
    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
    titleString="IRF4 targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- MSE- IRF4 HL long";
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');   
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now we model the drug targetting both
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     thisExperiment=experimentArray{i};
    drugs=thisExperiment{1}
    lines=thisExperiment{2}
    %% load all experimental data
    allExperimentalData=cell(length(drugs),length(lines));
    for cellIndex=1:length(lines)
        thisLine=lines{cellIndex};
        for drugIndex=1:length(drugs)
            thisDrug=drugs{drugIndex};
            %load experimental data
            experimentalDataFile=strcat('DataForFitting_',thisLine,'_',thisDrug,'.csv');
            fprintf('loading %s file of experimental data\n',experimentalDataFile);
            expDataStruct=loadExperimentalData(experimentalDataFile);
            allExperimentalData{drugIndex,cellIndex}=expDataStruct;
        end
    end 
    fprintf('all experimental data loaded\n');
    
    %lots plot an average of experimental results.
    fig=figure;
    hold on;
    transparency=0.1;
    irf4mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    irf4array=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycmrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    mycarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    blimp1mrnaarray=zeros(length(lines)*length(drugs),length(allExperimentalData{1,1}.time));
    currentIndex=1;

    for cellIndex=1:length(lines)
        for drugIndex=1:length(drugs)
            subplot(5,1,1);
            hold on;
            expDataStruct=allExperimentalData{drugIndex,cellIndex};
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_mrna,expDataStruct.irf4_mrna_std,'lineprops', {'color', colorArray(1,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4mrnaarray(currentIndex,:)=expDataStruct.irf4_mrna;

            subplot(5,1,2);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.irf4_protein,expDataStruct.irf4_protein_std,'lineprops', {'color', colorArray(2,:)},'transparent',1,'patchSaturation',transparency/2);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            irf4array(currentIndex,:)=expDataStruct.irf4_protein;

            subplot(5,1,3);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_mrna,expDataStruct.myc_mrna_std,'lineprops', {'color', colorArray(3,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycmrnaarray(currentIndex,:)=expDataStruct.myc_mrna;

            subplot(5,1,4);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.myc_protein,expDataStruct.myc_protein_std,'lineprops', {'color', colorArray(4,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            mycarray(currentIndex,:)=expDataStruct.myc_protein;

            subplot(5,1,5);
            hold on;
            s=shadedErrorBar(expDataStruct.time,expDataStruct.blimp1_mrna,expDataStruct.blimp1_mrna_std,'lineprops', {'color', colorArray(5,:)},'transparent',1,'patchSaturation',transparency);
            s.edge(1).LineStyle='none';
            s.edge(2).LineStyle='none';
            s.mainLine.LineStyle='none';
            blimp1mrnaarray(currentIndex,:)=expDataStruct.blimp1_mrna;
            %ModelDataStruct=allModelData{drugIndex,cellIndex};
            %plot([0:1:maxTime],modelDataStruct.irf4_mrna,'color',colorArray(1,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.irf4_protein,'color',colorArray(2,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_mrna,'color',colorArray(3,:),'linewidth',3);
            %plot([0:1:maxTime],modelDataStruct.myc_protein,'color',colorArray(4,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_mrna,'color',colorArray(5,:),'linewidth',3);
    %         plot([0:1:maxTime],modelDataStruct.blimp1_protein,'color',colorArray(6,:),'linewidth',3);

        currentIndex=currentIndex+1;
        end
    end
    allExperimentalDataWT=allExperimentalData;
    set(gcf,'color','w');
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    xlabel('hours');
    ylabel('concentration');
    
    %now lets run the model with presumed IRF4 HL
    fprintf("running model\n")
    thisDrug='avg';
    drugTC=loadDrug(thisDrug,maxTime,HL);

    %time course of model
    thisLine='IRF4 long hl';
    params=getAllParamsLine(thisLine);
    v=struct;
    v.params=params;
    initConds=getSSInitialConditions([0,0,0,0,0,0],v,absThresholdForSS,timeStepForSS);
    v.drugTC=drugTC;
    %solve ODE
    [t,x]=ode15s(@(t,x) mmODEDualTarget(t,x,[],v), [0:1:maxTime], initConds,v);

    x(:,1)=x(:,1)./x(1,1);
    x(:,2)=x(:,2)./x(1,2);
    x(:,3)=x(:,3)./x(1,3);
    x(:,4)=x(:,4)./x(1,4);
    x(:,5)=x(:,5)./x(1,5);
    x(:,6)=x(:,6)./x(1,6);

    %IRF4 mRNA
    figure(fig);
    subplot(5,1,1);
    plot([0:1:maxTime],x(:,1),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,1),'color',colorArray(1,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4 mRNA');

    %IRF4
    subplot(5,1,2);
    plot([0:1:maxTime],x(:,2),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,2),'color',colorArray(2,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('IRF4');

    %cMyc mRNA
    subplot(5,1,3);
    plot([0:1:maxTime],x(:,3),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,3),'color',colorArray(3,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc mRNA');

    %cMyc
    subplot(5,1,4);
    plot([0:1:maxTime],x(:,4),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,4),'color',colorArray(4,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('cMyc');

    %Blimp1 mRNA
    subplot(5,1,5);
    plot([0:1:maxTime],x(:,5),'color',[1,1,1,0.5],'linewidth',3);
    plot([0:1:maxTime],x(:,5),'color',colorArray(5,:),'linewidth',2);
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,1.5]);
    xlabel('hours');
    ylabel('concentration');
    title('Blimp1 mRNA');
    titleString="dual targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL long";
    sgtitle(titleString)
    set(gcf,'Position',[100 100 250 1000]);
    export_fig(titleString, '-m4','-png');   
    %%

    %Let's calculate MSEs
    irf4mrnaarray(irf4mrnaarray(:,2)==0,:)=[];
    irf4array(irf4array(:,2)==0,:)=[];
    modelirf4mrna=x([1,4*60,8*60,24*60,48*60],1)';
    modelirf4=x([2,4*60,8*60,24*60,48*60],2)';
    modelcMycmrna=x([3,4*60,8*60,24*60,48*60],3)';
    modelcMyc=x([4,4*60,8*60,24*60,48*60],4)';
    modelBlimp1mrna=x([5,4*60,8*60,24*60,48*60],5)';
    expirfmrnaarray=mean(irf4mrnaarray,1);
    expirf4array=mean(irf4array,1);
    expmycmrnaarray=mean(mycmrnaarray,1);
    expmycarray=mean(mycarray,1);
    expblimp1mrnaarray=mean(blimp1mrnaarray,1);
    irf4mrnadist=(expirfmrnaarray-modelirf4mrna).^2;
    irf4dist=(expirf4array-modelirf4).^2;
    mycmrnadist=(expmycmrnaarray-modelcMycmrna).^2;
    mycdist=(expmycarray-modelcMyc).^2;
    blimp1mrnadist=(expblimp1mrnaarray-modelBlimp1mrna).^2;

    % figure;
    % b=bar([irf4mrnadist;irf4dist;mycmrnadist;mycdist;blimp1mrnadist]','grouped','FaceColor','flat');
    % for k = 1:5
    %     b(k).CData = colorArray(k,:);
    % end

    figure;
    timePoints=[0,4*60,8*60,24*60,48*60];
    hold on;
    area(timePoints,irf4mrnadist,'faceColor',colorArray(1,:),'FaceAlpha',0.5,'EdgeColor',colorArray(1,:));
    area(timePoints,irf4dist,'faceColor',colorArray(2,:),'FaceAlpha',0.5,'EdgeColor',colorArray(2,:));
    area(timePoints,mycdist,'faceColor',colorArray(3,:),'FaceAlpha',0.5,'EdgeColor',colorArray(3,:));
    area(timePoints,mycmrnadist,'faceColor',colorArray(4,:),'FaceAlpha',0.5,'EdgeColor',colorArray(4,:));
    area(timePoints,blimp1mrnadist,'faceColor',colorArray(5,:),'FaceAlpha',0.5,'EdgeColor',colorArray(5,:));
    set(gcf,'color','w');
    ylabel('distance');
    legend({'IRF4 mRNA','IRF4','cMYC mRNA','cMYC','blimp1 mRNA'},'Location','southoutside')
    xticks([0:960:maxTime]);
    xticklabels([0:12:48]);
    ylim([0,0.55]);
    titleString="Dual targetting- Drugs  "+strjoin(drugs,"-")+"- Cells "+strjoin(lines)+ "- IRF4 HL long- MSE"
    title(titleString)
    set(gcf,'Position',[100 100 250 250]);
    export_fig(titleString, '-m4','-png');    
end

