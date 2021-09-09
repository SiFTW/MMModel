function dataStruct= loadExperimentalData(inputFileName)
    allData=csvread(inputFileName,1,1);
    dataStruct=struct;
    dataStruct.time=allData(:,1)./60;
    dataStruct.myc_mrna=allData(:,2);
    dataStruct.myc_protein=allData(:,3);
    dataStruct.irf4_mrna=allData(:,4);
    dataStruct.irf4_protein=allData(:,5);
    dataStruct.myc_mrna_std=allData(:,6);
    dataStruct.myc_protein_std=allData(:,7);
    dataStruct.irf4_mrna_std=allData(:,8);
    dataStruct.irf4_protein_std=allData(:,9);
    dataStruct.blimp1_mrna=allData(:,10);
    dataStruct.blimp1_mrna_std=allData(:,11); 
end