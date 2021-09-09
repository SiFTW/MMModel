function drugTC = loadDrug(drugName,maxTime,HL)
    fprintf('generating %s drug time course.\n',drugName);
    
    slopeOfDrug=0.8;
    maxValOfDrug=9;
    TimeOfMidDrug=5;
    drugScalar=1;
    switch drugName
        case 'SGC-CBP30'
            %IC50 = 3.47
            maxValOfDrug=drugScalar*1.527;
        case 'OTX015'
            %IC50 = 0.47uM
            maxValOfDrug=drugScalar*8.51;
        case 'JQ1'
            %IC50 = 0.31uM
            maxValOfDrug=drugScalar*12.90;
            
        case 'ISOX-DUAL'
            %IC50 4uM
            maxValOfDrug=drugScalar;
        case 'JQ1andSGC-CBP30'
            %ic50=0.28
            maxValOfDrug=drugScalar*14.2857;
        case 'avg'
            %ic50=0.28
            maxValOfDrug=drugScalar*6;
        otherwise
            fprintf('unrecognized drug\n');
    end
    
    if HL
        drugTC=generateDrugTC_hl(maxValOfDrug,TimeOfMidDrug,slopeOfDrug,maxTime)+1;
    else
        drugTC=generateDrugTC(maxValOfDrug,TimeOfMidDrug,slopeOfDrug,maxTime)+1;
    end
end
