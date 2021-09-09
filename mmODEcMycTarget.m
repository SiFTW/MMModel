function delta=mmODEcMycTarget(t,x,options,inputs)
    delta=zeros(size(x));
    IRF4mRNA=x(1);
    IRF4=x(2);
    cMYCmRNA=x(3);
    cMYC=x(4);
    Blimp1mRNA=x(5);
    Blimp1=x(6);

    drugTC=inputs.drugTC;

    drugThisTP=drugTC(max(1,floor(t)));
    v=inputs.params;
    %v1=basal cMyc transcription rate
    %v2=inducible cMyc transcription rate
    %v3=KD of cMyc transcription wrt IRF4
    %v4=basal IRF4 transcription rate
    %v5=inducible IRF4 transcription rate
    %v6=KD of IRF4 transcription wrt cMYC
    %v7=cMycTranslation Rate
    %v8=IRF4 translation rate
    %v9=cMycmRNA Degradation rate
    %v10=IRF4mRNA Degradation rate
    %v11=cMyc degradation rate
    %v12=IRF4 degradation rate

  
    MycNumerator = ((IRF4)/v(3));
    cMycTranscription = v(1)*(MycNumerator/(1+MycNumerator)*(v(2))+(1-v(2)));
    cMycTranscription=cMycTranscription/drugThisTP;


    IRF4Numerator = ((cMYC+Blimp1)/v(6))^v(22);
    IRF4Transcription = v(4)*(IRF4Numerator/(1+IRF4Numerator)*(v(5))+(1-v(5)));
    cMycTranslation=v(7)*cMYCmRNA;
    IRF4Translation=v(8)*IRF4mRNA;
    cMycmRNADegradation=v(9)*cMYCmRNA;
    IRF4mRNADegradation=v(10)*IRF4mRNA;
    cMycDegradation =v(11)*cMYC;
    IRF4Degradation =v(12)*IRF4;

    Blimp1Numerator = ((IRF4)/v(15))^v(21);
    Blimp1Transcription=v(13)*(Blimp1Numerator/(1+Blimp1Numerator)*(v(14))+(1-v(14)));
    Blimp1Translation=v(16)*Blimp1mRNA;
    Blimp1Degradation=v(17)*Blimp1;
    Blimp1mRNADegradation=v(20)*Blimp1mRNA;

    %delta IRF4mRNA
    delta(1)=IRF4Transcription-IRF4mRNADegradation;
    %delta IRF4
    delta(2)=IRF4Translation-IRF4Degradation;
    %delta cMycmRNA
    delta(3)=cMycTranscription-cMycmRNADegradation;
    %delta cMyc
    delta(4)=cMycTranslation-cMycDegradation;
    %delta Blimp1mRNA
    delta(5)=Blimp1Transcription-Blimp1mRNADegradation;
    %delta Blimp1
    delta(6)=Blimp1Translation-Blimp1Degradation;   

end