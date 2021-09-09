function params= getAllParamsLine(line)
    params=zeros(1,12);
    %v1=basal cMyc transcription rate
    %mitchell roy 2018 PNAS
    params(1)=0.000040;

    %v2=proportion of cMyc dependent on IRF4
    params(2)=0.8;

    %v3=KD of cMyc transcription wrt IRF4
    params(3)=0.5;
    %v4=basal IRF4 transcription rate
    params(4)=0.0000000927;
    %v5proportion of IRF4 dependent on Myc
    params(5)=0.8;
    %v6=KD of IRF4 transcription wrt cMYC
    params(6)=0.5;

    %v7=cMycTranslation Rate
    params(7)=12;

    %v8=IRF4 translation rate
    %12 proteins per mRNA per min (mitchell roy PNAS 2018)
    params(8)=12;

    %v9=cMycmRNA Degradation rate
    %myC mRNA hl = 60 minutes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2644350/
    params(9)=0.01155245;

    %v10=IRF4mRNA Degradation rate
    %4.315 hours https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2644350/
    %30 mins!
    params(10)=0.0231;

    %v11=cMyc degradation rate
    %30 minutes hl
    params(11)=0.023105;

    %v12=IRF4 degradation rate
    %60 hour half life
    params(12)=0.00165;

    %v13=Blimp1 transcription rate
    params(13)=0.0000000927;

    %v14=Proportion of blimp1 dependent on IRF4
    params(14)=0.8;

    %v15=KD of Blimp1 transcription wrt IRF4
    params(15)=0.01  ;

    %v16=blimp1 translation rate
    params(16)=12;

    %v17 = Blimp1 protein degradation rate
    %blimp 1 protein half life is 4 hours
    % doi:10.1038/s41467-017-00476-w
    % 10.1084/jem.20051611
    params(17)=0.00289;

    %v=18 proportion of IRF dependent on Blimp1
    params(18)=0.8;

    %v=19 KD of IRF transcription wrt to Blimp1
    %NOT USED WE JUST USE 0.9 for the whole thing cMyc+Blimp1
    %params(19)=0.5;

    %v=20 Blimp1 mRNA half-life
    %1 hour https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2644350/
    params(20)=0.01155245;
    
    %coefficient of Blimp1 promotion
    params(21)=4;

    %coefficient of IRF4 promotion
    params(22)=4;
    switch line
        case 'IRF4 long hl'
            params(12)=0.00019254;
            %because IRF concentrations have changed we need to change 
            %the concentration at which Blimp1 responds to IRF4
            params(15)=0.24; 
        case 'IRF4 short hl'
        % everything default            
    end

end
    


