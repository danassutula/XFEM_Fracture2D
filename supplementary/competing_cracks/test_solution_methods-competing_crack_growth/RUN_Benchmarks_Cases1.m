
%%% CHOOSE:
flag_loadCase = 1; % load a case
flag_saveCase = 1; % save case results

dir_loadCase = 'cases-1/inputs';
dir_saveCase = 'cases-1/results';

for type_Gs_ = {'uniform','nonunif'}
    type_Gs = type_Gs_{1};
    
    for type_Hs_ = {'indefn','posdef','negdef'}
        type_Hs = type_Hs_{1};
        
        for n_tips = [3,6,12]
            
            NAME = sprintf('GsType(%s)_HsType(%s)_nTips(%i)',...
                type_Gs,type_Hs,n_tips);
            
            % load case
            if flag_loadCase == 1
                load([dir_loadCase,'/',NAME,'_inputs']) % -> INPUTS
            end
             
            % run case
            TestMethod_GradientBasedOptim_MultiTrialCoarsened
            
            if flag_saveCase
                save([dir_saveCase,'/Gradient-MultiTrial-Coarsened/',...
                    NAME,'_results'],'INPUTS','RESULTS');
            end
            
            % run case
            TestMethod_GradientBasedOptim_SingleTrialFixedLength
            
            if flag_saveCase
                save([dir_saveCase,'/Gradient-SingleTrial-FixedLength/',...
                    NAME,'_results'],'INPUTS','RESULTS');
            end
            
        end
    end
end
