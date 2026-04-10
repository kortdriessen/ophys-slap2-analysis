function pipeLoCo(dr)
disp(['Script started at: ' char(datetime("now"))])

sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams.includeIntegrationROIs = false; 
sParams.sigma_px = 1.33;          
sParams.nmfIter = 2;              
sParams.dXY = 3;                  
sParams.lambda = [];              
sParams.denoiseWindow_s = 0.2;   
sParams.baselineWindow_Glu_s = 4; 
sParams.baselineWindow_Ca_s = 4;  
sParams.activityChannel = 2;      
sParams.tau_s = 0.03;             
sParams.tau2_s = 0.15;            
sParams.maxSynapseDensity = 0.01; 
sParams.nParallelWorkers = 32;    
sParams.drawUserRois = true;     
sParams.motionThresh = 2.5;       
sParams.analyzeHz = 200;          
sParams.nanThresh = 0.33;         
sParams.discardInitial_s = 0;     
sParams.operator = 'KD';       
sParams.makeJSON = false;             

%summarize
sParams.batchSize = 100;
summarize_LoCo_ISOLATED(dr, sParams);

