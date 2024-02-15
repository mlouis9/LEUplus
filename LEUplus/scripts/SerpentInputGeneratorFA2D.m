function SerpentInputGeneratorFA2D

% *************************************************************************
% Function that reads which type of fuel assemblies must be loaded
% and generates all the input files automatically (for XS generation)
%
% written by: Dan Kotlyar
%             1/19/2024
% Assisted by Matthew Louis
% *************************************************************************


% *************************************************************************
%                          USER DEFINED OPTIONS
% *************************************************************************
serpentTemplate = '../templateFile/serpentTemplate2D.txt';
dirXS           = '../no_pert_only_burnup/outputXSs';

mgWABAload = 6.03;      % mg/cm loading of B-10 into the WABA
mgIFBAload = 2.25/2.54; % mg/cm loading of B-10 into the WABA


% *************************************************************************



if exist(dirXS,'dir')==7
    disp(['The directory ' dirXS ' already axist and thus is deleted'])
    [status, message, ~] = rmdir(dirXS, 's');
    if ~status
        disp(message)
        return;
    end
    mkdir(dirXS); % create the directory
else
    mkdir(dirXS); % create the directory
end


% --------------------------------------------------------------------------
% Read all the different assemblies configurations
% --------------------------------------------------------------------------
[configFA, pinType] = assemblyConfigurations;




for iXS = 1:length(FAinCore) % loop that generates an input file for each cross-section type
    
    
    idx = identifyFA(keySpec, FAinCore{iXS}); % position/index of the assembly type in 'keySpec'
    %     keySpec(idx).ID    % assembly ID
    %     keySpec(idx).enr   % weight enrichment
    %     keySpec(idx).nFA   % number of this assembly loaded in the core
    %     keySpec(idx).nIFBA % number of IFBA rods in this assembly
    %     keySpec(idx).nWABA % number of WABA rods in this assembly
    %     keySpec(idx).kgU
    if ~idx, continue; end % if this type does not exist in the library
    
    % Serpent input file (for XS generation)
    SINP = ['./' dirXS '/' FAinCore{iXS}];
    copyfile(serpentTemplate, SINP);
    
    txt.title = {['set title ' '"Type= ' keySpec(idx).ID ' ; wt= ' num2str(keySpec(idx).enr) ' ; nIFBA= ' num2str(keySpec(idx).nIFBA) ' ; WABA= ' num2str(keySpec(idx).nWABA) ' ; kgU= ' num2str(keySpec(idx).kgU) '"']};
    printToEOF(SINP, txt.title); % PRINT-TITLE
    
    % Generate a unique fuel assembly configurations with the required #IFBA and
    % #WABA rods
    [currConfig, pin] = generateCustomAssemblyConfig(configFA, pinType, keySpec(idx).nIFBA, keySpec(idx).nWABA);
    txt.array = arrayTxt(currConfig.array);
    printToEOF(SINP, txt.array); % PRINT-ARRAY
    
    % Calculate the weight fraction for each isotope
    [Ui, Oi] = calculateUO2(keySpec(idx).enr);
    
    % Obtain the fuel mat definition (volume and density)
    txt.fuelMAt = defineFuelMat(Ui, Oi, keySpec(idx).kgU, currConfig.volFP);
    txt.fuelPin = pinTxt(pin,'fp'); % obtain the fuel pin universe definition
    
    printToEOF(SINP, txt.fuelMAt); % PRINT-FUEL MATERIAL
    printToEOF(SINP, txt.fuelPin); % PRINT-FUEL PIN/UNIVERSE
    
    % Obtain the IFBA mat definition (volume, density and isotopics)
    if keySpec(idx).nIFBA % only if the number of IFBA rods not zero
        ifbaNDi = calculateIFBA_ND(mgIFBAload, pin); % B-10, Zr concentrations
        txt.ifbaMat = defineIFBaMat(currConfig.volIB,ifbaNDi);
        txt.ifbaPin = pinTxt(pin,'if'); % obtain the ifba pin universe definition
        
        printToEOF(SINP, txt.ifbaMat); % PRINT-MATERIAL
        printToEOF(SINP, txt.ifbaPin); % PRINT-PIN/UNIVERSE
    end
    % Obtain the WABA mat definition (volume, density and isotopics)
    if keySpec(idx).nWABA % only if the number of IFBA rods not zero
        wabaNDi = calculateWABA_ND(mgWABAload, pin); % B-10, C, Al, O concentrations
        txt.wabaMat = defineWABaMat(currConfig.volWB,wabaNDi);
        txt.wabaPin = pinTxt(pin,'wb'); % obtain the ifba pin universe definition
        
        printToEOF(SINP, txt.wabaMat); % PRINT-MATERIAL
        printToEOF(SINP, txt.wabaPin); % PRINT-PIN/UNIVERSE
    end
    
    
end













function [configFA, pinType] = assemblyConfigurations

% -------------------------------------------------------------------------
% Function that stores all the different fuel assembly configurations
% #IFBA and #WABA are different from file to file
% -------------------------------------------------------------------------

% legend for the different type of pins within the assembly
f = 'fp';  % fuel pin
gt= 'gt';  % guide tube
in ='it';  % instrumental tube
wb= 'wb';  % waba tube
i = 'if';  % ifba tube
pinType.f = f; pinType.gt = gt; pinType.in = in; pinType.wb = wb; pinType.i = i;


idx = 0;
idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  in f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f};

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  in f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f};

idx = idx + 1;
configFA(idx).array = {
    i  f  i  f  i  f  i  i  f  i  i  f  i  f  i  f  i
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    i  f  i  i  i  gt i  i  gt i  i  gt i  i  i  f  i
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    i  f  i  i  i  i  f  i  i  i  f  i  i  i  i  f  i
    f  i  gt i  i  gt i  i  gt i  i  gt i  i  gt i  f
    i  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  i
    i  f  i  f  i  i  f  i  i  i  f  i  i  f  i  f  i
    f  i  gt i  i  gt i  i  in i  i  gt i  i  gt i  f
    i  f  i  f  i  i  f  i  i  i  f  i  i  f  i  f  i
    i  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  i
    f  i  gt i  i  gt i  i  gt i  i  gt i  i  gt i  f
    i  f  i  i  i  i  f  i  i  i  f  i  i  i  i  f  i
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    i  f  i  i  i  gt i  i  gt i  i  gt i  i  i  f  i
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    i  f  i  f  i  f  i  i  f  i  i  f  i  f  i  f  i};

idx = idx + 1;
configFA(idx).array = {
    i  f  f  f  i  f  f  i  f  i  f  f  i  f  f  f  i
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  i  f  i  i  gt i  i  gt i  i  gt i  i  f  i  f
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    i  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  i
    f  i  gt i  i  gt i  i  gt i  i  gt i  i  gt i  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    i  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  i
    f  i  gt i  i  gt i  i  in i  i  gt i  i  gt i  f
    i  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  i
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  i  gt i  i  gt i  i  gt i  i  gt i  i  gt i  f
    i  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  i
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    f  i  f  i  i  gt i  i  gt i  i  gt i  i  f  i  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    i  f  f  f  i  f  f  i  f  i  f  f  i  f  f  f  i};

idx = idx + 1;
configFA(idx).array = {
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    f  f  f  i  i  gt i  i  gt i  i  gt i  i  f  f  f
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    f  i  i  i  f  i  f  f  i  f  f  i  f  i  i  i  f
    f  f  gt i  i  gt i  i  gt i  i  gt i  i  gt f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  i  gt i  i  gt i  i  in i  i  gt i  i  gt i  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  gt i  i  gt i  i  gt i  i  gt i  i  gt f  f
    f  i  i  i  f  i  f  f  i  f  f  i  f  i  i  i  f
    f  f  i  gt i  i  f  f  i  f  f  i  i  gt i  f  f
    f  f  f  i  i  gt i  i  gt i  i  gt i  i  f  f  f
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i};

idx = idx + 1;
configFA(idx).array = {
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    f  f  i  f  i  gt f  i  gt i  f  gt i  f  i  f  f
    f  f  f  gt i  i  f  f  f  f  f  i  i  gt f  f  f
    f  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  f
    f  i  gt i  i  gt f  i  gt i  f  gt i  i  gt i  f
    f  f  f  f  f  f  i  f  f  f  i  f  f  f  f  f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  i  gt f  i  gt f  i  in i  f  gt i  f  gt i  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  f  f  f  f  i  f  f  f  i  f  f  f  f  f  f
    f  i  gt i  i  gt f  i  gt i  f  gt i  i  gt i  f
    f  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  f
    f  f  f  gt i  i  f  f  f  f  f  i  i  gt f  f  f
    f  f  i  f  i  gt f  i  gt i  f  gt i  f  i  f  f
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i};

idx = idx + 1;
configFA(idx).array = {
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    f  f  f  f  i  gt f  i  gt i  f  gt i  f  f  f  f
    f  f  f  gt i  f  f  f  f  f  f  f  i  gt f  f  f
    f  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  f
    f  i  gt f  i  gt f  i  gt i  f  gt i  f  gt i  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  i  gt f  i  gt f  i  in i  f  gt i  f  gt i  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  i  gt f  i  gt f  i  gt i  f  gt i  f  gt i  f
    f  f  i  i  f  i  f  f  i  f  f  i  f  i  i  f  f
    f  f  f  gt i  f  f  f  f  f  f  f  i  gt f  f  f
    f  f  f  f  i  gt f  i  gt i  f  gt i  f  f  f  f
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i};

idx = idx + 1;
configFA(idx).array = {
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  i  gt f  i  gt i  f  gt i  f  f  f  f
    f  f  f  gt i  f  f  f  f  f  f  f  i  gt f  f  f
    f  f  i  i  f  i  f  f  f  f  f  i  f  i  i  f  f
    f  f  gt f  i  gt f  i  gt i  f  gt i  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  gt f  f  gt f  i  in i  f  gt f  f  gt f  f
    f  f  i  f  f  i  f  f  i  f  f  i  f  f  i  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  i  gt f  i  gt i  f  gt i  f  gt f  f
    f  f  i  i  f  i  f  f  f  f  f  i  f  i  i  f  f
    f  f  f  gt i  f  f  f  f  f  f  f  i  gt f  f  f
    f  f  f  f  i  gt f  i  gt i  f  gt i  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    i  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  i};

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  i  f  f  f  f  f  i  f  f  f  f  f
    f  f  f  f  i  gt f  f  gt f  f  gt i  f  f  f  f
    f  f  f  gt f  f  f  f  i  f  f  f  f  gt f  f  f
    f  f  f  i  f  f  f  f  f  f  f  f  f  i  f  f  f
    f  i  gt f  f  gt i  f  gt f  i  gt f  f  gt i  f
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt i  f  gt i  f  in f  i  gt f  i  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  i  f  f  i  f  f  i  f  f  f  f  f
    f  i  gt f  f  gt i  f  gt f  i  gt f  f  gt i  f
    f  f  f  i  f  f  f  f  f  f  f  f  f  i  f  f  f
    f  f  f  gt i  f  f  f  i  f  f  f  i  gt f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  f  f  i  f  f  f  f  f  i  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f};

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  wb f  f  wb f  f  wb f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  wb f  f  wb f  f  wb f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  wb f  f  in f  f  wb f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  wb f  f  wb f  f  wb f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  wb f  f  wb f  f  wb f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  wb f  f  wb f  f  wb f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  wb f  f  in f  f  wb f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  wb f  f  wb f  f  wb f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  wb f  f  gt f  f  wb f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  wb f  f  in f  f  wb f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  wb f  f  gt f  f  wb f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  gt i  f  gt f  i  gt f  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  i  f  f  f  f  f  f  f  i  f  f  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  i  f  f  f  f  f  f  f  f  f  f  f  i  f  f
    f  f  f  f  f  f  f  i  f  i  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  in f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  i  f  i  f  f  f  f  f  f  f
    f  f  i  f  f  f  f  f  f  f  f  f  f  f  i  f  f
    f  f  gt f  f  gt f  f  gt f  f  gt f  f  gt f  f
    f  f  f  f  i  f  f  f  f  f  f  f  i  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  f  gt i  f  gt f  i  gt f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  wb f  f  gt f  f  wb f  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  wb f  f  in f  f  wb f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  wb f  f  gt f  f  wb f  f  gt f  f  wb f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  gt f  f  f  f  f  f  f  f  f  gt f  f  f
    f  f  f  f  f  wb f  f  gt f  f  wb f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };

idx = idx + 1;
configFA(idx).array = {
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  wb f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  wb f  f  in f  f  wb f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  gt f  f  gt f  f  wb f  f  gt f  f  gt f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  wb f  f  f  f  f  f  f  f  f  wb f  f  f
    f  f  f  f  f  gt f  f  gt f  f  gt f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f
    f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f  f };


% ///////////////////////////////////////////////////
% Calculate the number of pins of each type
% ///////////////////////////////////////////////////

for iFA = 1:length(configFA)
    configFA(iFA).nFP = sum(sum(ismember(configFA(iFA).array,f)));
    configFA(iFA).nGT = sum(sum(ismember(configFA(iFA).array,gt)));
    configFA(iFA).nIT = sum(sum(ismember(configFA(iFA).array,in)));
    configFA(iFA).nWB = sum(sum(ismember(configFA(iFA).array,wb)));
    configFA(iFA).nIB = sum(sum(ismember(configFA(iFA).array,i)));
end




function [currConfig, pin] = generateCustomAssemblyConfig(configFA, pinType, nIFBA, nWABA)

% -------------------------------------------------------------------------
% Function that generates a fuel assembly configuration with a specific set
% #IFBA
% #WABA rods
% -------------------------------------------------------------------------

% legend for the different type of pins within the assembly
f = pinType.f;   % fuel pin
wb= pinType.wb;  % waba tube
i = pinType.i;   % ifba tube

nFP = 264;       % total number of fuel pins

% define each type of pin separately
idx = 0;  % counter for new pin types
idx = idx + 1;  % new universe/pin type --> fuel pin
pin(idx).mat = {'fuel','helium', 'clad','water'};
pin(idx).rad = [0.392176 0.40005 0.4572 inf];
pin(idx).type = 'fp';
idx = idx + 1;  % new universe/pin type --> ifba pin
pin(idx).mat = {'fuel','zrb2', 'helium','clad','water'};
pin(idx).rad = [0.392176 0.393176 0.40005 0.4572 inf];
pin(idx).type = 'if';
idx = idx + 1;  % new universe/pin type --> waba pin
pin(idx).mat = {'water','clad', 'helium','waba', 'helium','clad','water'};
pin(idx).rad = [0.29 0.34 0.35 0.40386 0.41783 0.48387 inf];
pin(idx).type = 'wb';
idx = idx + 1;  % new universe/pin type --> guide tube
pin(idx).mat = {'water','clad', 'water'};
pin(idx).rad = [0.56134 0.61214 inf];
pin(idx).type = 'gt';
idx = idx + 1;  % new universe/pin type --> inst. tube
pin(idx).mat = {'air','clad', 'water'};
pin(idx).rad = [0.56134 0.61214 inf];
pin(idx).type = 'it';



if sum([configFA.nIB] == nIFBA)
    ifbaConfig = configFA([configFA.nIB] == nIFBA).array;
else
    disp(['A configuration with ' num2str(nIFBA) ' IFBA rods does NOT exist'])
    return;
end
if sum([configFA.nWB] == nWABA)
    wabaConfig = configFA([configFA.nWB] == nWABA).array;
else
    disp(['A configuration with ' num2str(nWABA) ' WABA rods does NOT exist'])
    return;
end

% Create a custom assembly configuraion from knowing the number of IFBA and
% WABA rods:
currConfig.array = ifbaConfig;
for i = 1:size(ifbaConfig,1)
    for j = 1:size(ifbaConfig,2)
        if strcmpi(wabaConfig{i, j},pinType.wb), currConfig.array{i, j}=pinType.wb; end
    end
end

% ///////////////////////////////////////////////////
% Calculate the volume of each material
% ///////////////////////////////////////////////////

currConfig.volFP = 0; % volume of fuel
currConfig.volIB = 0; % volume of IFBA
currConfig.volWB = 0; % volume of WABA
for iPin = 1:length(pin) % find a specific pin type to calculate its volume
    
    if strcmpi(pin(iPin).type, pinType.f)
        currConfig.volFP = currConfig.volFP+...
            (nFP-nIFBA)*pi*(pin(iPin).rad(1))^2;
    end
    if strcmpi(pin(iPin).type, pinType.i)
        currConfig.volIB = currConfig.volIB+...
            nIFBA*pi*(pin(iPin).rad(2)^2-pin(iPin).rad(1)^2);
        
        currConfig.volFP = currConfig.volFP+...
            nIFBA*pi*(pin(iPin).rad(1))^2;
    end
    
    if strcmpi(pin(iPin).type, pinType.wb)
        currConfig.volWB = currConfig.volWB+...
            nWABA*pi*(pin(iPin).rad(4)^2-pin(iPin).rad(3)^2);
    end
    
end




function [keySpec] = fuelAssemblyKeysReader(filename)

% -------------------------------------------------------------------------
% Function that reads all the specs for each fuel assembly type
% -------------------------------------------------------------------------
fidi = fopen(filename,'r');
idx = 1; % counter to calculate the number of different fuel assemblies
keySpec = struct(); % structure to contain all specs
fgetl(fidi);
while ~feof(fidi)
    tline = fgetl(fidi);
    [str1, strVals]=strtok(tline);
    if strcmpi(str1,''), continue; end % if this is an empty line
    vals = str2num(strVals);  % convert to numeric values
    if length(vals)~=5, disp(['Data is missing in the keys fiel for assembly type ' str1]); return; end
    
    keySpec(idx).ID    = str1;
    keySpec(idx).enr   = vals(1); % weight enrichment
    keySpec(idx).nFA   = vals(2); % number of this assembly loaded in the core
    keySpec(idx).nIFBA = vals(3); % number of IFBA rods in this assembly
    keySpec(idx).nWABA = vals(4); % number of WABA rods in this assembly
    keySpec(idx).kgU   = vals(5); % kg of HM
    idx = idx + 1;
end
fclose(fidi);


function [idx] = identifyFA(keySpec, currFA)
% ----------------------------------------------------
% Function that identifies for a specific assembly type
% ----------------------------------------------------
idx = 0;
for idx0 = 1:length(keySpec)
    if strcmpi(keySpec(idx0).ID, currFA), idx = idx0; break; end
end


function [Ui, Oi] = calculateUO2(wf)
% ----------------------------------------------------
% Calculate the isotopic composition for UO2
% ----------------------------------------------------
format long
mu4 = 234.040952;
mu5 = 235.043930;
mu8 = 238.050788;
mo16 = 15.9949146196;
mo17 = 16.9991317;

% --- Calculate mass fraction of each element in fuel ---
x4 = 0.008 * wf;
x5 = wf;
x8 = 100 - (x5 + x4);

mU = (1/((x4/mu4) + (x5/mu5) + (x8/mu8))) * 100;
mO = (mo16 * 0.99757) + (mo17 * 0.00038);

wU = mU / (mU + 2*mO);
wO = 1 - wU;

u4 = (x4 * wU) / 100;
u5 = (x5 * wU) / 100;
u8 = (x8 * wU) / 100;
Ui = [u4, u5, u8];

o16 = wO * 0.99757;
o17 = wO * 0.00243;
Oi = [o16, o17];


function [Ni] = calculateWABA_ND(mgWABAload, pin)

% Isotopic molar masses (gr/mol)
A.B10  = 10.01294;
A.B11  = 11.00931;
A.C12  = 12.00000;
A.O16  = 15.99990;
A.Al27 = 26.98154;

avogadro  = 0.602214086;
densAl2O3 = 3.95; % gr/cm**3
posAl2O3  = 0.95;  % porosity factor
wtAl2O3   = 0.86;  % weight fraction of alominia in the mixture
awAl2O3   = 2*A.Al27 + 3*A.O16;

idx = []; % store the identifier for the
S = 0;    % radial area (cm**2) of the Al2O3-B4C material
for iMat=1:length(pin)
    idx = find(strcmpi(pin(iMat).mat,'waba'));
    if ~isempty(idx)
        S = pi*(pin(iMat).rad(idx)^2-pin(iMat).rad(idx-1)^2);
        break;
    end
end

% Calculate the atomic density #/b/cm for B10 and C
NB10 = (mgWABAload/1000)*avogadro / (S*A.B10);
% Since the poison matrial is B4C then N(C)=N(B10)/4
Nc = NB10 / 4;

% evaluate the concentration of Al27 and O16
N_Al2O3 = wtAl2O3*densAl2O3*posAl2O3*avogadro/awAl2O3;
Nal = 2*N_Al2O3;
No  = 3*N_Al2O3;

Ni = [NB10;Nc;Nal;No];


function [Ni] = calculateIFBA_ND(mgIFBAload, pin)

% Isotopic molar masses (gr/mol)
A.B10  = 10.012937;
avogadro  = 0.602214086;

idx = []; % store the identifier for the
S = 0;    % radial area (cm**2) of the Al2O3-ZrB2 material
for iMat=1:length(pin)
    idx = find(strcmpi(pin(iMat).mat,'zrb2'));
    if ~isempty(idx)
        S = pi*(pin(iMat).rad(idx)^2-pin(iMat).rad(idx-1)^2);
        break;
    end
end

% Calculate the atomic density #/b/cm for B10 and Zr
NB10 = (mgIFBAload/1000)*avogadro / (S*A.B10);
% Since the poison matrial is ZrB2 then N(Zr)=N(B10)/2
Nzr = NB10 / 2;

Ni = [NB10;Nzr];


function [matTxt] = defineFuelMat(Ui, Oi, kgU, unitVol)

% -------------------------------------------------------------------------
% Define the material density and volume
% -------------------------------------------------------------------------
format long
Height = 144 * 2.54; % cm - assembly height
massUO2 = kgU/sum(Ui);

volUO2 = unitVol * Height; % total volume is only used to calculate the density of the fuel
densUO2 = massUO2*1E+3 / volUO2; % weight density in g/cm**3

matTxt = {['mat fuel -' num2str(densUO2) ' vol ' num2str(unitVol) ' burn 1'];
    ['92234.09c  -' num2str(Ui(1),'%6.6e')];
    ['92235.09c  -' num2str(Ui(2),'%6.6e')];
    ['92238.09c  -' num2str(Ui(3),'%6.6e')];
    [' 8016.09c  -' num2str(Oi(1),'%6.6e')];
    [' 8017.09c  -' num2str(Oi(2),'%6.6e')];};


function [matTxt] = defineIFBaMat(unitVol,ND)
% -------------------------------------------------------------------------
% Define the material volume in the IFBA material card
% -------------------------------------------------------------------------
matTxt = {['mat zrb2 sum vol ' num2str(unitVol) ' burn 1'];
    [' 5010.06c  ' num2str(ND(1),'%6.6e')];
    ['40000.06c  ' num2str(ND(2),'%6.6e')]};

function [matTxt] = defineWABaMat(unitVol,ND)
% -------------------------------------------------------------------------
% Define the material volume in the IFBA material card
% -------------------------------------------------------------------------
matTxt = {['mat waba sum vol ' num2str(unitVol) ' burn 1'];
    [' 5010.06c  ' num2str(ND(1),'%6.6e')];
    [' 6012.06c  ' num2str(ND(2),'%6.6e')];
    ['13027.06c  ' num2str(ND(3),'%6.6e')];
    [' 8016.06c  ' num2str(ND(4),'%6.6e')]};


function [txtPrint] = pinTxt(pin,type)

for idx=1:length(pin)
    if strcmpi(pin(idx).type,type), break; end
end

txtPrint{1,:} = ['pin ' pin(idx).type];
for j=1:length(pin(idx).mat)
    if ~isinf(pin(idx).rad(j))
        txtPrint{j+1,:} = [pin(idx).mat{j} '  ' num2str(pin(idx).rad(j),'%7.8f')];
    else
        txtPrint{j+1,:} = [pin(idx).mat{j}];
    end
end

function [txtPrint] = arrayTxt(arrayFA)

txtPrint{1,:} = 'lat 1 1 0.0 0.0 17 17 1.25984';
for j=1:size(arrayFA,1)
    currLine = [];
    for i=1:size(arrayFA,2)
        currLine = [currLine ' ' arrayFA{j,i}];
    end
    txtPrint{j+1,:} = currLine;
end


function printToEOF(filename, txt)

% -------------------------------------------------------------------------
% Function that prints a specific text to the end of the file
% -------------------------------------------------------------------------
fid = fopen(filename,'a');
fprintf(fid, '\n');
for i=1:length(txt)
    fprintf(fid, '%s\n', txt{i,:});
end
fprintf(fid, '\n');
fclose(fid);

