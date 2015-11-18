

a = dir(fullfile('.','Data/*.atts')); % Creates a list of atts files
tic
startTime = datestr(now)

MinGap = 80

%check if result directory exist, else create a resultdirectory
resDir = ['Results_mingap_' num2str(MinGap) ];
if (exist(resDir) ~= 7)
    mkdir(resDir);
end

parfor i=1:length(a), 
%     main(sprintf('Data/%s', a(i).name)); 
    if exist(['Results_mingap_' num2str(MinGap) filesep a(i).name '_AllOutput.mat']) == 0
       if (~(strcmp(a(i).name ,'airfield_edges.pgm.pgm.atts') == 1 || strcmp(a(i).name, '220.pgm.atts') == 1))
          a(i).name
          main(['Data' filesep a(i).name], MinGap)
       end 
    end

end;

RocAnalysis (['Results_mingap_' num2str(MinGap) filesep]);
% RocAnalysis (['Results_mingap_' num2str(MinGap) filesep]);
startTime
endTIme = datestr(now)
toc
%disp(['Time taken for 1 full exp run at mingap ' num2str(MinGap) 'is ' num2str(elapsedTime) ' (sec)']);