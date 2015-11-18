MinGap = [20 30 40 50 60 70 80 90 100 110 120 130 140]';

imgDataDir = ['Data'];
imgAttsFileNames = dir([ imgDataDir filesep '*.atts']);

for ii=1:size(MinGap)
	for j=1:size(imgAttsFileNames)
		imgName = imgAttsFileNames(j).name;
		if (~(strcmp(imgName ,'airfield_edges.pgm.pgm.atts') == 1 || strcmp(imgName, '220.pgm.atts') == 1))

			resMatFile = ['Results_mingap_' num2str(MinGap(ii)) filesep imgName '_AllOutput.mat' ];
			if (exist(resMatFile) == 0)
				jobId = [ '_' num2str(j) '_mingap_' num2str(MinGap(ii)) ];
				img2 = '';
				disp(['qsub -N McaVision' jobId ' rp_bashRun2Images.sh ''' imgName ''' ''' img2 ''' ' num2str(MinGap(ii)) ]);
			end
			
		end
	end
end