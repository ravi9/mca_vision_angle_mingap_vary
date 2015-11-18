function rp_gridMain_Run2Img(img1, img2, minGap)

img1
img2
minGap
tic
startTime = datestr(now)
if (~(strcmp(img1 ,'airfield_edges.pgm.pgm.atts') == 1 || strcmp(img1, '220.pgm.atts') == 1))
	main(['Data' filesep img1], minGap);
end
startTime
endTIme = datestr(now)
toc

tic
startTime = datestr(now)
if (~(strcmp(img2 ,'airfield_edges.pgm.pgm.atts') == 1 || strcmp(img2, '220.pgm.atts') == 1))
	main(['Data' filesep img2], minGap);
end
startTime
endTIme = datestr(now)
toc