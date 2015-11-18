function y = RocAnalysis (Directory)

%% This is a standaloine function that is supposed to run AFTER "main.m"
%% creates all the *.mat files, each of which stores the ROC information
%% for an image along with other variables
%% The default directory to look at is passed an argument.
%% The function will plot an aggregate ROC, which is the median ROC with
%% the median of the absolute differences.

Files = dir(fullfile('.', sprintf('%s/*.mat', Directory)));

M_Time = []; T_Time = []; N =[];

FalseAlarm = [0:0.01:1.0];
DetectionRates = [];
for i=1:length(Files)
    load(sprintf('%s/%s', Directory, Files(i).name)); 
    %plot(Pf, Pd); hold on; 
    
    %Remove duplicate points, i.e. different Pd for same Pf -- replace by
    %maximum Pd for those points
    n = length(Pf);
    Pff = [Pf(1) ]; Pdd = [Pd(1)];
    for (j=2:length(Pf))
        if (Pf(j) ~= Pf(j-1))
            Pff = [Pff Pf(j) ]; Pdd = [Pdd Pd(j)];
        end;
    end;
    Pff = [Pff 1]; Pdd = [Pdd 1.0];
    Dr = interp1(Pff, Pdd, FalseAlarm,'linear');
    DetectionRates = [DetectionRates; Dr];
    T_Time = [T_Time t_vision];
    M_Time = [M_Time t_magnet];
    N = [N length(MCoupled)];
    %plot(Pff, Pdd,'g'); plot(FalseAlarm, Dr, 'y');
    %roc = @(x)norm(Pdd - SemiparametricRoc(x(1), x(2), Pff))^2;
    %[x,fval,exitflag] = fminsearch(roc,[1, 1]);
    %DetectionRates = [DetectionRates; SemiparametricRoc(x(1), x(2), FalseAlarm)];
    %plot(Pff, Pdd,'g'); plot(FalseAlarm, SemiparametricRoc(x(1), x(2), FalseAlarm), 'y');
    %pause;
end;


y = DetectionRates;

errorbar(FalseAlarm, mean(DetectionRates), std(DetectionRates));
xlabel('False Positive Rate'); ylabel('True Detection Rate');
axis([0 1 0 1]);

figure;
[N I] = sort(N, 'ascend');
M_Time = M_Time(I); T_Time = T_Time(I);
p = polyfit(N, M_Time, 2), Mp_Time = polyval(p, N);
p = polyfit(N, T_Time, 3), Tp_Time = polyval(p, N);

plot(N, M_Time, 'bo'); hold on; plot(N, Mp_Time, 'b', 'LineWidth', 2);
plot (N, T_Time, 'rx'); plot(N, Tp_Time, 'r', 'LineWidth', 2);
xlabel('Problem Size - # of image primitives');
ylabel('Wall clock time (seconds)');

fprintf(1, '\n Time Stats: Traditional Vision = %f +- %f (sec), Magnets: %f +- %f (sec)', mean(T_Time), std(T_Time), mean(M_Time), std(M_Time));

function Pd = SemiparametricRoc(mu, sigma, Pf)
Pd = normcdf(mu/sigma + norminv(Pf)/sigma);
    