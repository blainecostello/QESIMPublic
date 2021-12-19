
% Write 10 pbs scripts and write shell script to submit all jobs with the.
baseFolder = "RS_Davis_1"
mkdir(baseFolder)

runCase = 2;
seedLB = 1;
seedUB = 4;

NumSteps = 20
BaseRunID = 11000

% ~ Code Outline: ~
% First case looks only at permittivity convergence.
% for all different permittivity values,
% and for all volume fractions with each of these values...
% generate a pbs script
% then generate a shell script to run all pbs scripts.

% generate array of parameter to vary.
perms = logspace(-12,9,NumSteps);
vfs = 0.01; %[0.05,0.1,0.15,0.2,0.3];

basefnm = ""+baseFolder+"/QESIM_S1_C" + runCase + "_VF"

fnms = strings(1,length(vfs)*NumSteps);
for j = 1:NumSteps
    for i = 1:length(vfs)
        runID = BaseRunID + j
        thisPerm = perms(j)
        
        runLine = "matlab -r "+char(39)+"QESIM_Exp_NL(90, 1, "+ vfs(i) +", 10E-9, [10, 1], [1E-7,"+thisPerm+"], 50, "+runID+")"+char(39)+"\n";
    
        fnms((i-1)*NumSteps+j) = basefnm +""+ i + "_"+ j +".pbs"
        fid = fopen(fnms((i-1)*NumSteps+j),'wt');
        fprintf(fid, "# QE NL_COMSOL Compare V1 - C"+runCase+"  Nano Copper Fillers\n");
        fprintf(fid, "#PBS -N QESIM1_NL_C"+runCase+"_VF"+i+"\n");
        fprintf(fid, "#PBS -l nodes=1:ppn=1\n");
        fprintf(fid, "#PBS -l pmem=512gb\n");
        fprintf(fid, "#PBS -l walltime=12:00:00\n");
        fprintf(fid, "#PBS -q inferno\n");
        fprintf(fid, "#PBS -j oe\n");
        %fprintf(fid, "#PBS -t "+seedLB+"-"+seedUB+"\n");
        fprintf(fid, "\n");
        fprintf(fid, "cd $PBS_O_WORKDIR\n");
        fprintf(fid, "module purge\n");
        fprintf(fid, "module load matlab/r2020b\n");
        fprintf(fid, "\n");
        fprintf(fid, runLine);
        fclose(fid);
    end
end

% write shell script
fid = fopen(""+baseFolder+"/RunJobs1_CV_Case"+runCase+".sh",'wt');
fprintf(fid, "#!/bin/bash\n");
for ii = 1:length(vfs)
    for jj = 1:NumSteps
        lineString = "qsub -A GT-jd87 ./" + fnms((ii-1)*NumSteps+jj) +"\n";
        fprintf(fid, lineString);
    end
end
fprintf(fid, "qstat -u bcostello3\n");
fclose(fid);



