
% Write 10 pbs scripts and write shell script to submit all jobs with the.

runCase = 2;
seedLB = 1;
seedUB = 4;



basefnm = "QESIM_S3_OC" + runCase + "_VF"

fnms = strings(1,10);
vfs = [0.03:0.03:0.3];
for i = 1:10
    
    if runCase == 1
        runLine = "matlab -r "+char(39)+"QESIM_Exp(90, $PBS_ARRAYID, "+ vfs(i) +", 7E-9, [2.5, 1], [1.5E-7, 5.87E7], 300, 20031)"+char(39)+"\n";
    elseif runCase == 2
        runLine = "matlab -r "+char(39)+"QESIM_Exp_ovlp(90, $PBS_ARRAYID, "+ vfs(i) +", 7E-9, [2.5, 1], [1.5E-8, 5.87E4], -2, 10000, 65003)"+char(39)+"\n";
    elseif runCase == 3
        runLine = "matlab -r "+char(39)+"QESIM_Exp(90, $PBS_ARRAYID, "+ vfs(i) +", 7E-9, [2.5, 1], [1.5E-7, 3.69E7], 300, 20033)"+char(39)+"\n";
    elseif runCase == 4
        runLine = "matlab -r "+char(39)+"QESIM_Exp(90, $PBS_ARRAYID, "+ vfs(i) +", 7E-9, [2.5, 1], [1.5E-7, 3.69E7], 10000, 20034)"+char(39)+"\n";
    end
    
    fnms(i) = basefnm +""+ i + ".pbs"
    fid = fopen(fnms(i),'wt');
    fprintf(fid, "# QE MSTest V3 - C"+runCase+"  Nano Copper Fillers\n");
    fprintf(fid, "#PBS -N QESIM3_OC"+runCase+"_VF"+i+"\n");
    fprintf(fid, "#PBS -l nodes=1:ppn=1\n");
    fprintf(fid, "#PBS -l pmem=512gb\n");
    fprintf(fid, "#PBS -l walltime=12:00:00\n");
    fprintf(fid, "#PBS -q inferno\n");
    fprintf(fid, "#PBS -j oe\n");
    fprintf(fid, "#PBS -t "+seedLB+"-"+seedUB+"\n");
    fprintf(fid, "\n");
    fprintf(fid, "cd $PBS_O_WORKDIR\n");
    fprintf(fid, "module purge\n");
    fprintf(fid, "module load matlab/r2020b\n");
    fprintf(fid, "\n");
    fprintf(fid, runLine);
    fclose(fid);
end

% write shell script
fid = fopen("RunJobs3_Case"+runCase+".sh",'wt');
fprintf(fid, "#!/bin/bash\n");
for i = 1:10
    lineString = "qsub -A GT-jd87 ./" + fnms(i) +"\n";
    fprintf(fid, lineString);
end
fprintf(fid, "qstat -u bcostello3\n");
fclose(fid);



