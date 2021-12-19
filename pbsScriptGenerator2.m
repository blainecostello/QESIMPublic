
Write 10 pbs scripts and write shell script to submit all jobs with the.

runCase = 1;
seedLB = 1;
seedUB = 3;




basefnm = "" + dirnm + "/QESIM_S2_C"+runCase+"_VF";

vfs = [0.01:0.01:0.05];
fnms = strings(1,length(vfs));
for i = 1:length(vfs)
    runLine = "matlab -r "+char(39)+"QESIM_Exp(90, $PBS_ARRAYID, "+ vfs(i) +", 10E-9, [5.4, 1], [2.88E-9, 5.87E4], 10000, 444000)"+char(39)+"\n";
    fnms(i) = basefnm +""+ i + ".pbs"
    fid = fopen(fnms(i),'wt');
    fprintf(fid, "# QE MSTest V4 - C"+runCase+"  Nano Copper Fillers\n");
    fprintf(fid, "#PBS -N QESIM_C"+runCase+"_VF"+i+"\n");
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

write shell script
fid = fopen("" + dirnm + "/RunJobs2_Case"+runCase+".sh",'wt');
fprintf(fid, "#!/bin/bash\n");
for i = 1:length(vfs)
    lineString = "qsub -A GT-jd87 ./" + fnms(i) +"\n";
    fprintf(fid, lineString);
end
fprintf(fid, "qstat -u bcostello3\n");
fclose(fid);



