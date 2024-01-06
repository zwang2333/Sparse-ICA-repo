#!/bin/bash

cd /home/zwan873/Real_Data_Application/Data/deconfound_sparseICA

#mkdir -p ./LOGS # create directory to save results or output
#mkdir -p ./result

submit_mark=1   # counter, count the index of the nodes

        for a in $(seq 1 400); do   # index for the submitted

                    sed -e s:seedID:"${a}":g </home/zwan873/Real_Data_Application/Data/deconfound_sparseICA/07_deconfounded_singleseed.R >bash_outputs/07_deconfounded_singleseed${a}.R
                    sed -e s:seedID:"${a}":g </home/zwan873/Real_Data_Application/Data/deconfound_sparseICA/07_deconfound_singleseed.sh >bash_outputs/deconfound_singleseed${a}.sh

        
                    # set the limit for the total number of jobs sumitted by YOUR_NAME
                    num=`qstat -u zwan873 | wc -l`
                    while [ $num -gt 80 ]; do   # set the maximum number of jobs to 50
                        sleep 60    # sleep 60 seconds before check the number of jobs again
                        num=`qstat -u zwan873 | wc -l`
                    done
                    
                    
                    
                    if [ ! -e ./DRTMLE_AIPTW/DRTMLE_AIPTW_ASD_TD_z_stat_"$a".RData ]; then    # check to see if there is already an output or result

                    while [ $submit_mark -gt 0 ]; do

                        job_num=`qstat -u zwan873 | grep node4 | wc -l` #190G  # count the number of jobs on node 4
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 1 ]]; then # if the job_num is less than 3 (i.e. $job_num -lt 3), then you can submit jobs to it; if it is larger than 3, you will move on to the next node
                            # "-pe smp 5" request 5 computing thread
                            # "-l h=node4.cluster" request node4 
                            # i am not sure how to set the computing thread for R, make sure the computing thread your job is using matches the computing thread you requested
                            # for now, you just focus on setting the job_num limit at each node, the computing_thread limit at each node, and the totol number of job limit
                            # if you don't want to submit the job to some node, you can change the name to some other node you prefer, but do not delete the submission block
                            # if you want to skip node5, just change node5 to node6 or something else. it is ok if two blocks are sumitting to the same node
                            # you can use "nohup ./run_job.sh &" to let it run in the background, but there is a risk if your jobs are not arranged properly.
                            # if you are absolutely sure, the job arrangement or setup is correct, you can go ahead and use nohup command.
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node4.cluster -l mem_free=15G -l h_vmem=15G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=2
                            # sleep 15 seconds to make sure the submission is complete before submitting the next job, sometimes the cluster is busy, it might take longer to complete the submission, then you need to increase the waiting time accordingly
                            sleep 15
                            break
                        fi
                        # the following block checks if the job submitted to node4 is successful. If it is successful, submit_mark would have been set to 2, then you can move on to submit jobs to node5. If it is not successful, submit_mark would still be 1, in order to move on to node5, you need to set submit_mark to 2, so that it can move on to the second block
                        # pay attension to the last block, the beginning block and the end block should match each other to form a loop
                        if [ $submit_mark -eq 1 ]; then
                            submit_mark=2
                        fi

                        job_num=`qstat -u zwan873 | grep node5 | wc -l` #190G
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 2 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node5.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=3
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 2 ]; then
                            submit_mark=3
                        fi

                        job_num=`qstat -u zwan873 | grep node6 | wc -l` #190G
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 3 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node6.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=4
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 3 ]; then
                            submit_mark=4
                        fi

                        job_num=`qstat -u zwan873 | grep node8 | wc -l` #190G
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 4 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node8.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=5
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 4 ]; then
                            submit_mark=5
                        fi

                        job_num=`qstat -u zwan873 | grep master2 | wc -l` #120G
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 5 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=master2.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=6
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 5 ]; then
                            submit_mark=6
                        fi

                        job_num=`qstat -u zwan873 | grep node15 | wc -l` #120G
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 6 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node15.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=7
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 6 ]; then
                            submit_mark=7
                        fi


                        job_num=`qstat -u zwan873 | grep node16 | wc -l`  #120G
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 7 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node16.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=8
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 7 ]; then
                            submit_mark=8
                        fi

                        job_num=`qstat -u zwan873 | grep node17 | wc -l`  #120G
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 8 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node17.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=9
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 8 ]; then
                            submit_mark=9
                        fi


                        job_num=`qstat -u zwan873 | grep node18 | wc -l`  #120G
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 9 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node18.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=10
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 9 ]; then
                            submit_mark=10
                        fi
                        
                        job_num=`qstat -u zwan873 | grep node19 | wc -l`  #120G 40CPU
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 10 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node19.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=11
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 10 ]; then
                            submit_mark=11
                        fi

                        job_num=`qstat -u zwan873 | grep node20 | wc -l`  #120G 40CPU
                        if [[ $job_num -lt 3 ]] && [[ $submit_mark -eq 11 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node20.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=12
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 11 ]; then
                            submit_mark=12
                        fi

                        job_num=`qstat -u zwan873 | grep node21 | wc -l`  #190G 56CPU
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 12 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node21.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=13
                            sleep 15
                            break
                        fi
                        if [ $submit_mark -eq 12 ]; then
                            submit_mark=13
                        fi

                        # in the last block, you need to reset submit_mark to 1, so that the submission process could restart from node4 at the beginning
                        job_num=`qstat -u zwan873 | grep node20 | wc -l`  #120G 40CPU
                        if [[ $job_num -lt 2 ]] && [[ $submit_mark -eq 13 ]]; then
                            qsub  -R y -N seed${a} -w n -pe smp 1 -cwd -l h=node20.cluster -l mem_free=10G -l h_vmem=10G -j y -o ./bash_outputs/log_"$a" ./bash_outputs/deconfound_singleseed${a}.sh
                            submit_mark=1
                            sleep 15
                            break
                        fi
                        # this is the last block, if the submission to node20 is successful, submit_mark would have been reset to 1. If the submission to node20 is not successful, i.e. submit_mark is still 13, you need to reset it to 1, so that the submission process could restart froom node4 at the beginning
                        # you can also add more blocks to it, make sure that the begining and the end blocks should match each other so that they form a loop
                        if [ $submit_mark -eq 13 ]; then
                            submit_mark=1
                        fi


                        sleep 15

                    done

                    fi
        done

