#!/bin/bash

USERNAME=hduser
HOSTS="hadoop1 hadoop2 hadoop3 hadoop4"
#SCRIPT="sudo apt-get install ncbi-blast+; mkdir /home/hduser/proteins; rm -r /home/hduser/proteins/*"
SCRIPT="mkdir /home/hduser/proteins; rm -r /home/hduser/proteins/*"
for HOSTNAME in ${HOSTS} ; do
    echo ${HOSTNAME}
    ssh -t -o StrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"
    scp /home/hduser/Lab/prosim/proteins.* ${USERNAME}@${HOSTNAME}:/home/hduser/proteins/
    scp /home/hduser/Lab/blast.sh ${USERNAME}@${HOSTNAME}:/home/hduser/proteins/blast.sh
done