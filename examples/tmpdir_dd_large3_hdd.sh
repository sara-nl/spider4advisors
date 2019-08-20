#!/bin/bash
#SBATCH -c 4 -t 4:00:00 --constraint=ivy,hdd
echo "start script: dd read/write tests"
echo "host:"
hostname
echo ""

#--------------------------------------------------------------------------------------------
# CREATE UNIQ LOCAL TEMP DIRECTORIES FOR THE DD COMMANDS
# - $TMPDIR IS TO /tmp on SPIDER this is insufficiently unique when the global cephfs is used
# - Using mktemp in combination with $TMPDIR to create unique temporary directories
#
RUNDIR_UNIQ_LOCAL=`mktemp -d -p ${TMPDIR}`
echo "CREATING A UNIQUE LOCAL TEMPORARY DIRECTORY: ", $RUNDIR_UNIQ_LOCAL
echo ""

# setup shared directory
mkdir -p /home/${USER}${TMPDIR}
RUNDIR_UNIQ_GLOBAL=`mktemp -d -p /home/${USER}${TMPDIR}`
echo "CREATING A UNIQUE GLOBAL TEMPORARY DIRECTORY: ", $RUNDIR_UNIQ_GLOBAL
echo ""
echo ""
#
#exit 0
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# TEST-A :: SMALL 40 GB FILES
# - Note that we should not use the oflag=direct for the global cephfs
#
# LOCAL SCRATCH (Slurm job automatically creates: /scratch/slurm.<jobid>.0)
echo "A - write 10 GB to local HDD scratch space on WN: ", $RUNDIR_UNIQ_LOCAL
date
dd if=/dev/zero of=$RUNDIR_UNIQ_LOCAL/test10GB_4MB_di.img bs=4M count=2500 conv=fdatasync # oflag=direct
echo ""
echo "A - read 10 GB from local HDD scratch space on WN"
date
dd if=$RUNDIR_UNIQ_LOCAL/test10GB_4MB_di.img of=/dev/null bs=4M count=2500
date
#echo "A - clean up the local HDD scratch on WN"
#rm $RUNDIR_UNIQ_LOCAL/test10GB_4MB_di.img
#
echo ""
echo ""
# CEPHFS /HOME
echo "A - write 10 GB to global cephfs (/home/<user>/tmp/): ", $RUNDIR_UNIQ_GLOBAL
date
dd if=/dev/zero of=$RUNDIR_UNIQ_GLOBAL/test10GB_4MB_di.img bs=4M count=2500 conv=fdatasync
echo ""
echo "A - read 10 GB from global cephfs (/home/<user>/tmp/)"
date
dd if=$RUNDIR_UNIQ_GLOBAL/test10GB_4MB_di.img of=/dev/null bs=4M count=2500
date
#echo "A - clean up the global cephfs (/home/<user>/tmp/)"
#rm $RUNDIR_UNIQ_GLOBAL/test40GB_4MB_di.img
#
echo ""
echo ""
echo ""
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# TEST-B :: LARGE 100 GB FILES
#
# LOCAL SCRATCH
echo "B - write 100 GB to local HDD scratch space on WN: ", $RUNDIR_UNIQ_LOCAL
date
dd if=/dev/zero of=$RUNDIR_UNIQ_LOCAL/test100GB_4MB_di.img bs=4M count=25000 conv=fdatasync
echo ""
echo "B - read 100 GB from local HDD scratch space on WN"
date
dd if=$RUNDIR_UNIQ_LOCAL/test100GB_4MB_di.img of=/dev/null bs=4M count=25000
date
echo "B - clean up the local HDD scratch on WN"
rm $RUNDIR_UNIQ_LOCAL/test100GB_4MB_di.img
#
echo ""
echo ""
# CEPHFS /HOME
echo "B - write 100 GB to global cephfs (/home/<user>/tmp/): ", $RUNDIR_UNIQ_GLOBAL
date
dd if=/dev/zero of=$RUNDIR_UNIQ_GLOBAL/test100GB_4MB_di.img bs=4M count=25000 conv=fdatasync
echo ""
echo "B - read 100 GB from global cephfs (/home/<user>/tmp/)"
date
dd if=$RUNDIR_UNIQ_GLOBAL/test100GB_4MB_di.img of=/dev/null bs=4M count=25000
date
echo "B - clean up the global cephfs (/home/<user>/tmp/)"
rm $RUNDIR_UNIQ_GLOBAL/test100GB_4MB_di.img
#
echo ""
echo ""
echo ""
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# TEST-C :: FROM A RE-READ SMALL 40 GB FILES
#
# LOCAL SCRATCH
echo "C - read 10 GB from local HDD scratch space on WN"
date
dd if=$RUNDIR_UNIQ_LOCAL/test10GB_4MB_di.img of=/dev/null bs=4M count=2500
date
echo "C - clean up the local HDD scratch on WN"
rm $RUNDIR_UNIQ_LOCAL/test10GB_4MB_di.img
#
echo ""
echo ""
# USE CEPHFS /HOME
echo "C - read 10 GB from global cephfs (/home/<user>/tmp/)"
date
dd if=$RUNDIR_UNIQ_GLOBAL/test10GB_4MB_di.img of=/dev/null bs=4M count=2500
date
echo "C - clean up the global cephfs (/home/<user>/tmp/)"
rm $RUNDIR_UNIQ_GLOBAL/test10GB_4MB_di.img
#
echo ""
echo ""
echo ""
#--------------------------------------------------------------------------------------------


date
echo ""
echo "end script"
