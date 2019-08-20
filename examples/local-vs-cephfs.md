## Performance of staging and scratch area

### Background information

The internal project that delivers the Ceph storage to Spider is Apollo and it uses HDD (spinning) disk technology to provide Ceph as a storage solution. The read/write speeds of modern HDD disks are typically around 150 Mbyte/sec. Ceph utilizes an array of disks to increase the write/read speeds of data via CephFS to obtain a performance that is higher than that of single disk. This works particularly well for writing data as his is mainly governed by CephFS. It works less well for reading as this is mainly governed by your application. Also note that on Ceph your data has been setup with a 3-fold redundancy to prevent data loss as consequency of disk failures.

The disk storage on the Spider WNs is either HDD or SSD and can be selected via the Slurm `--constraint` option. Here we focus on the local storage in terms of the SSDs only as an example to clarify the difference in speeds obtained when running your data processing on a shared resource vs. a local resource. The SSDs in Spider are of the NVMe type. These NVMe SSD disks provide write speeds up to 3500 Mbyte/sec, whereas standard SATA SSD disks provide write/read speeds of 500/530 Mbyte/sec. The local resource mode for data storage (and processing) for some of you may be reminiscent of the way in which the high-throughput Grid platform is utilized and this is intended as such.

### Testing write/read speeds for data on Spider

We have written a small script to test the write/read speeds on Spider for the two different storage systems. This script uses the lunix command-line utility `dd` (see [wikipedia](https://en.wikipedia.org/wiki/Dd_(Unix)) and [man pages](http://man7.org/linux/man-pages/man1/dd.1.html)) to write data to a file followed by reading this file. It reports the obtained speed(s) in stdout which is caught in your slurm-[jobid].out file. Please note that the results are somewhat idealized.

* Let's run the SSD example:

```sh
cd $HOME
mkdir performance-tests
cd performance-tests/
#download the script and inspect it
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/tmpdir_dd_large3_ssd.sh
chmod u+x tmpdir_dd_large3_ssd.sh
#submit the job for execution
sbatch tmpdir_dd_large3_ssd.sh
```
  
What happened? The job performs 3 tests; A, B and C. In test A it writes/reads a small file (tens of Gbyte) first to/from a local SSD on a WN and then to/from a folder in your /home folder on Ceph. Test B is identical to A, except that it uses files that are much larger (*note this part takes about 20-30m*). Test C repeats only the reading action on the small files from test A.  
  
The expected (approximate) output is:

|Type          | Speed |  
|:-------------|------:|  
|A1 write SSD  | 1 GB/s|  
|A1 read SSD   | 4 GB/s|  
|A2 write Ceph | 1 GB/s|  
|A2 read Ceph  | 5 GB/s|  
|B1 write SSD  | 1 GB/s|  
|B1 read SSD   | 2 GB/s|  
|B2 write Ceph | 1 GB/s|  
|B2 read Ceph  | 0.5 GB/s|  
|C read SSD    | 1 GB/s|  
|C read Ceph   | 0.5 GB/s|  
  
  
> *Food for Brain*    
>  
> Q: Why is there a difference between A2 and B2 ?    
> A: Linux file caching, via the node-based memory, means that the file is still in memory after B2 and hence it is very fast to read in B2. Note that we do not use caching for Ceph iself (this is possible but the benefit is not clear).     
>  
> Q: Ceph has normal HDD disks and not SSD disks, how is it possible that can files be written at speeds greater than about 150 MB/s ?    
> A: The writing to Ceph is optimized (mainly) by Cephfs itself and hence the file writing is distributed over multiple disks.   
>  
> Q: Ceph has normal HDD disks and not SSD, how can files be read at speeds >~200 MB/s?           
> A: By tuning the read-aheads one can can obtain better performance for file reads. As an example a 32 MB (very large) read-ahead implies that when reading 4 MB blocks (standard Ceph block size for our configuration) we can theoretically obtain a performance that is 32/4=8 times faster than the speed of single disk. In practice this condition is (almost) never obtained as the reading procedure is typically more tightly controlled by your application than by CephFS and hence this tuning may not (fully) work.    
>
> Q: Why is the Ceph read speed in C much lower than in A2 ?   
> A: Due to the large B2 file, now being (partially) in the memory cache of the WN, the A2 file has been evicted from the cache and hence the reading can not benefit from this in file-in-memory cache. When evaluating data write/speeds it is very important to account for caching effects.  
>    
> Q: What type of storage/processing solution (local scratch vs. Ceph) would you advise to a new user who needs to process 100's of TB of data with typical data files of about 100 GB in size ?    
> A: TBD.  

* Let's replace the SSD disks by HDD disks and run the example again:

We can replace the local scratch SSD disks in the above dd tests by HDD disks and check their write/read perfomance relative to Ceph. 

```sh
cd $HOME/performance-tests
#download the script and inspect it
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/tmpdir_dd_large3_hdd.sh
chmod u+x tmpdir_dd_large3_hdd.sh
#submit the job for execution
sbatch tmpdir_dd_large3_hdd.sh
```  

What happened now? The script again performs 3 tests; A, B and C. In test A it writes/reads a small file (tens of Gbyte) first to/from a local HDD on a WN and then to/from a folder in your /home folder on Ceph. Test B is identical to A, except that it uses files that are much larger (*note this part takes about 20-30m*). Test C repeats only the reading action on the small files from test A.  
  
The expected (approximate) output is:

|Type          | Speed |  
|:-------------|------:|  
|A1 write HDD  | 0.1 GB/s|  
|A1 read HDD   | 4 GB/s|  
|A2 write Ceph | 1 GB/s|  
|A2 read Ceph  | 5 GB/s|  
|B1 write HDD  | 0.1 GB/s|  
|B1 read HDD   | 0.2 GB/s|  
|B2 write Ceph | 1 GB/s|  
|B2 read Ceph  | 0.5 GB/s|  
|C read HDD    | 4 GB/s|  
|C read Ceph   | 0.5 GB/s|  
  
  
> *Food for Brain*   
> 
> Q: What type of storage/processing solution (local scratch vs. Ceph) would you advise to a new user who needs to process 100's of TB of data with typical data files of about 100 GB in size ?    
> A: TBD
  
