## Interoperability with existing platforms 

We will be using a genomics pipeline example to test the interoperability features of Spider with other platforms, such as HPC Cloud.

In this example we have pre-configured a Virtual Machine with the necessary components to integrate the analysis with Spider. The components are: `CVMFS` & `Singularity`. In order to access the VM, ask the instructors for your credentials.

> In case that you want to setup the VM yourself you need to setup a machine (e.g. Centos7.6/ 1core/ 4GB RAM):
>  - Install and configure CVMFS - follow the instructions [here](http://doc.grid.surfsara.nl/en/latest/Pages/Advanced/softdrive_on_laptop.html#softdrive-on-laptop)
>  - Install Singularity - follow the instructions [here](https://sylabs.io/guides/3.0/user-guide/installation.html) and install from source Singularity v3.1.0

### Run the analysis

The data we are going 
to use is part of a long-term evolution [experiment](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment)  led by Richard Lenski
to assess adaptation in E. coli. A population was propagated for more than 50,000 
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this 
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the 
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download a small set of data that has already been trimmed and will run the variant calling workflow.

* Login to the pre-configured VM on HPC Cloud (ask for your username/password if you still don't have this). Once logged in download the analysis script and execute it:

```sh
ssh username@spider.usersupport-cloud.surf-hosted.nl 
cd $HOME
mkdir ecoli-analysis-cloud
cd ecoli-analysis-cloud/
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/run-variant-calling-cloud.sh
chmod u+x run-variant-calling-cloud.sh
./run-variant-calling-cloud.sh
```

* While your analysis is running inspect the script `run-variant-calling-cloud.sh`. Where is your input data fetched from? Where is your software installed?

> Food for brain:
> - Run the same analysis on Spider. How would you submit the `run-variant` script to the cluster? Do you have to make any changes to your software or input data paths?  
> - Think of other platforms to port the same analysis, e.g. Cartesius. Sketch your solution.


