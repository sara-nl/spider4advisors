**Contents of the page:**

* [1. Introduction](#1-introduction)
* [2. Welcome to spider](#2-welcome-to-spider)
  * [i. Getting familiar with the platform](#i-getting-familiar-with-the-platform)
  * [ii. Getting started](#ii-getting-started)
* [3. Spider features](#3-spider-features)
  * [Collaboration within your project](#collaboration-within-your-project)
  * [High throughput data processing model](#high-throughput-data-processing-model)
  * [Interactive analysis with Jupyter Notebooks](#interactive-analysis-with-jupyter-notebooks)
  * [Software portability with containers](#software-portability-with-containers)
  * [Integration with scalable external storage](#integration-with-scalable-external-storage)
  * [Interoperability with existing platforms](#interoperability-with-existing-platforms)
  * [Performance of staging and scratch area](#performance-of-staging-and-scratch-area)
  * [Bring your own example](#bring-your-own-example)
* [4. Feedback](#4-feedback)


---

# 1. Introduction

* Spider backbone:

![Spider backbone](/images/Spider_backbone.png)

* Spider at a glance:

![Spider poster](/images/Spider_poster.png)   


# 2. Welcome to spider

This section includes the prerequisites to start using the platform. It consists of two parts: i) Getting familiar with the platform and ii) Getting started. We ask to finish both parts before moving to the next section with advanced features. By the end of this section you will have:

- A good grasp for the purpose of the platform
- Logged in to the system
- Submitted your first Spider job

### i. Getting familiar with the platform

As an advisor you often search for the best suited platform(s) to match users' requirements.
In order to understand what Spider can offer, we ask you to read through our wiki introduction page [here](http://doc.spider.surfsara.nl/en/latest/Pages/about.html).

### ii. Getting started

Let's start using the platform.  Your project is called `surfadvisors` and your personal account has been sent to your email.
To setup your account, access the platform and run your first Spider job, please follow the instructions in our wiki page [here](http://doc.spider.surfsara.nl/en/latest/Pages/getting_started.html).

# 3. Spider features

This section includes several different features of the platform. It consists of several parts and each part contains an example that can be run independently. You can run the examples in the proposed sequence or simply pick your favourite flavor to start exploring the Spider features!

### Collaboration within your project

The Spider environment is optimised for collaboration. Every user of the platform is member of a project and every project in Spider gets an individual workspace where the project members can collaborate by sharing data, software or workflows.

You may have already made yourself familiar with the project spaces on Spider. If not, don't worry. In this example you will:

- run a data analysis in project spaces
- make use of the project spaces collaborative features  
- understand how role-based access works on Spider  

Interested? Try out the example [here](examples/cephfs-usage.md).

### High throughput data processing model

Spider is a high-throughput data-processing platform which means enabling processing of large structured data sets in short time spans. A method to achieve efficient data I/O is to split up the data processing pipelines into many parallel independent jobs where each job retrieves a chunk of data to process on the *local scratch* storage of a worker node (e.g. SSDs). This data processing model is called `embarrassingly parallel` jobs, or else the known `Grid processing model`.

Even if you never used the Grid before, this example will show you how to increase the I/O performance of your jobs by using  *local scratch* and fast network connections between the processing nodes and your input/output data storage locations. Particularly, you will:

- run a data analysis with input/output data located on your project space (on CephFS; Ceph File System)
- copy input/output data to/from a large scratch area on local SSD
- make use of the globally defined variable `$TMPDIR`

Interested? Try out the example [here](examples/tmpdir-usage.md).

### Interactive analysis with Jupyter Notebooks

One of the great things about Spider is that -opposite to the Grid- it allows for interactive analysis of large volumes of data. For this purpose we offer Jupyter Notebooks that can be launched on the same powerful high-throughput infrastructure of Spider. The existing Spider users can use Notebooks to run interactively some analysis with data stored/produced on Spider, or debug/prototype their pipelines with software already installed on their project spaces before submitting production runs, or even as a way to pack their work for replication and sharing with their colleagues. The Spider Notebooks is a supporting tool for all this, but not offered as a standalone service for non Spider users or training in classrooms.

Let's see how to use Jupyter Notebooks on Spider. In this example, you will:

- launch a Notebook and inspect the environment
- install packages within your Notebook or use existing software to run an analysis
- display and publish your results

Interested? Try out the example [here](examples/jupyter-usage.md).

### Software portability with containers

Your analysis on Spider can be run with software that was installed either by the software manager in your project or our sys admin of the system or yourself. What if you want to run the same analysis on another system(s)? Or you want to simply test some workflow on Spider but don't want to install the necessary software from scratch or mess up an existing installation? This is where containers can come in extremely handy. As you do not have admin rights on the system, we do not support Docker containers. But the good news is that we do support `Singularity` containers! In this example, you will:

- use a Singularity image with pre-installed software
- employ the container into your workflows
- run the analysis by using the software from a Singularity container

Interested? Try out the example [here](examples/apptainer-usage.md).

### Integration with scalable external storage

You can download your raw data on Spider before you start your analysis. However, if you need to analyse data in excess of hundreds of TBs, your data is probably stored elsewhere, most likely on a scalable external storage system. Such a storage system for storing and retrieving huge amounts of data, is `dCache`. `dCache` consists of magnetic tape storage and hard disk storage and both are addressed by a common file system. Wouldn't it be convenient to simply download the data to be analyzed on the fly from dCache, similar to the `Grid processing model` but without using Grid certificates?

This can be achieved thanks to `dCache macaroons` and the high-bandwidth connection between dCache and Spider (up to 1200 Gbit/s). In this example, you will:

- make use of `dCache macaroons` within Spider
- run your analysis by fetching data on the fly directly from dCache
- evaluate options for large scale data analysis automation

Interested? Try out the example [here](examples/macaroons-usage.md).

### Interoperability with existing platforms

Recompiling your software every time you switch processing platforms or moving data around different systems is both time-consuming and makes reproducibility of your work difficult. Spider aims to be a connecting platform and our answer to the interoperability challenges between systems is: `Singularity` for software portability (see [Software portability with containers](#software-portability-with-containers)), `CVMFS/Softdrive` for software distribution and `dCache macaroons` as a user-friendly interface to large storage systems (see [Integration with scalable external storage](#integration-with-scalable-external-storage)). We will combine all this in one example that can run without modifications on HPC Cloud and Spider. In particular, we will:

- wrap our software in a Singularity container
- install the container in a central place and distribute it automagically across different systems with CVMFS/Softdrive
- get our input from dCache with macaroons (no Grid certificates)
- fetch our code from Github
- run the same analysis on multiple platforms, no vendor lock-in

Interested? Try out the example [here](examples/cloud-usage.md).

### The Advanced dCache API (ADA)

dCache has multiple interfaces to approach the data stored on the disks and tapes in the grid storage system. To make it easier for users to access their data SURF has developed the [Advanced dCache API](https://spiderdocs.readthedocs.io/en/latest/Pages/storage/ada-interface.html) (ADA). ADA combines different way to approach your data into a single program. The following example only contains a small set of the commands of ADA. To see the full of ADA scope, read the [ADA documentation](https://spiderdocs.readthedocs.io/en/latest/Pages/storage/ada-interface.html). In this example, we will  use an already existing macaroon to:

- list files on dCache 
- stage a file we want to download
- copy the file to local storage
- unstage the file

Interested? Try out the example [here](examples/ada-usage.md).

### Performance of staging and scratch area

For a high-throughput platform such as Spider important aspects are; (i) the network bandwidth for the transport of data sets, (ii) the I/O operations per second that the storage device can handle, but also (iii) the speed with which we write data and read data from storage. On Spider you have two options for storing files; (a) globally mounted [Ceph](https://ceph.io/) disk storage that is accessed through the shared [CephFS](https://ceph.io/ceph-storage/file-system/) filesystem and (b) local disk storage (often referred to as local scratch space) on the worker node (WN) that is accessed within your jobs only. Your files in /home/[username] and /project/[projectname] are provisioned on Ceph, while the local scratch space for a batch job is provisioned via the number of cores `-c` as 80 Gbyte per core. In this example we will investigate the performance of staging (or CephFS) and scratch area. This is an advanced example if you are interested in performance evaluation as we will:

- discuss technical information on HDD and SSD disk types
- test write/read speeds on Spider on HDDs and SSDs

Interested? Try out the example [here](examples/local-vs-cephfs.md).

### Bring your own example

If you already understood the basics of Spider usage and you have a good example in mind to fit in this platform, please try it out! You can port an application from your own system to Spider or vice-versa or even think of ways to integrate with other existing platforms at SURFsara or public clouds.

Need help? Ask for help to the Spider advisors in the room ;)

# 4. Feedback

We hope that you enjoyed your first journey with Spider! You probably have already gathered some feedback for the platform and we would like to hear this :)

Please fill in the survey [here](https://natalieda.wufoo.com/forms/s1h649s50nu702b/) *by Friday 24th August, end of the day*. It will take you less than 5 minutes but help us a lot!

Next, we invite you for coffee on *Monday at 16:00* to discuss any other comments from your experience with the platform.
