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

Interested? Try out the example [here](https://github.com/sara-nl/spidercourse/blob/master/extras/cephfs-usage-adv.md).

### High throughput data processing model 

Spider is a high-throughput data-processing platform which means enabling processing of large structured data sets in short time spans. A method to achieve efficient data I/O is to split up the data processing pipelines into many parallel independent jobs where each job retrieves a chunk of data to process on the *local scratch* storage of a worker node (e.g. SSDs). This data processing model is called `embarrassingly parallel` jobs, or else the known `Grid processing model`.

Even if you never used the Grid before, this example will show you how to increase the I/O performance of your jobs by using  *local scratch* and fast network connections between the processing nodes and your input/output data storage locations. Particularly, you will:

- run a data analysis with input/output data located on your project space (on CephFS; Ceph File System)
- copy input/output data to/from a large scratch area on local SSD
- make use of the globally defined variable `$TMPDIR`

Interested? Try out the example [here](https://github.com/sara-nl/spidercourse/blob/master/extras/tmpdir-usage-adv.md).

### Interactive analysis with Jupyter Notebooks

One of the great things about Spider is that -opposite to the Grid- it allows for interactive analysis of large volumes of data. For this purpose we offer Jupyter Notebooks that can be launched on the same powerful high-throughput infrastructure of Spider. The existing Spider users can use Notebooks to run interactively some analysis with data stored/produced on Spider, or debug/prototype their pipelines with software already installed on their project spaces before submitting production runs, or even as a way to pack their work for replication and sharing with their colleagues. The Spider Notebooks is a supporting tool for all this, but not offered as a standalone service for non Spider users or training in classrooms. 

Lets see how to use Jupyter Notebooks on Spider. In this example, you will:

- launch a Notebook and inspect the environment
- install packages within your Notebook or use existing software to run an analysis
- display and publish your results

Interested? Try out the example [here](https://github.com/sara-nl/spidercourse/blob/master/extras/jupyter-usage.md).

### Software portability with containers 

Your analysis on Spider can be run with software that was installed either by the software manager in your project or our sys admin of the system or yourself. What if you want to run the same analysis on another system(s)? Or you want to simply test some workflow on Spider but don't want to install the necessary software from scratch or mess up an existing installation? This is where containers can come in extremely handy. As you do not have admin rights on the system, we do not support Docker containers. But the good news is that we do support `Singularity` containers! In this example, you will:

- use a Singularity image with pre-installed software
- employ the container into your workflows
- run the analysis by using the software from a Singularity container

Interested? Try out the example [here](https://github.com/sara-nl/spidercourse/blob/master/extras/singularity-usage-adv.md).

### Integration with scalable external storage 
(Accessing data from external storage systems)

### Interoperability with existing platforms 
(Software distribution and analysis on Cloud)

### Performance of staging and scratch area

### Bring your own example

If you already understood the basics of Spider usage and you have a good example in mind to fit in this platform, please try it out! You can port an application from your own system to Spider or vice-versa or even think of ways to integrate with other existing platforms at SURFsara or public clouds. 

Need help? Ask for help to the Spider advisors in the room ;)
