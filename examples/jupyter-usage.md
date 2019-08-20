## Interactive analysis with Jupyter Notebooks

Let's start our first Spider notebook! To make it easy, we have created a tool to allow users 
launch Jupyter Notebooks. This is called `startnotebook` and can be invoked from the Spider User 
Interface (UI).

* Login to Spider and inspect the `startnotebook`  tool:

```sh
$ startnotebook -h
```

As you see `startnotebook` accepts the following parameters:

| option | description |
| ------ | ----------- |
|`--flavor FLAVOR` | The FLAVOR determines the amount of resources assigned to the Jupyter notebook. Currently the  The available flavours are `small` (1core/8GBmem) and `medium` (4cores/32GBmem). If not specified, the deafult is `small`. |
|`--check` | This is a dry run to inspect the resources to be used by the notebook. It doesn't launch a notebook. |
|` --kernel-timeout KERNEL_TIMEOUT` | Idle time (s) before the kernel is killed. This counts when a kernel is idle. If not specified, the deafult is 10 min. |
|`--notebook-timeout NOTEBOOK_TIMEOUT`| Time (s) before the notebook is killed. This counts when no kernels are running. If not specified, the deafult is 10 min. |
|`--name JOBNAME`| The JOBNAME is the name you want to give to your Notebook folder. If not specified, a random name is created as 'notebook-[random string]' |
|` --verbose`| Verbose mode to display output from background processes that are triggered. | 

* Start a notebook with name 'my-awesome-research`

```sh
$ startnotebook --name my-awesome-notebook
# Starting Notebook.
# Determining Notebook url......
# Notebook sucessfully launched!
# Access your notebook using the following parameters:
#
#	url : http://145.38.255.23:8827
#	pass: OdL@@F>yQ:nQ
```

> What happened? Behind the scenes, `startnotebook` submits a job to the SLURM cluster. The job executes a Singularity container, containing the Jupyter notebook. Once the job starts running the tool
retrieves the host and port where the Jupyter notebook is running and provides you with the corresponding URL, and credentials to allow you access the notebook from your laptop.

* Open a browser in your laptop and copy the provided URL. Hit 'Enter' and paste the provided password. Once successfully logged in, start a new 'Python 3' notebook.

* Check the hidden configuration files created locally in your Spider home account. You should see something like:

```sh
ls -la my-awesome-notebook/
# drwxrwxr-x .config
# drwxr-xr-x .ipynb_checkpoints
# drwxr-xr-x .ipython
# drwxrwxr-x .jupyter
# drwxr-xr-x .local
# -rw-r--r-- 1 Untitled.ipynb
```

> Where is your Notebook running? Can you guess what is the temporary and persistent storage used by your notebooks?

* Let's assume that one of your colleagues has prepared a Notebook for you to reproduce his results. Download `my-research-notebook` to your home for further use:

```
cd $HOME/my-awesome-notebook
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/my-research-notebook.ipynb
ls -l
```

* Run the examples in `my-research-notebook` within your notebook. Each cell can be selected by clicking on it, and can be executed by clicking on the 'Run' icon on the top of the page, or by pressing Shift + Enter.

> Recap and Food for brain:  
> - Spider Notebooks can be used by any user of the Spider platform. Can this feature be requested separetely (i.e. without having a Spider project)?  
> - We are using the Spider's SLURM cluster to power the Jupyter Notebooks. The resources used as long as the Notebook runs, but how are they accounted?  
> - The Notebooks are disposable and not meant to be used for production runs. Idle Notebooks longer that 10 minutes would be forced to shutdown by default. Can this parameter be controled by the user, how?  
> - The user can access his home, project-space data and cvmfs directories within the notebooks. Where should any persistent data  be saved to use after the notebook has been shut down?  
> - The aim is to help Spider users easily launch Notebooks for testing, post-processing, visualisation, tutorials and collaboration. What possible use cases can you think of?  




