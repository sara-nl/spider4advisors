## Interactive analysis with Jupyter Notebooks

Let's start our first Spider notebook! Before getting started, we have to set up the environment
to be able to run Jupyter notebooks on the Spider User Interface (UI).

### Pre-requisites

Connect to Spider and install the latest version of mini-conda by fetching and running 
(this script)[https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh]

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

Now create a new conda environment for the notebooks and activate it by doing

```sh
conda create -n notebook
conda activate notebook
```

If you have forgotten your environment name from a previous installation, type

```sh
conda info --envs
```

to get a list of environments.

Now you have a specific environment running in which you can install packages
to run. To run this example, we need to install `jupyterlab` and `matplotlib`.
To do this, run:

```sh
conda install -c conda-forge jupyterlab
conda install matplotlib
```

Before you can start the notebook server, a tunnel has to be created to the
Spider cluster, to be able to forward the notebook contents to your browser.
This is done by opening a new terminal and executing:

```sh
ssh -N -L 8899:UI_NAME:8899 USERNAME@spider.surfsara.nl
```

where `UI_NAME` is the name of the interactive node, `ui-01` or `ui-02`. The 
port is the where the kernel is forwarded through to your browser, and for
`USERNAME` use your own. 

Now start the notebook server on Spider by running in your terminal

```sh
jupyter notebook --no-browser --ip=0.0.0.0 --port=8899
```

and open in your browser the address given, or if the formatting is not 
familiar, go to:
`http://localhost:8899/?token=abc123`
where the token value is given in your terminal.


* Start a notebook with name 'my-awesome-research`

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

* Let's assume that one of your colleagues has prepared a Notebook for you to reproduce his results. Download `my-research-notebook` to your home for further use:

```
mkdir $HOME/my-awesome-notebook
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




