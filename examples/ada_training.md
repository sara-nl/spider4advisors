# ADA demo session

## Preparation

### 1. What you need to run the examples

- Access to a terminal (Spider UI or any other machine) and internet connection
- ADA wrapper (is available on Spider. On other machines you need to download it from: `git clone https://github.com/sara-nl/SpiderScripts.git`)
- Rclone (available on any SURFsara UI. On other machines you can get it from: https://rclone.org/install/)
- Token file to authenticate on a specific dCache project allocation. In the next sections we assume that the Data Manager already created the necessary token files for the project members.

## Authentication

### 2. Get the config file that contains your token:

In the examples below we are using a token file provided by a Data Manager that has permissions on a specific project allocation. The token file allows full permissions to any subfolfers under the following `disk` directory: `/pnfs/grid.sara.nl/data/users/anatolid/disk/ada-demo/` 

```sh
# retrieve the ada token
wget -O ada-demo.conf https://surfdrive.surf.nl/files/index.php/s/R6JMdHQ9f3a7saq/download
# the Data Manager created this config file with `get-macaroon --url  https://webdav.grid.surfsara.nl:2880/pnfs/grid.sara.nl/data/users/anatolid/disk/ada-demo/ --duration P7D --chroot --user anatolid --permissions DOWNLOAD,UPLOAD,DELETE,MANAGE,LIST,READ_METADATA,UPDATE_METADATA --output rclone ada-demo`
# inspect the file and find your data path and privileges
cat ada-demo.conf 
view-macaroon ada-demo.conf  # available on any SURFsara UI
```

## Data Transfers

### 3. Download a dataset locally

```sh
rclone --config=ada-demo.conf copy ada-demo:/foo ./ada-demo-folder -P
ls ada-demo-folder/
```

### 4. Upload the dataset and create your own working dir

```sh
rclone --config=ada-demo.conf copy ./ada-demo-folder ada-demo:/<your-name> -P
```

## ADA in practice

### 5. Check if the wrapper is available

```sh
ada --help
```

### 6. Test your token

```sh
ada --tokenfile ada-demo.conf --whoami
#on your laptop you may get this ERROR: no API specified. Use --api <api> or specify a default API in one of the configuration files (/etc/ada.conf /home/<username>/.ada/ada.conf).
#this error means that we need to specify the api address. On Spider we have a default config file in /etc/ada.conf for such settings. On other machines, you need to create this config file in your `~/.ada/ada.conf` home folder with the content here: https://github.com/sara-nl/SpiderScripts/blob/master/ada/etc/ada.conf
#once you create the `~/.ada/ada.conf` file, test again the command: `ada --tokenfile ada-demo.conf --whoami`
```

### 7. List the subdirectories and files

```sh
ada --tokenfile ada-demo.conf --list /<your-name>
ada --tokenfile ada-demo.conf --longlist /<your-name>
```

### 8. Show all details of a file or directory

```sh
ada --tokenfile ada-demo.conf --stat /<your-name>/flowers.jpg
```

### 9. Get checksum

```sh
ada --tokenfile ada-demo.conf --checksum /your-name/ --recursive
# check your debug folder while the cmd is running with: watch "ls -l ~/.ada/headers/"
```

### 10. Create a directory

```sh
ada --tokenfile ada-demo.conf --mkdir /<your-name>/<your-dir>`
```

### 11. Delete a file or directory

```sh
ada --tokenfile ada-demo.conf --delete /<your-name>/flowers.jpg
ada --tokenfile ada-demo.conf --delete /<your-name>/
#WARNING: directory '/<your-name>/' is not empty. If you want to remove it and its contents, you can add the --recursive argument.
ada --tokenfile ada-demo.conf --delete /<your-name>/ --recursive [--force]
```

### 12. Events and stage/unstage operations 

In the examples below we are using a token file provided by a Data Manager that has permissions on a specific project allocation. The token file allows full permissions to any subfolfers under the following `tape` directory: `/pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/` 

```sh
# retrieve the ada token
wget -O ada-demo-tape.conf https://surfdrive.surf.nl/files/index.php/s/2u7zjD3bt7htjxP/download
# the Data Manager created this config file with `get-macaroon --url https://webdav.grid.surfsara.nl:2880/pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/ --duration PT168H --user anatolid --permissions DOWNLOAD,UPLOAD,DELETE,MANAGE,LIST,READ_METADATA,UPDATE_METADATA --output rclone ada-demo-tape`
# inspect the file and find your data path and privileges
cat ada-demo-tape.conf 
view-macaroon ada-demo-tape.conf  # available on any SURFsara UI
```

Let's create some channels to start listening in events.

## 13. Subscribe to any changes in a given directory

```sh
# Here we create a channel to catch any event that happens in this directory
ada --tokenfile ada-demo-tape.conf --longlist /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/
ada --tokenfile ada-demo-tape.conf --channels #available channels
ada --tokenfile ada-demo-tape.conf --events changes-in-folder /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/ --recursive

#open/close dirs in dcache views
#read/write a file
rclone --config=ada-demo-tape.conf sync ./ada-demo-folder ada-demo-tape:/pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/ -P
#files copied to tape display ATTRIB events. This is not telling match on the staging status of the event, so in the next step we will create a channel specifically to track the locality status
```

## 14. Subscribe to all locality and QoS changes in a given directory

```sh 
# Here we create a channel to catch staging events only that happens in this directory
# When you start it, all files in the scope will be listed, including their locality and QoS
ada --tokenfile ada-demo-tape.conf --report-staged changes-in-qos-tape /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/ --recursive

# Delete a dir and recreate it
./ada --tokenfile ada-demo-tape.conf --delete /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/ --recursive
rclone --config=ada-demo-tape.conf sync ./ada-demo-folder ada-demo-tape:/pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/ -P

# Stage/Unstage
#on a different terminal try the following commnads and inspect the changes in the `report-staged` channel
./ada --tokenfile ada-demo-tape.conf --longlist /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/flowers.jpg
./ada --tokenfile ada-demo-tape.conf --stage /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/flowers.jpg
./ada --tokenfile ada-demo-tape.conf --unstage /pnfs/grid.sara.nl/data/users/anatolid/tape/ada-demo-tape/<your-name>/flowers.jpg
```

### 15. Other Authentication options

```sh
# Proxy authentication 
ada --api https://dcacheview.grid.surfsara.nl:22882/api/v1 --proxy --longlist /pnfs/grid.sara.nl/data/lofar/user/sksp/distrib/

# Username/password authentication, via a .netrc file
ada --netrc --longlist /pnfs/grid.sara.nl/data/users/anatolid/
$ cat .netrc
#machine dcacheview.grid.surfsara.nl
#login anatolid
#password XXX
```
