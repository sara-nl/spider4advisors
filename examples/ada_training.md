# ADA demo session

## Preparation

### 1. What you need to run the examples

- Access to a terminal (Spider UI or any other machine)
- ADA wrapper (is available on Spider. On other machines you need to download it from: `git clone https://github.com/sara-nl/SpiderScripts.git`)
- Rclone (available on any SURFsara UI, On other machines you can get it from: https://rclone.org/install/)
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
```

This means that we need to specify the api address. On Spider we have a default config file in /etc/ada.conf for such settings.
On other machines, use this config file in `~/.ada/ada.conf`: https://github.com/sara-nl/SpiderScripts/blob/master/ada/etc/ada.conf

Test again:

```sh
ada --tokenfile ada-demo.conf --whoami
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

### 12. Events and stage/unstage demo 




```sh
ada --tokenfile tokenfile_dolphin_no_chroot.conf --list /users/anatolid/ --api https://dolphin12.grid.surfsara.nl:20443/api/v1
ada --tokenfile tokenfile_dolphin_no_chroot.conf --channels --api https://dolphin12.grid.surfsara.nl:20443/api/v1
ada --tokenfile tokenfile_dolphin_no_chroot.conf --events changes-in-tape /users/anatolid/tape --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1

#open/close dirs in dcache views
#read/write a file
rclone --config=tokenfile_dolphin_no_chroot.conf sync ./ada-demo-folder tokenfile_dolphin_no_chroot:/users/anatolid/tape/natalie/ -P
#files copied to tape display ATTRIB events

#  When you start it, all files in the scope will be listed, including their locality and QoS
ada --tokenfile tokenfile_dolphin_no_chroot.conf --report-staged changes-in-qos-tape /users/anatolid/tape --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1

# Delete dir and recreate it
ada --tokenfile tokenfile_dolphin_no_chroot.conf --delete /users/anatolid/tape/natalie/ --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1
rclone --config=tokenfile_dolphin_no_chroot.conf sync ./ada-demo-folder tokenfile_dolphin_no_chroot:/users/anatolid/tape/natalie/ -P

# Stage/Unstage
ada --tokenfile tokenfile_dolphin_no_chroot.conf --longlist /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
ada --tokenfile tokenfile_dolphin_no_chroot.conf --stage /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
ada --tokenfile tokenfile_dolphin_no_chroot.conf --unstage /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
```

### 13. Other Authentication

```sh
# Proxy
ada --api https://dcacheview.grid.surfsara.nl:22882/api/v1 --proxy --longlist /pnfs/grid.sara.nl/data/lofar/user/sksp/distrib/

# Username/password, via a .netrc file
ada --netrc --longlist /pnfs/grid.sara.nl/data/users/anatolid/
$ cat .netrc
#machine dcacheview.grid.surfsara.nl
#login anatolid
#password XXX
```
