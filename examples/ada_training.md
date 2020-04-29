# ADA demo session

## Preparation

### 1. What you need to run the examples

- Access to a terminal (Spider or any other machine)
- ADA wrapper (will be avail on Spider in production). Get it with:

```sh
# on Spider or any other machine
git clone https://github.com/sara-nl/SpiderScripts.git
```
- Rclone (avail on Spider, on laptop get it from https://rclone.org/install/)

- Token to authenticate on dCache. Get the config file that contains your token:

```sh
# on Spider
cp /tmp/ada-demo.conf .
# on other machine
wget -O ada-demo.conf https://surfdrive.surf.nl/files/index.php/s/3BbhD2eluhOEKGF/download
# the Data Manager created this config file with `get-macaroon --url  https://webdav.grid.surfsara.nl:2880/pnfs/grid.sara.nl/data/users/anatolid/disk/ada-demo/ --duration P7D --chroot --user anatolid --permissions DOWNLOAD,UPLOAD,DELETE,MANAGE,LIST,READ_METADATA,UPDATE_METADATA --output rclone ada-demo`
# inspect the file and find your data path and privileges
cat ada-demo.conf
view-macaroon ada-demo.conf  # on Spider
```

## Data Transfers

### 2. Browse into Ada directory (optional, but makes copy-paste easy)

```sh
cd SpiderScripts/ada/
```sh

### 3. Download a dataset locally

```sh
rclone --config=ada-demo.conf copy ada-demo:/foo ./ada-demo-folder -P
ls ada-demo-folder/
```

### 4. Upload the dataset and create your own working dir

```sh
rclone --config=ada-demo.conf copy ./ada-demo-folder ada-demo:/<your-name> -P
```

Nothing new so far ;)

## ADA in practice

### 5. Check if the wrapper is available

```sh
./ada --help
```

### 6. Test your token

```sh
./ada --tokenfile ada-demo.conf --whoami
#ERROR: no API specified. Use --api <api> or specify a default API in one of the configuration files (/etc/ada.conf /home/<username>/.ada/ada.conf).
```sh

This means that we need to specify the api address. In production we will have a default config file in /etc/ada.conf for such settings.
For now, make your own config file:


```sh
mkdir ~/.ada
nano ~/.ada/ada.conf
# Default settings for the ADA script
   api=https://dcacheview.grid.surfsara.nl:22880/api/v1
   #api=https://dolphin12.grid.surfsara.nl:20443/api/v1
```sh

Test again:

```sh
./ada --tokenfile ada-demo.conf --whoami
```

### 7. List the subdirectories and files

```sh
./ada --tokenfile ada-demo.conf --list /<your-name>
./ada --tokenfile ada-demo.conf --longlist /<your-name>
```

### 8. Show all details of a file or directory

```sh
./ada --tokenfile ada-demo.conf --stat /<your-name>/flowers.jpg
```

### 9. Get checksum

```sh
./ada --tokenfile ada-demo.conf --checksum /your-name/ --recursive
# check your debug folder while the cmd is running with: watch "ls -l ~/.ada/headers/"
```

### 10. Create a directory

```sh
./ada --tokenfile ada-demo.conf --mkdir /<your-name>/<your-dir>`
```

### 11. Delete a file or directory

```sh
./ada --tokenfile ada-demo.conf  --delete /<your-name>/flowers.jpg
./ada --tokenfile ada-demo.conf --delete /<your-name>/
#WARNING: directory '/<your-name>/' is not empty. If you want to remove it and its contents, you can add the --recursive argument.
./ada --tokenfile ada-demo.conf --delete /<your-name>/ --recursive [--force]
```

### 12. Events and stage/unstage demo on Dophin tape pool 

We run this test on our test cluster because the production version lacks behind.

```sh
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --list /users/anatolid/ --api https://dolphin12.grid.surfsara.nl:20443/api/v1
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --channels --api https://dolphin12.grid.surfsara.nl:20443/api/v1
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --events changes-in-tape /users/anatolid/tape --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1

#open/close dirs in dcache views
#read/write a file
rclone --config=tokenfile_dolphin_no_chroot.conf sync ./ada-demo-folder tokenfile_dolphin_no_chroot:/users/anatolid/tape/natalie/ -P
#files copied to tape display ATTRIB events

#  When you start it, all files in the scope will be listed, including their locality and QoS
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --report-staged changes-in-qos-tape /users/anatolid/tape --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1

# Delete dir and recreate it
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --delete /users/anatolid/tape/natalie/ --recursive --api https://dolphin12.grid.surfsara.nl:20443/api/v1
rclone --config=tokenfile_dolphin_no_chroot.conf sync ./ada-demo-folder tokenfile_dolphin_no_chroot:/users/anatolid/tape/natalie/ -P

# Stage/Unstage
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --longlist /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --stage /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
./ada --tokenfile tokenfile_dolphin_no_chroot.conf --unstage /users/anatolid/tape/natalie/flowers.jpg --api https://dolphin12.grid.surfsara.nl:20443/api/v1
```

### 13. Other Authentication

```sh
#proxy
./ada --api https://dcacheview.grid.surfsara.nl:22882/api/v1 --proxy --longlist /pnfs/grid.sara.nl/data/projects.nl/tropomi/natalie/
./ada --netrc --longlist /pnfs/grid.sara.nl/data/users/anatolid/
$ cat .netrc
#machine dcacheview.grid.surfsara.nl
#login anatolid
#password XXX
```
