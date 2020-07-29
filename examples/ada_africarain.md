# Africarain introduction to ADA

## dCache Overview

dCache is a powerful data storage platform used for many data intensive projects, that can be interfaced with in a number of ways. Our ADA (Advanced dCache API) interface is based on the dCache API and the webdav protocol to access and process your data on dCache from any platform and with various authentication methods.

## Use cases

### Generic workflow on Spider

![Spider_ADA_example_simple](https://user-images.githubusercontent.com/12894031/86337875-f5c62d80-bc51-11ea-87cd-142881ea3289.png)

### Automated workflow on Spider

![Spider_ADA_example_automated](https://user-images.githubusercontent.com/12894031/86337814-e2b35d80-bc51-11ea-865d-7da677bf7987.png)

## Authentication

The Data manager has direct access to dCache and is able to use the Spider credentials to view the data and create tokens for other project members. 
The project members can work with dCache as long as they have access to a token.

### Access data (Data manager)

dCache storage can be viewed both through the Ada tools or through your browser using the web client, this is just one additional way you can explore the storage space. You can log-in using you Spider credentials [africarain-username]

NOTE: you may be asked for a browser certificate, just select cancel and you will be asked for your credentials

`https://webdav.grid.surfsara.nl:2880/pnfs/grid.sara.nl/data/africarain/disk/ `

### Create a macaroon (Data manager)

The project members need a token file to authenticate on a specific dCache project allocation. These tokens are called macaroons and can be created by the Data Manager.
Macaroons can be used to give access to dCache data in a very granular way. This enables data managers autonomously share their data in dCache without having to reach out to SURFsara to request access. 

Obtain a macaroon with the following command:

``` /bash
get-macaroon \ 
    --url [project storage url] \
    --duration [validity time] \ 
    --chroot \ 
    --user [africarain-username] \ 
    --permissions [allowed activities] \
    --output rclone [<name>.conf tokenfile]
```

Create a tokenfile on the Spider UI:

``` /bash
get-macaroon --url https://webdav.grid.surfsara.nl:2880/pnfs/grid.sara.nl/data/africarain/disk/ --duration P7D --chroot  --user africarain-clecoz --permissions DOWNLOAD,UPLOAD,DELETE,MANAGE,LIST,READ_METADATA,UPDATE_METADATA --output rclone africarain
```

### Share macaroons (Data manager)

The config file generated in the step above can be shared with project members and collaborators for them to access their data. The holder of this config file can operate on the dCache project data directly and thus, the config file should be shared with the project team in a non-public space, for example user's home directories, or the 'Shared' or 'Data' project space directories on Spider.

`cp africarain.conf /project/africarain/Data/`


### Access data (Non authenticated project members)

dCache storage can be viewed both through the Ada tools or through your browser using the bearer token inside the tokenfile. Paste in the browser the following url:

`https://webdav.grid.surfsara.nl:2880/?authz=[bearer_token]`


## rClone / webdav

**rclone** is a webdav client that supports by default 4 parallel streams of data, and is installed on the Spider platform. 

### Write data to dCache

`rclone --config=africarain.conf copy ./[SOURCE]/ africarain:[DESTINATION] -P`

Example:

`rclone --config=africarain.conf copy ./ada-demo-folder/ africarain:/ada-demo-remote/ -P`


### Copy data from dCache

`rclone --config=africarain.conf copy africarain:/[SOURCE] ./[DESTINATION] -P`

Example, copy an existing test folder to Spider

`rclone --config=africarain.conf copy africarain:/ada-demo-remote/ ./ada-demo-local/ -P`


## Ada

ADA is a wrapper of tools created by SURFsara to simplify your interactions with dCache. Rclone can support uploading and downloading data but other operations such as listing or deleting files and directories can be performed directly on the dCache API. ADA wraps all of this functionality into one clean package saving you the hassle of having to download and troubleshoot multiple packages and dependencies. ADA is installed on Spider

### Check your access to the system

`ada --tokenfile africarain.conf --whoami`
``` json
{
  "status": "AUTHENTICATED",
  "uid": 51539,
  "gids": [
    51181
  ],
  "username": "africarain-clecoz",
  "rootDirectory": "/pnfs/grid.sara.nl/data/africarain/disk",
  "homeDirectory": "/"
}
```

### Create a directory on dCache

`ada --tokenfile africarain.conf --mkdir [DIRECTORY]`

Example, create a data structure:

`ada --tokenfile africarain.conf --mkdir output-data`

### Move files on dCache

`ada --tokenfile africarain.conf --mv [SOURCE] [DESTINATION]`

Example, move existing remote data to a different folder:

`ada --tokenfile africarain.conf --mv /ada-demo-remote /output-data/ada-demo-remote`

### List files on dCache

`/ada --tokenfile africarain.conf --longlist [DIRECTORY]`

Example, list all properties in a folder:

`ada --tokenfile africarain.conf --longlist /output-data/ada-demo-remote`

### Recursively remove folders 

`ada --tokenfile africarain.conf --delete [DIRECTORY] --recursive --force`

Example, remove the created folders:

`ada --tokenfile africarain.conf --delete /output-data --recursive --force`





