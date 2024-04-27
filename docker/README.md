## Build the docker image for risefl

Steps for building the risefl docker image from `risefl.Dockerfile`

```bash
docker build --network=host -t risefl:latest -f ./risefl.Dockerfile .
```

Useful docker build flags:
- `-f`, file flag, for specifying the location of the dockerfile (By default, the Dockerfile is assumed to be located here, but you can specify a different location with the file flag (-f))
- `-t`, tag flag, for specifying a tag for the image `name:tag` format

Run the docker container and test the environment
```bash
docker run -it risefl:latest /bin/bash
```
