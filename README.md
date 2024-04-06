# Monomerizer
Take a PDB structure of a symmetric homo-oligomeric barrel, and monomerize without occluding the pore.

Monomerization of shortened version of PDB 2XQS is shown with only small loops added to monomerize.
![image](images/2XQS_431_combined.png)

Larger modifications can be made aswell such as the large helical linkers for shortened version of PDB 4JPP.
![image](images/4JPP_640_combined.png)

## Running Locally as a Single Command
Running requires only a PDB file as input.  In the example case we use a shortened version of the `2XQS` structure
from the RCSB.  This file can be found in the repository as `tests/pdbs/2XQS_short.pdb`.

### Make a working directory and move over the input file
1. `mkdir docker_tmp`
2. `cp monomerizer/tests/pdbs/2XQS_short.pdb ./docker_tmp`

### Make sure the docker deamon is started
1. `sudo systemctl start docker`

### Interactively start the docker image
1. `sudo docker run -it --entrypoint /bin/bash -v $(pwd)/tmp_docker/:/workspace --gpus all antiquatedarachnid/monomerizer:latest`

### Run `monomerizer`
1. `python ../monomerizer/scripts/make_monomer.py 2XQS_short.pdb output --generator-steps 2 --num-designs 1 --linker-lengths 11`

## Running on the Cloud
In this example we'll run using `runpod` (https://www.runpod.io/) which will allow easy access to larger and or faster
GPUs than you might have locally.

### Setup a template for monomerizer
1. Make a new runpod template with the Container Image as `antiquatedarachnid/monomerizer:latest`
2. Set the Docker Command to `/bin/bash`

### Start the container
1. Launch the container using the `runpod` WebGUI
2. ssh into the pod using terminal or web-connect options (see RunPod Web GUI buttons)

### Copy over your desired input file
1. `runpodctl send myfile.pdb` on e.g. your local computer and get the `code`
2. `runpodctl receive {code}` on the pod

### Run `monomerizer`
1. `python ../monomerizer/scripts/make_monomer.py 2XQS_short.pdb output --generator-steps 2 --num-designs 1 --linker-lengths 11`

### Compress the results and send to your local machine
1. `tar -czvf outputs.tar.gz outputs`
2. `runpodctl send outputs.tar.gz`
3. `runpotctl receive {code}`

## Development
If you make changes to the code, you can locally use a the new docker image or publish your own

### Building the Container Image
`sudo docker build -t monomerizer:latest .`
- add `--no-cache` if needing a full rebuild

### Publishing the Docker Image
```
sudo docker login
sudo docker tag monomerizer:latest your_docker_id/monomerizer:latest
sudo docker push your_docker_id/monomerizer:latest
```
