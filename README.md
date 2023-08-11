# Monomerizer
Take a PDB structure of a symmetric homo-oligomeric barrel, and monomerize without occluding the pore.

## Building the Container Image
sudo docker build -t monomerizer:latest .

## Running the Container
sudo docker run -it -v $(pwd)/tmp_docker:/data --gpus all monomerizer:latest "/data/inputs/2XQS_short.pdb" "/data/outputs" "--generator-steps" "2" "--num-designs" "1" "--linker-lengths" "11"
