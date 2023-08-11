# Install base utilities
FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get install -y git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN git clone https://github.com/RosettaCommons/protein_generator.git
WORKDIR /protein_generator
RUN wget http://files.ipd.uw.edu/pub/sequence_diffusion/checkpoints/SEQDIFF_230205_dssp_hotspots_25mask_EQtasks_mod30.pt
RUN wget http://files.ipd.uw.edu/pub/sequence_diffusion/checkpoints/SEQDIFF_221219_equalTASKS_nostrSELFCOND_mod30.pt

RUN conda env create -f environment.yml
RUN pip install poetry

WORKDIR /
RUN git clone https://github.com/broomsday/monomerizer.git
WORKDIR /monomerizer

SHELL ["/bin/bash", "--login", "-c"]
RUN poetry config virtualenvs.create false \
    && poetry install 

ENV PROTEIN_GENERATOR="/protein_generator"
ENV CONDA_SOURCE="/opt/conda/etc/profile.d/conda.sh"

# TODO: remove this once dockerize branch is merged
RUN git fetch --all
RUN git checkout dockerize

ENTRYPOINT ["python", "scripts/make_monomer.py"]