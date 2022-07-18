FROM continuumio/miniconda3
COPY environment.yml .
RUN conda install -c conda-forge mamba
RUN /bin/bash -c "mamba env create -f environment.yml"

