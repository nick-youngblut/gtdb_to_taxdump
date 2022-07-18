FROM condaforge/mambaforge
RUN apt-get update
RUN apt-get install sudo
COPY environment.yml .
RUN mamba env create -n dev_env -f environment.yml

