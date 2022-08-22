FROM continuumio/miniconda3

COPY docker_conda_env_file_DAenv.txt ./

RUN conda update --name base conda &&\
    conda create --name s2lai --file docker_conda_env_file_DAenv.txt

COPY Lynch_potpredict_v2_MJB.R docker-wrapper.py MODIS.py utils.py MaxWet1.* UKCP18*.csv ./

# CMD does not execute anything at build time, 
# but specifies the intended command for the image.
CMD conda run --name s2lai python docker-wrapper.py


##################################
#
#
#RUN apt-get update
#RUN apt-get install default-jdk -y
#RUN mkdir -p /data/outputs/
#RUN pip install --no-cache-dir --compile -r requirements.txt
#
# CMD does not execute anything at build time, but specifies the intended command for the image.
#CMD python test.py
#
# this changes the default shell that is used in any subsequent RUN
# commands, and possibly also any commands we run ourselves in
# a running image built from this container
# Here any shell commands we run will have this prefixed to it
# Doesn't seem to be needed
#SHELL ["conda", "run", "--name", "s2lai", "/bin/bash", "-c"]