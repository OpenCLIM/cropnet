FROM continuumio/miniconda3

COPY *.R *.py *.txt ./

RUN conda update --name base conda &&\
    conda create --name s2lai --file docker_conda_env_file_DAenv.txt

# this changes the default shell that is used in any subsequent RUN
# commands, and possibly also any commands we run ourselves in
# a running image built from this container
# Here any shell commands we run will have this prefixed to it
SHELL ["conda", "run", "--name", "s2lai", "/bin/bash", "-c"]

# CMD does not execute anything at build time, 
# but specifies the intended command for the image.
CMD conda run --name s2lai python docker_test.py


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