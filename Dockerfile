FROM continuumio/miniconda3:latest

RUN mkdir /opt/taranis/
ADD utils /opt/taranis/utils
ADD test /opt/taranis/test
ADD *.py /opt/taranis/
ADD environment.yml /opt/taranis/
ADD README.md /opt/taranis/
ADD LICENSE /opt/taranis/

RUN cd /opt/taranis
RUN /opt/conda/bin/conda env create -f /opt/taranis/environment.yml && /opt/conda/bin/conda clean -a
RUN /opt/conda/bin/conda env export --name taranis > taranis.yml
ENV PATH /opt/conda/envs/taranis:/opt/conda/envs/taranis/utils:$PATH
