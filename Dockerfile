FROM mambaorg/micromamba:1.4.3

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
  maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /foyer

WORKDIR /foyer

RUN apt-get update && apt-get install -y git

RUN micromamba create --file environment-dev.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN  micromamba install -c conda-forge nomkl jupyter python="3.10"
RUN  python setup.py install
RUN  echo "source activate foyer-dev" >> /home/.bashrc
RUN  micromamba clean -afy
RUN  mkdir -p /home/micromamba/data
RUN  chown -R micromamba:micromamba /foyer
RUN  chown -R micromamba:micromamba /opt
RUN  chown -R micromamba:micromamba /home


ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["jupyter"]
