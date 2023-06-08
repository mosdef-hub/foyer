FROM mambaorg/micromamba:1.4.3

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
  maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /foyer

WORKDIR /foyer

RUN micromamba install --yes --file environment-dev.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN  micromamba activate foyer-dev
RUN  micromamba install -c conda-forge nomkl jupyter python="3.10"
RUN  python setup.py install
RUN  echo "source activate foyer-dev" >> /home/anaconda/.profile
RUN  micromamba clean -afy
RUN  mkdir -p /home/anaconda/data
RUN  chown -R anaconda:anaconda /foyer
RUN  chown -R anaconda:anaconda /opt
RUN  chown -R anaconda:anaconda /home/anaconda


ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["jupyter"]
