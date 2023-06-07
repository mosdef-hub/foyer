ARG PY_VERSION=3.10
FROM continuumio/miniconda3:4.10.3-alpine AS builder

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
  maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /foyer

WORKDIR /foyer
ARG PY_VERSION=3.10

# Create a group and user
RUN addgroup -S anaconda && adduser -S anaconda -G anaconda

RUN /sbin/apk add --no-cache git
RUN conda update conda -yq && conda install -c conda-forge mamba
RUN conda config --set always_yes yes --set changeps1 no
RUN . /opt/conda/etc/profile.d/conda.sh
RUN  sed -i -E "s/python.*$/python="$(PY_VERSION)"/" environment-dev.yml
RUN  mamba env create nomkl -f environment-dev.yml
RUN  conda activate foyer-dev
RUN  mamba install -c conda-forge nomkl jupyter python="$PY_VERSION"
RUN  python setup.py install
RUN  echo "source activate foyer-dev" >> /home/anaconda/.profile
RUN  conda clean -afy
RUN  mkdir -p /home/anaconda/data
RUN  chown -R anaconda:anaconda /foyer
RUN  chown -R anaconda:anaconda /opt
RUN  chown -R anaconda:anaconda /home/anaconda

WORKDIR /home/anaconda

COPY devtools/docker-entrypoint.sh /entrypoint.sh

RUN chmod a+x /entrypoint.sh

USER anaconda

ENTRYPOINT ["/entrypoint.sh"]
CMD ["jupyter"]
