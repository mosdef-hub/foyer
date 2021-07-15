Using foyer with Docker
========================

As much of scientific software development happens in unix platforms, to avoid the quirks of development dependent on system you use, a recommended way is to use docker or other containerization technologies. This section is a how to guide on using foyer with docker.

Prerequisites
-------------
A docker installation in your machine. Follow this `link <https://docs.docker.com/get-docker/>`_ to get a docker installation working on your machine. If you are not familiar with docker and want to get started with docker, the Internet is full of good tutorials like the ones `here <https://docker-curriculum.com/>`_ and `here <https://www.youtube.com/watch?v=zJ6WbK9zFpI&feature=youtu.be>`_.

Quick Start
-----------
After you have a working docker installation, please use the following command to use run a jupyter-notebook with all the dependencies for `foyer` installed:

.. code-block:: bash

    $ docker pull mosdef/foyer:latest
    $ docker run -it --name foyer -p 8888:8888 mosdef/foyer:latest su anaconda -s\
      /bin/sh -l -c "jupyter-notebook --no-browser --ip="0.0.0.0" --notebook-dir\
      /home/anaconda/foyer-notebooks"


If every thing happens correctly, you should a be able to start a `jupyter-notebook` server running in a python environment with all the dependencies for `foyer` installed.

Alternatively, you can also start a Bourne shell to use python from the container's terminal:

.. code-block:: bash

    $ docker run -it --name foyer mosdef/foyer:latest

.. important::

    The instructions above will start a docker container but containers by nature are ephemeral, so any filesystem changes (like adding a new notebook) you make will only persist till the end of the container's lifecycle. If the container is removed, any changes or code additions will not persist.

Persisting User Volumes
-----------------------
If you will be using `foyer` from a docker container, a recommended way is to mount what are called user volumes in the container. User volumes will provide a way to persist all filesystem/code additions made to a container regardless of the container lifecycle. For example, you might want to create a directory called `foyer-notebooks` in your local system, which will store all your `foyer` notebooks/code. In order to make that accessible to the container(where the notebooks will be created/edited), use the following steps:


1. Create a directory in your filesystem

.. code-block:: bash

    $ mkdir -p /path/to/foyer-notebooks
    $ cd /path/to/foyer-notebooks

2. Define an entry-point script. Inside `foyer-notebooks` in your local file system create a file called :code:`dir_entrypoint.sh` and paste the following content.

.. code-block:: bash

    #!/bin/sh

    chown -R anaconda:anaconda /home/anaconda/foyer-notebooks

    su anaconda -s /bin/sh -l -c "jupyter-notebook --no-browser --ip="0.0.0.0" --notebook-dir /home/anaconda/foyer-notebooks"

3. Run docker image for `foyer`

.. code-block:: bash

    $ docker run -it --name foyer -p 8888:8888 --entrypoint /home/anaconda/foyer-notebooks/dir_entrypoint.sh -v $HOME/foyer-notebooks:/home/anaconda/foyer-notebooks mosdef/foyer:latest


Cleaning Up
-----------
You can remove the created container by using the following command:

.. code-block:: bash

    $ docker container rm foyer

.. note::

    Instead of using `latest`, you can use the image :code:`mosdef/foyer:stable` for most recent stable release of `foyer` and run the tutorials.
