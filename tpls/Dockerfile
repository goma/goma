FROM ubuntu:22.04

RUN apt update
RUN apt install -y git build-essential m4 zlib1g-dev libx11-dev gfortran pkg-config autoconf python3-dev vim tmux nano gdb valgrind gcc-12 g++-12 gfortran-12 libssl-dev
RUN apt autoremove && apt clean

RUN git clone https://github.com/wortiz/goma -b build-script /opt/goma

RUN PYTHONUNBUFFERED=1 /opt/goma/tpls/install-tpls.py --cc=gcc-12 --cxx=g++-12 --fc=gfortran-12 --build-shared=yes --download-dir=/tmp/downloads -j 32 --extract-dir=/tmp/extract /opt/goma-libs

RUN rm -r /opt/goma

RUN { \
      echo '#!/usr/bin/env bash' \
      && echo ". /opt/goma-libs/config.sh" \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh


ARG USER=goma
RUN adduser -m $USER
RUN groupadd sudo
RUN usermod -a -G sudo $USER
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# make sure everything is in place
USER $USER
ENV HOME /home/$USER
ENV USER $USER
ENV OMPI_MCA_btl_base_warn_component_unuse "0"
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1
WORKDIR $HOME

ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

