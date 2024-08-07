FROM ubuntu:22.04

RUN apt update
RUN apt install -y git build-essential m4 zlib1g-dev libx11-dev gfortran pkg-config autoconf python3-dev vim tmux nano gdb valgrind gcc-12 g++-12 gfortran-12 libssl-dev
RUN apt autoremove && apt clean

RUN git clone https://github.com/wortiz/goma -b build-script /opt/goma

RUN /opt/goma/tpls/install-tpls.py --cc=gcc-12 --cxx=g++-12 --fc=gfortran-12 --build-shared=yes --download-dir=/tmp/downloads -j 32 --extract-dir=/tmp/extract /opt/goma-libs

RUN { \
      echo '#!/usr/bin/env bash' \
      && echo ". /opt/goma-libs/config.sh" \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh

ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

