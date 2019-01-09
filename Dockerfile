
FROM mlip4tracy_stable as mtp
FROM pslibrary4tracy as ps
FROM dft_qe4tracy as qe

FROM ubuntu:18.04

RUN apt update && apt install -yq \
    sudo \
    python3.6 \
    python3-pip \
    python3-setuptools \
    autoconf  \
    #ca-certificates  \
    #gawk \
    cmake \
    g++ \
    gdb \
    gfortran-5  \
    make  \
    nano \
    vim \
    wget \
    #--no-install-recommends \
    && apt autoremove   -y \
    && rm -rf /var/apt/lists/* 

ENV DEBIAN_FRONTEND=noninteractive
RUN ln -s /usr/bin/python3 /usr/bin/python

# we create the user 'tracy' and add it to the list of sudoers
RUN adduser -q --disabled-password --gecos tracy tracy          \
    && printf "\ntracy ALL=(ALL:ALL) NOPASSWD:ALL\n" >>/etc/sudoers.d/tracy \
    && (echo "tracy:tracy" | chpasswd)

ENV PS_ROOT "/home/tracy/pslibrary/pz/PSEUDOPOTENTIALS"
ENV MTP_ROOT "/home/tracy/mtp"
ENV QE_ROOT  "/home/tracy/q-e-qe-6.3"
ENV QUIP_ROOT "/home/tracy/quip"
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/Co.* ${PS_ROOT}/
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/W.* ${PS_ROOT}/
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/V.* ${PS_ROOT}/
COPY --from=mtp --chown=tracy:tracy ${MTP_ROOT} ${MTP_ROOT}
COPY --from=qe --chown=tracy:tracy ${QE_ROOT}/PW/src/pw.x ${QE_ROOT}/bin/

# Remove unwanted files
# For non-dev mode, remove source code
ARG DEV_MODE
RUN if [ "${DEV_MODE}" != "YES" ] ; then \
        rm -rf "${QUIP_ROOT}" > /dev/null 2>&1 ; \
    fi

ENV USER_NAME="root"
ENV HOME_DIR="/root"

USER ${USER_NAME}

RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install pytest \
    && python3 -m pip install tqdm \
    && python3 -m pip install pandas

RUN mkdir $HOME_DIR/matdb
COPY requirements.txt $HOME_DIR/matdb
RUN python3 -m pip install -r $HOME_DIR/matdb/requirements.txt

COPY matdb $HOME_DIR/matdb/matdb
COPY tests $HOME_DIR/matdb/tests
COPY support $HOME_DIR/matdb/support
COPY setup.py $HOME_DIR/matdb
COPY setup.cfg $HOME_DIR/matdb
COPY MANIFEST.in $HOME_DIR/matdb
COPY examples $HOME_DIR

RUN cd $HOME_DIR/matdb \
    && python3 -m pip install .

# This is to fix the annoying false error messages when the matdb trying to check for depends.
COPY one_off_fix/pip /usr/bin/pip3

# ASE library assumes each attribute of PP_HEADER element in the ps file is in a separated line
# but pslibrary generates all the attributes in a single line. The below file fixed this issue.
COPY one_off_fix/espresso.py /usr/local/lib/python3.6/dist-packages/ase/io/espresso.py.new

RUN echo export PATH=$PATH:${MTP_ROOT}/bin >> ${HOME_DIR}/.bashrc
WORKDIR "$HOME_DIR"
RUN mkdir -p /root/compute/MTP/CoWV
RUN chmod +x /root/*.sh
