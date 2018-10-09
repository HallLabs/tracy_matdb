
FROM quip4tracy as quip
FROM mtp4tracy as mtp  

FROM pslibrary4tracy

ENV QUIP_ROOT "/home/tracy/quip"
ENV MTP_ROOT "/home/tracy/mtp"
COPY --from=quip --chown=tracy:tracy ${QUIP_ROOT} ${QUIP_ROOT}
COPY --from=mtp --chown=tracy:tracy ${MTP_ROOT} ${MTP_ROOT}

ENV USER_NAME="root"
ENV HOME_DIR="/root"

USER ${USER_NAME}

RUN apt-get update \
    && apt-get -yq upgrade \
    && apt-get install -yq python-tk gdb

RUN python -m pip install pytest \
    && python -m pip install tqdm \
    && python -m pip install pandas

RUN mkdir $HOME_DIR/matdb
COPY requirements.txt $HOME_DIR/matdb
RUN python -m pip install -r $HOME_DIR/matdb/requirements.txt

COPY matdb $HOME_DIR/matdb/matdb
COPY tests $HOME_DIR/matdb/tests
COPY support $HOME_DIR/matdb/support
COPY setup.py $HOME_DIR/matdb
COPY setup.cfg $HOME_DIR/matdb
COPY MANIFEST.in $HOME_DIR/matdb
COPY examples/Qe_CoWV.yml $HOME_DIR

RUN cd $HOME_DIR/matdb \
    && python -m pip install .

ENV QUIP_ARCH linux_x86_64_gfortran
RUN cd ${QUIP_ROOT} \
    && make install-quippy > /dev/null 

# This is to fix the annoying false error messages when the matdb trying to check for depends.
COPY one_off_fix/pip /usr/bin/pip

# ASE library assumes each attribute of PP_HEADER element in the ps file is in a separated line
# but pslibrary generates all the attributes in a single line. The below file fixed this issue.
COPY one_off_fix/espresso.py /usr/local/lib/python2.7/dist-packages/ase/io/espresso.py

RUN echo export PATH=$PATH:${MTP_ROOT}/bin >> ${HOME_DIR}/.bashrc
WORKDIR "$HOME_DIR"
