
FROM mtp4tracy as mtp  
FROM pslibrary4tracy as ps
FROM dft_qe4tracy as qe

FROM quip4tracy

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
COPY examples $HOME_DIR

RUN cd $HOME_DIR/matdb \
    && python -m pip install .

# This is to fix the annoying false error messages when the matdb trying to check for depends.
COPY one_off_fix/pip /usr/bin/pip

# ASE library assumes each attribute of PP_HEADER element in the ps file is in a separated line
# but pslibrary generates all the attributes in a single line. The below file fixed this issue.
COPY one_off_fix/espresso.py /usr/local/lib/python2.7/dist-packages/ase/io/espresso.py

RUN echo export PATH=$PATH:${MTP_ROOT}/bin >> ${HOME_DIR}/.bashrc
WORKDIR "$HOME_DIR"
RUN mkdir -p /root/compute/MTP/CoWV
RUN chmod +x /root/*.sh
