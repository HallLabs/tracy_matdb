
FROM pslibrary4tracy as ps
FROM dft_qe4tracy as qe
FROM mlip4tracy_stable

ENV PS_ROOT "/home/tracy/pslibrary/pz/PSEUDOPOTENTIALS"
ENV MTP_ROOT "/home/tracy/mtp"
ENV QE_ROOT  "/home/tracy/q-e-qe-6.3"
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/Co.* ${PS_ROOT}/
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/W.* ${PS_ROOT}/
COPY --from=ps --chown=tracy:tracy ${PS_ROOT}/V.* ${PS_ROOT}/
COPY --from=qe --chown=tracy:tracy ${QE_ROOT}/PW/src/pw.x ${QE_ROOT}/bin/

ENV USER_NAME="root"
ENV HOME_DIR="/root"

USER ${USER_NAME}
WORKDIR "$HOME_DIR"

RUN apt update && apt install -yq \
    python3 \
    python3-pip \
    python3-setuptools \
    && apt autoremove   -y \
    && rm -rf /var/apt/lists/* 

RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install pytest \
    && python3 -m pip install tqdm \
    && python3 -m pip install pandas \
    && python3 -m pip install aflow    # this is needed for unit test "tests/calculators/test_aflux.py"

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
COPY one_off_fix/espresso.py /usr/local/lib/python3.5/dist-packages/ase/io/espresso.py.new

RUN echo export PATH=$PATH:${MTP_ROOT}/bin >> ${HOME_DIR}/.bashrc
RUN mkdir -p /root/compute/MTP/CoWV
RUN chmod +x /root/*.sh
