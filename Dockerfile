
FROM pslibrary4tracy as ps
FROM dft_qe4tracy as qe
FROM mlip4tracy

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
    && python3 -m pip install pytest-cov \
    && python3 -m pip install tqdm \
    && python3 -m pip install pandas \
    && python3 -m pip install aflow    # this is needed for unit test "tests/calculators/test_aflux.py"

RUN mkdir -p $HOME_DIR/codes/matdb
COPY requirements.txt $HOME_DIR/codes/matdb
RUN python3 -m pip install -r $HOME_DIR/codes/matdb/requirements.txt
RUN ln -s /usr/local/bin/phonopy /usr/bin/phonopy # this is to fix unit test issue for test_vasp.py on local device.

COPY matdb $HOME_DIR/codes/matdb/matdb
COPY tests $HOME_DIR/codes/matdb/tests
COPY support $HOME_DIR/codes/matdb/support
COPY setup.py $HOME_DIR/codes/matdb
COPY setup.cfg $HOME_DIR/codes/matdb
COPY MANIFEST.in $HOME_DIR/codes/matdb
COPY examples $HOME_DIR/codes
COPY scripts/* $HOME_DIR/codes/

RUN cd $HOME_DIR/codes/matdb \
    && python3 -m pip install -e .

# This is to fix the annoying false error messages when the matdb trying to check for depends.
COPY one_off_fix/pip /usr/bin/pip3

# ASE library assumes each attribute of PP_HEADER element in the ps file is in a separated line
# but pslibrary generates all the attributes in a single line. The below file fixed this issue.
COPY one_off_fix/espresso.py /usr/local/lib/python3.5/dist-packages/ase/io/espresso.py.new

RUN echo export PATH=$PATH:${MTP_ROOT}/bin >> ${HOME_DIR}/.bashrc
RUN mkdir -p /root/codes/compute/MTP/CoWV
RUN chmod +x /root/codes/*.sh

RUN sudo mkdir -p -m 0755 /var/run/sshd
CMD ["sudo","/usr/sbin/sshd","-D"]
