# FROM wsmorgan/matdb
FROM tracy_science:flatten


ENV USER_NAME="root"
ENV HOME_DIR="/root"

USER root

# important for matdb
# ENV PATH="/root/.local/bin:${PATH}"
# ENV PATH /root/.local/bin:${PATH}

RUN \
    apt-get update \
    && apt-get -y upgrade \
    && apt-get install -yq python-pip \
    && apt-get install python-tk

RUN apt-get install nano

RUN \
    python -m pip install --upgrade pip \
    && python -m pip install pytest \
    && python -m pip install numpy \
    && python -m pip install tqdm \
    && python -m pip install pandas

RUN mkdir $HOME_DIR/matdb
COPY requirements.txt $HOME_DIR/matdb
RUN python -m pip install --user -r $HOME_DIR/matdb/requirements.txt

COPY matdb $HOME_DIR/matdb/matdb
COPY tests $HOME_DIR/matdb/tests
COPY support $HOME_DIR/matdb/support
COPY setup.py $HOME_DIR/matdb
COPY setup.cfg $HOME_DIR/matdb

RUN \
    cd $HOME_DIR/matdb \
    && python -m pip install --user .
