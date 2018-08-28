FROM wsmorgan/matdb

RUN python -m pip install --upgrade pip
RUN python -m pip install pytest
RUN python -m pip install numpy
RUN python -m pip install tqdm

RUN mkdir /root/matdb
COPY requirements.txt /root/matdb
RUN python -m pip install -r /root/matdb/requirements.txt

COPY matdb /root/matdb/matdb
COPY support /root/matdb/support
COPY setup.py /root/matdb
COPY setup.cfg /root/matdb

WORKDIR /root/matdb
RUN python -m pip install .

COPY . /root/matdb
