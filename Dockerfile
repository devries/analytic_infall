FROM ubuntu:18.04

RUN apt-get update && apt-get install -y build-essential libgsl0-dev libcfitsio-dev
ADD . /app
WORKDIR /app

RUN ./configure && make && make install

CMD bash
