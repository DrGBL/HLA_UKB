FROM ubuntu

#Sys::Hostname
#/usr/lib/libstdc++.so.6: no version information available

RUN apt-get update && \
	apt-get install build-essential -y && \
	apt-get install wget -y && \
	apt-get install -y perl-doc && \
	apt-get install zip -y && \
	apt-get install libz-dev && \
	apt-get install python3 -y && \
	DEBIAN_FRONTEND=noninteractive  apt-get -y install default-jre-headless  && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* 

WORKDIR /HLA

COPY bowtie2-2.4.4-source.zip /HLA/bowtie2-2.4.4-source.zip
RUN unzip /HLA/bowtie2-2.4.4-source.zip
RUN rm /HLA/bowtie2-2.4.4-source.zip
RUN make -C bowtie2-2.4.4/

ENV PATH="/HLA/bowtie2-2.4.4/:${PATH}"

COPY hlahd.1.4.0.tar.gz /HLA/hlahd.1.4.0.tar.gz
RUN tar -zxf /HLA/hlahd.1.4.0.tar.gz
RUN rm /HLA/hlahd.1.4.0.tar.gz
RUN cd /HLA/hlahd.1.4.0/ && \
	sh install.sh && \
	sh update.dictionary.sh

ENV PATH="/HLA/hlahd.1.4.0/bin:${PATH}"

